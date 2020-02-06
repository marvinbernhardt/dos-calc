#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdbool.h>
#include <argp.h>
#include <cblas.h>
#include <math.h>
#include <chemfiles.h>
#include "velocity-decomposition.c"
#include "fft.c"
#include "verbPrintf.c"
#include "parse-dosparams.c"
#include "write-dos.c"
#include "structs.h"

// CLI stuff
const char* argp_program_version = VERSION;
const char* argp_program_bug_address = "<bernhardt@cpc.tu-darmstadt.de>";
static char doc[] = "dos-calc -- a programm to calculate densities of states from trajectories";
static char args_doc[] = "dosparams trajectory";

static struct argp_option options[] = {
    {"verbose",         'v', 0,      0, "Produce verbose output", 0},
    {"no-pbc",          'p', 0,      0, "Do not recombine molecules seperated by periodic boundary conditions", 0},
    {"outfile",         'o', "FILE", 0, "Output json file. Default: dos.json", 0},
    {"framelength",     'f', "FL",   0, "Lenght of each frame in trajectory (in ps). Taken from trajectory if possible, but can be overwritten with this argument.", 0},
    { 0 }
};

struct arguments
{
    char *(input_files[2]);
    bool verbosity;
    bool no_pbc;
    char *outfile;
    float framelength;
};

static error_t parse_opt (int key, char *arg, struct argp_state *state)
{
    /* Get the input argument from argp_parse, which we
       know is a pointer to our arguments structure. */
    struct arguments *arguments = state->input;

    switch (key)
    {
        case 'v':
            arguments->verbosity = true;
            break;
        case 'p':
            arguments->no_pbc = true;
            break;
        case 'o':
            arguments->outfile = arg;
            break;
        case 'f':
            arguments->framelength = strtof(arg, NULL);
            break;

        case ARGP_KEY_ARG:
            /* Too many arguments. */
            if (state->arg_num > 2)
            {
                argp_usage(state);
            }
            arguments->input_files[state->arg_num] = arg;
            break;

        case ARGP_KEY_END:
            if (state->arg_num < 2)
                /* Not enough arguments. */
                argp_usage (state);
            break;

        default:
            return ARGP_ERR_UNKNOWN;
    }
    return 0;
}


/* Our argp parser. */
static struct argp argp = {options, parse_opt, args_doc, doc, 0, 0, 0};

int main( int argc, char *argv[] )
{
    // command line arguments
    struct arguments arguments;

    // Default values.
    arguments.verbosity = false;
    arguments.no_pbc = false;
    arguments.outfile = "dos.json";
    arguments.framelength = 0.0;

    // parse command line arguments
    argp_parse(&argp, argc, argv, 0, 0, &arguments);

    // define convenience aliases
    bool verbosity = arguments.verbosity;
    char *dosparams_file = arguments.input_files[0];
    char *trajectory_file = arguments.input_files[1];
    verbPrintf(arguments.verbosity, "dosparams file: %s\n", dosparams_file);
    verbPrintf(arguments.verbosity, "trajectory file: %s\n", trajectory_file);

    // input that will be scanned
    size_t nsamples;
    size_t nblocks;
    unsigned long nblocksteps;
    size_t nmoltypes;
    size_t *moltypes_nmols;
    size_t *moltypes_natomspermol;
    float **moltypes_atommasses;
    char *moltypes_rot_treat;
    int **moltypes_abc_indicators;
    size_t ncross_spectra;
    cross_spectrum_def *cross_spectra_def;

    // parse numbers from json
    // NOTE: this calls alloc() functions, use free_dosparams_arrays() in the end
    parse_dosparams(dosparams_file,
            &nsamples,  // output
            &nblocks,
            &nblocksteps,
            &nmoltypes,
            &moltypes_nmols,
            &moltypes_natomspermol,
            &moltypes_atommasses,
            &moltypes_rot_treat,
            &moltypes_abc_indicators,
            &ncross_spectra,
            &cross_spectra_def);
    verbPrintf(verbosity, "nsamples: %zu\n", nsamples);
    verbPrintf(verbosity, "nblocks: %zu\n", nblocks);
    verbPrintf(verbosity, "nblocksteps: %zu\n", nblocksteps);
    verbPrintf(verbosity, "nmoltypes: %zu\n", nmoltypes);
    for (size_t h=0; h<nmoltypes; h++)
    {
        verbPrintf(verbosity, "moltype %zu nmols: %zu\n", h, moltypes_nmols[h]);
        verbPrintf(verbosity, "moltype %zu natomspermol: %zu\n", h, moltypes_natomspermol[h]);
        verbPrintf(verbosity, "moltype %zu atommasses: ", h);
        for (size_t j=0; j<moltypes_natomspermol[h]; j++)
            verbPrintf(verbosity, "%f ", moltypes_atommasses[h][j]);
        verbPrintf(verbosity, "\n");
        verbPrintf(verbosity, "moltype %zu rot_treat: %c\n", h, moltypes_rot_treat[h]);
        verbPrintf(verbosity, "moltype %zu abc_indicators: ", h);
        for (size_t j=0; j<4; j++)
            verbPrintf(verbosity, "%d ", moltypes_abc_indicators[h][j]);
        verbPrintf(verbosity, "\n");
    }

    // generate convenience variables and arrays
    // convenience variables
    size_t natoms = 0;
    size_t nmols = 0;
    unsigned long nfrequencies = nblocksteps / 2 + 1;
    size_t *moltypes_firstmol = calloc(nmoltypes, sizeof(size_t));
    size_t *moltypes_firstatom = calloc(nmoltypes, sizeof(size_t));
    for (size_t h=0; h<nmoltypes; h++)
    {
        moltypes_firstmol[h] = nmols;
        moltypes_firstatom[h] = natoms;
        nmols += moltypes_nmols[h];
        natoms += moltypes_nmols[h] * moltypes_natomspermol[h];
    }
    size_t *mols_moltypenr = calloc(nmols, sizeof(size_t));
    size_t *mols_natoms = calloc(nmols, sizeof(size_t));
    float *mols_mass = calloc(nmols, sizeof(float));
    for (size_t i=0; i<nmols; i++)
    {
        for (size_t h=0; h<nmoltypes; h++)
        {
            if (i >= moltypes_firstmol[h] && i < moltypes_firstmol[h] + moltypes_nmols[h])
                mols_moltypenr[i] = h;
        }
        mols_natoms[i] = moltypes_natomspermol[mols_moltypenr[i]];
        mols_mass[i] = 0;
        for (size_t j=0; j<mols_natoms[i]; j++)
            mols_mass[i] += moltypes_atommasses[mols_moltypenr[i]][j];
    }
    size_t *mols_firstatom = calloc(nmols, sizeof(size_t));
    for (size_t h=0; h<nmoltypes; h++)
    {
        for (size_t ii=0; ii<moltypes_nmols[h]; ii++)
        {
            mols_firstatom[moltypes_firstmol[h]+ii] = moltypes_firstatom[h] + moltypes_natomspermol[h] * ii;
        }
    }

    // open file once for tests
    verbPrintf(verbosity, "starting with file %s\n", trajectory_file);
    CHFL_TRAJECTORY* file = chfl_trajectory_open(trajectory_file, 'r');
    CHFL_FRAME* frame = chfl_frame();
    int result = chfl_trajectory_read(file, frame);

    if ((file == NULL) || (result != CHFL_SUCCESS)) {
        fprintf(stderr, "ERROR: Reading trajectory failed.\n");
        return 1;
    }

    // check for number of atoms
    uint64_t natoms_traj = 0;
    chfl_vector3d* positions = NULL;
    chfl_frame_positions(frame, &positions, &natoms_traj);
    if (natoms_traj < natoms) {
        fprintf(stderr, "ERROR: The topology you give has more atoms than first frame of the trajectory\n");
        return 1;
    }
    else if (natoms_traj > natoms) {
        fprintf(stderr, "WARNING: The topology you give has less atoms than first frame of the trajectory\n");
        fprintf(stderr, "         Some atoms are ignored in every frame\n");
    }

    // check for velocities
    bool has_velocities = false;
    chfl_frame_has_velocities(frame, &has_velocities);
    if (!has_velocities) {
        fprintf(stderr, "ERROR: No velocities in trajectory.\n");
        return 1;
    }

    // check for orthorombic box
    CHFL_CELL* cell = chfl_cell_from_frame(frame);
    chfl_cellshape shape;
    chfl_cell_shape(cell, &shape);
    if((arguments.no_pbc == false) && (shape != CHFL_CELL_ORTHORHOMBIC)) {
        fprintf(stderr, "ERROR: can not do recombination on non orthorombic box.\n");
        return 1;
    }
    chfl_free(cell);

    // get framelength
    float framelength;
    if (arguments.framelength == 0.0) {
        double time0, time1;
        int result0, result1;
        // get time0
        CHFL_PROPERTY *property = chfl_frame_get_property(frame, "time");
        result0 = chfl_property_get_double(property, &time0);
        verbPrintf(verbosity, "time of first frame is %f ps\n", time0);
        chfl_free(property);
        // get time1
        chfl_trajectory_read(file, frame);
        property = chfl_frame_get_property(frame, "time");
        result1 = chfl_property_get_double(property, &time1);
        verbPrintf(verbosity, "time of second frame is %f ps\n", time1);
        chfl_free(property);
        framelength = (float) (time1 - time0);
        if ((framelength == 0.0)
            || (result0 != CHFL_SUCCESS)
            || (result1 != CHFL_SUCCESS)) {
            fprintf(stderr, "ERROR: Reading framelength from trajectory failed. You can provide the framelength with command line arguments.\n");
            return 1;
        }
    }
    else {
        framelength = arguments.framelength;
    }
    chfl_free(frame);
    chfl_trajectory_close(file);
    verbPrintf(verbosity, "framelength is %f ps\n", framelength);

    // open file for calculations
    file = chfl_trajectory_open(trajectory_file, 'r');

    // output arrays
    const size_t ndos = 12;
    const char *dos_names[12] =
    { "trn_x",
      "trn_y",
      "trn_z",
      "rot_x",
      "rot_y",
      "rot_z",
      "vib_x",
      "vib_y",
      "vib_z",
      "roto_a",
      "roto_b",
      "roto_c"
    };
    // order is: trn_xyz, rot_xyz, vib_xyz, rot_omega_abc
    float* moltypes_dos_samples = calloc(nmoltypes*ndos*nsamples*nfrequencies, sizeof(float));
    float* cross_spectra_samples = calloc(ncross_spectra*nsamples*nfrequencies, sizeof(float));
    // moments of inertia
    float* moltypes_samples_moments_of_inertia = calloc(nmoltypes*nsamples*3, sizeof(float));

    // start samples loop
    verbPrintf(verbosity, "going through %zu samples\n", nsamples);
    for (size_t sample=0; sample<nsamples; sample++)
    {
        verbPrintf(verbosity, "now doing sample %zu\n", sample);

        float* mol_moments_of_inertia = calloc(nmols*3, sizeof(float));

        // start block loop
        verbPrintf(verbosity, "going through %zu blocks\n", nblocks);
        for (size_t block=0; block<nblocks; block++)
        {
            verbPrintf(verbosity, "now doing block %zu\n", block);
            verbPrintf(verbosity, "start decomposition\n");

            float* mol_velocities_sqrt_m_trn = calloc(nmols*3*nblocksteps, sizeof(float));
            float* mol_omegas_sqrt_i_rot = calloc(nmols*3*nblocksteps, sizeof(float));
            float* atom_velocities_sqrt_m_vib = calloc(natoms*3*nblocksteps, sizeof(float));
            float* atom_velocities_sqrt_m_rot = calloc(natoms*3*nblocksteps, sizeof(float));

            decomposeVelocities (file,
                    nblocksteps,
                    nmols,
                    nmoltypes,
                    mols_firstatom,
                    mols_natoms,
                    mols_moltypenr,
                    moltypes_atommasses,
                    mols_mass,
                    moltypes_natomspermol,
                    moltypes_rot_treat,
                    moltypes_abc_indicators,
                    arguments.no_pbc,
                    mol_velocities_sqrt_m_trn,
                    mol_omegas_sqrt_i_rot,
                    atom_velocities_sqrt_m_vib,
                    atom_velocities_sqrt_m_rot,
                    mol_moments_of_inertia);


            verbPrintf(verbosity, "start DoS calculation (FFT)\n");
            dos_calculation(nmoltypes,
                    nblocksteps,
                    nfrequencies,
                    moltypes_firstmol,
                    moltypes_firstatom,
                    moltypes_nmols,
                    moltypes_natomspermol,
                    mol_velocities_sqrt_m_trn,
                    mol_omegas_sqrt_i_rot,
                    atom_velocities_sqrt_m_vib,
                    atom_velocities_sqrt_m_rot,
                    ndos,
                    nsamples,
                    sample,
                    ncross_spectra,
                    cross_spectra_def,
                    moltypes_dos_samples, //output
                    cross_spectra_samples
                    );

            //free velocities
            free(mol_velocities_sqrt_m_trn);
            free(mol_omegas_sqrt_i_rot);
            free(atom_velocities_sqrt_m_vib);
            free(atom_velocities_sqrt_m_rot);
        }
        verbPrintf(verbosity, "finished all blocks\n");

        // divide moi by number of blocks and number of blocksteps
        cblas_sscal(nmols*3, 1.0 / (float)nblocks / (float)nblocksteps, mol_moments_of_inertia, 1);

        // moi summation over all molecules (this sample)
        for (size_t i=0; i<nmols; i++)
        {
            for (size_t abc=0; abc<3; abc++)
            {
                size_t moi_index =
                    + mols_moltypenr[i] * nsamples *   3
                    +                       sample *   3
                    +                                abc;
                moltypes_samples_moments_of_inertia[moi_index] += mol_moments_of_inertia[3*i+abc];
            }
        }
        free(mol_moments_of_inertia);
    }
    verbPrintf(verbosity, "finished all samples\n");
    chfl_trajectory_close(file);

    // divide moi by nmols
    for (size_t h=0; h<nmoltypes; h++)
    {
        cblas_sscal(3*nsamples, 1.0 / (float)moltypes_nmols[h], &moltypes_samples_moments_of_inertia[h*nsamples*3], 1);
    }

    // normalize dos
    for (size_t h=0; h<nmoltypes; h++)
    {
        float norm_factor = 1.0 / (float)nblocks;  // normalize for blocks
        norm_factor *= framelength / (float)nblocksteps;  // normalize DFT
        norm_factor /= (float)moltypes_nmols[h];  // normalize nmols
        DPRINT("norm_factor for moltype %zu: %f\n", h, norm_factor);
        size_t dos_index = h*ndos*nsamples*nfrequencies;
        cblas_sscal(ndos*nsamples*nfrequencies, norm_factor, &moltypes_dos_samples[dos_index], 1);
    }

    // normalize cross spectra
    float norm_factor = 1.0 / (float)nblocks;
    norm_factor *= framelength / (float)nblocksteps;
    cblas_sscal(ncross_spectra*nsamples*nfrequencies, norm_factor, &cross_spectra_samples[0], 1);

    // write dos.json
    result = write_dos(arguments.outfile,
        nsamples,
        nblocksteps,
        nfrequencies,
        framelength,
        ndos,
        ncross_spectra,
        dos_names,
        nmoltypes,
        moltypes_dos_samples,
        cross_spectra_samples,
        moltypes_samples_moments_of_inertia,
        cross_spectra_def);
    if (result != 0) {
        fprintf(stderr, "ERROR: Could not write json to file.\n");
        return 1;
    }

    // free output
    free(moltypes_dos_samples);
    free(cross_spectra_samples);
    free(moltypes_samples_moments_of_inertia);

    // free input arrays
    free_dosparams_arrays(nmoltypes,
            &moltypes_nmols,
            &moltypes_natomspermol,
            &moltypes_atommasses,
            &moltypes_rot_treat,
            &moltypes_abc_indicators,
            &ncross_spectra,
            &cross_spectra_def);

    // free convenience arrays
    free(moltypes_firstmol);
    free(moltypes_firstatom);
    free(mols_moltypenr);
    free(mols_natoms);
    free(mols_mass);
    free(mols_firstatom);

    return 0;
}
