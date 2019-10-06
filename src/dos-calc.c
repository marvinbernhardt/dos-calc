#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdbool.h>
#include <argp.h>
#include <cblas.h>
#include <math.h>
#include "xdrfile.h"
#include "xdrfile_trr.h"
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
    { 0 }
};

struct arguments
{
    char *(input_files[2]);
    bool verbosity;
    bool no_pbc;
    char *outfile;
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
    unsigned long nfftsteps = nblocksteps / 2 + 1;
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

    // check atom number
    // warning/error when more or less natoms in traj
    verbPrintf(verbosity, "starting with file %s\n", trajectory_file);
    int natoms_traj_int;
    if (read_trr_natoms(trajectory_file , &natoms_traj_int) != 0) {
        fprintf(stderr, "ERROR: could not read number of atoms from %s\n", trajectory_file);
        return 1;
    }
    size_t natoms_traj = (size_t) natoms_traj_int;
    if (natoms > natoms_traj) {
        fprintf(stderr, "ERROR: The topology you give has more atoms than first frame of the trajectory\n");
        return 1;
    }
    else if (natoms < natoms_traj) {
        fprintf(stderr, "WARNING: The topology you give has less atoms than first frame of the trajectory\n");
        fprintf(stderr, "         Some atoms are ignored in every frame\n");
    }

    // open file and check first two frames
    XDRFILE* traj = xdrfile_open(trajectory_file, "r");
    if (traj == NULL) {
        fprintf(stderr, "ERROR: could not open file %s\n", trajectory_file);
        return 1;
    }
    verbPrintf(verbosity, "reading framelength from first two frame's timestamps: ");
    float framelength;
    int step;
    float time0, time1;
    float lambda;
    matrix box;
    rvec* r;
    rvec* v;
    int has_prop;
    r = calloc(natoms_traj, sizeof(*r));
    v = calloc(natoms_traj, sizeof(*v));
    int result = 0;
    result = read_trr(traj, natoms_traj, &step, &time0, &lambda, box, r, v, NULL, &has_prop);
    if (result != 0) {
        fprintf(stderr, "ERROR: First frame of trajectory broken\n");
        return 1;
    }
    result = read_trr(traj, natoms_traj, &step, &time1, &lambda, box, r, v, NULL, &has_prop);
    if (result != 0) {
        fprintf(stderr, "ERROR: Second frame of trajectory broken\n");
        return 1;
    }
    framelength = time1 - time0;
    verbPrintf(verbosity, "%f\n", framelength);
    free(r);
    free(v);
    // go back to start of file
    xdrfile_close(traj);
    traj = xdrfile_open(trajectory_file , "r");

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
      "rot_a",
      "rot_b",
      "rot_c"
    };
    // order is: trn_xyz, rot_xyz, vib_xyz, rot_omega_abc
    float* moltypes_dos_samples = calloc(nmoltypes*ndos*nsamples*nfftsteps, sizeof(float));
    float* cross_spectra_samples = calloc(ncross_spectra*nsamples*nfftsteps, sizeof(float));
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

            decomposeVelocities (traj,
                    nblocksteps,
                    natoms_traj,
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
                    nfftsteps,
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
    xdrfile_close(traj);

    // divide moi by nmols
    for (size_t h=0; h<nmoltypes; h++)
    {
        cblas_sscal(3*nsamples, 1.0 / (float)moltypes_nmols[h], &moltypes_samples_moments_of_inertia[h*nsamples*3], 1);
    }

    // normalize dos
    for (size_t h=0; h<nmoltypes; h++)
    {
        float norm_factor = 1.0 / (float)nblocks;
        norm_factor *= 2.0 * framelength / (float)nblocksteps / (float)moltypes_nmols[h];
        DPRINT("norm_factor for moltype %zu: %f\n", h, norm_factor);
        size_t dos_index = h*ndos*nsamples*nfftsteps;
        cblas_sscal(ndos*nsamples*nfftsteps, norm_factor, &moltypes_dos_samples[dos_index], 1);
    }

    // normalize cross spectra
    float norm_factor = 1.0 / (float)nblocks;
    norm_factor *= 2.0 * framelength / (float)nblocksteps;
    cblas_sscal(ncross_spectra*nsamples*nfftsteps, norm_factor, &cross_spectra_samples[0], 1);

    // write dos.json
    result = write_dos(arguments.outfile,
        nsamples,
        nblocksteps,
        nfftsteps,
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
