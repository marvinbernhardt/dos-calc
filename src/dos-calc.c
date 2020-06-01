#include "fft.c"
#include "parse-dosparams.c"
#include "structs.h"
#include "trajectory-functions.c"
#include "velocity-decomposition.c"
#include "verbPrintf.c"
#include "write-dos.c"
#include <argp.h>
#include <cblas.h>
#include <chemfiles.h>
#include <math.h>
#include <omp.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

// CLI stuff
const char *argp_program_version = VERSION;
const char *argp_program_bug_address = "<bernhardt@cpc.tu-darmstadt.de>";
static char doc[] =
    "dos-calc -- a programm to calculate densities of states from trajectories";
static char args_doc[] = "dosparams trajectory";

static struct argp_option options[] = {
    {"verbose", 'v', 0, 0, "Produce verbose output", 0},
    {"no-pbc", 'p', 0, 0,
     "Do not recombine molecules seperated by periodic boundary conditions", 0},
    {"outfile", 'o', "FILE", 0, "Output json file. Default: dos.json", 0},
    {"framelength", 'f', "FL", 0,
     "Lenght of each frame in trajectory (in ps). Taken from trajectory if "
     "possible, but can be overwritten with this argument.",
     0},
    {"skip-frames", 's', "SF", 0,
     "Number of frames to be skipped at the beginning of the trajectory.", 0},
    {"refconf", 'e', "FILE", 0,
     "When using Eckart frame, take reference configuration from this file", 0},
    {0}};

struct arguments {
    char *(input_files[2]);
    bool verbosity;
    bool no_pbc;
    char *outfile;
    float framelength;
    unsigned long long skip_frames;
    char *refconf;
};

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
    /* Get the input argument from argp_parse, which we
       know is a pointer to our arguments structure. */
    struct arguments *arguments = state->input;

    switch (key) {
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
    case 's':
        arguments->skip_frames = strtoull(arg, NULL, 10);
        break;
    case 'e':
        arguments->refconf = arg;
        break;

    case ARGP_KEY_ARG:
        /* Too many arguments. */
        if (state->arg_num > 2) {
            argp_usage(state);
        }
        arguments->input_files[state->arg_num] = arg;
        break;

    case ARGP_KEY_END:
        if (state->arg_num < 2)
            /* Not enough arguments. */
            argp_usage(state);
        break;

    default:
        return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

/* Our argp parser. */
static struct argp argp = {options, parse_opt, args_doc, doc, 0, 0, 0};

int main(int argc, char *argv[]) {
    // command line arguments
    struct arguments arguments;

    // TIMING: timings array (parse, read_trj, vel_decomp, fft, write)
    double timings[5] = {0, 0, 0, 0, 0};
    double begin = omp_get_wtime();

    // Default values.
    arguments.verbosity = false;
    arguments.no_pbc = false;
    arguments.outfile = "dos.json";
    arguments.framelength = 0.0;
    arguments.skip_frames = 0;
    arguments.refconf = NULL;

    // parse command line arguments
    argp_parse(&argp, argc, argv, 0, 0, &arguments);

    // define convenience aliases
    bool verbosity = arguments.verbosity;
    char *dosparams_file = arguments.input_files[0];
    char *trajectory_file = arguments.input_files[1];
    char *refconf_file = arguments.refconf;
    verbPrintf(arguments.verbosity, "dosparams file: %s\n", dosparams_file);
    verbPrintf(arguments.verbosity, "trajectory file: %s\n", trajectory_file);
    if (refconf_file) {
        verbPrintf(arguments.verbosity, "refconf file: %s\n", refconf_file);
    }

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
    // NOTE: this calls alloc() functions, use free_dosparams_arrays() in the
    // end
    parse_dosparams(dosparams_file,
                    &nsamples, // output
                    &nblocks, &nblocksteps, &nmoltypes, &moltypes_nmols,
                    &moltypes_natomspermol, &moltypes_atommasses,
                    &moltypes_rot_treat, &moltypes_abc_indicators,
                    &ncross_spectra, &cross_spectra_def);

    if (arguments.verbosity) {
        print_dosparams(nsamples, nblocks, nblocksteps, nmoltypes,
                        moltypes_nmols, moltypes_natomspermol,
                        moltypes_atommasses, moltypes_rot_treat,
                        moltypes_abc_indicators);
    }

    // generate convenience variables and arrays
    size_t natoms = 0;
    size_t nmols = 0;
    unsigned long nfrequencies = nblocksteps / 2 + 1;
    size_t *moltypes_firstmol = calloc(nmoltypes, sizeof(size_t));
    size_t *moltypes_firstatom = calloc(nmoltypes, sizeof(size_t));
    size_t *mols_moltypenr;
    size_t *mols_natoms;
    float *mols_mass;
    size_t *mols_firstatom;
    calc_convenience_variables(
        nmoltypes, moltypes_nmols, moltypes_natomspermol, moltypes_atommasses,
        &natoms, // output
        &nmols, moltypes_firstmol, moltypes_firstatom, &mols_moltypenr,
        &mols_natoms, &mols_mass, &mols_firstatom);

    // open trajectory and one frame for tests
    verbPrintf(verbosity, "testing file %s\n", trajectory_file);
    CHFL_TRAJECTORY *file = chfl_trajectory_open(trajectory_file, 'r');
    CHFL_FRAME *frame = chfl_frame();
    int result = chfl_trajectory_read(file, frame);
    if ((file == NULL) || (result != CHFL_SUCCESS)) {
        fprintf(stderr, "ERROR: Reading trajectory failed.\n");
        return 1;
    }
    // check for number of atoms
    check_frame_natoms(frame, natoms);
    // check for velocities
    check_frame_velocities(frame);
    // check for orthorombic box
    check_frame_orthorombic_box(frame, arguments.no_pbc);
    // get framelength
    float framelength;
    if (arguments.framelength == 0.0) {
        framelength = get_traj_framelength(file, frame);
    } else {
        framelength = arguments.framelength;
    }
    chfl_trajectory_close(file);
    verbPrintf(verbosity, "framelength is %f ps\n", framelength);

    // test if refconf file is needed but not given and vice versa
    bool refconf_needed = false;
    for (size_t h = 0; h < nmoltypes; h++) {
        if (moltypes_rot_treat[h] == 'e' || moltypes_rot_treat[h] == 'p' ||
            moltypes_rot_treat[h] == 'E' || moltypes_rot_treat[h] == 'P') {
            refconf_needed = true;
        }
    }
    if (!refconf_file && refconf_needed) {
        fprintf(stderr,
                "ERROR: No refconf file given but needed for some molecules\n");
        return 1;
    }
    if (refconf_file && !refconf_needed) {
        fprintf(
            stderr,
            "WARNING: Refconf file provided but not needed for any moltype\n");
    }

    // open, test, and evaluate refconf file
    float *refconf_pos = calloc(natoms * 3, sizeof(float));
    float *refconf_box = calloc(3, sizeof(float));
    float *atom_refpos_principal_components = calloc(natoms * 3, sizeof(float));
    if (refconf_file) {
        verbPrintf(verbosity, "start reading refconf\n");
        file = chfl_trajectory_open(refconf_file, 'r');
        int result = chfl_trajectory_read(file, frame);
        if ((file == NULL) || (result != CHFL_SUCCESS)) {
            fprintf(stderr, "ERROR: Reading refconf file failed.\n");
            return 1;
        }
        // check for number of atoms
        check_frame_natoms(frame, natoms);
        // check for orthorombic box
        check_frame_orthorombic_box(frame, arguments.no_pbc);
        get_frame_pos_box(frame, natoms, refconf_pos, refconf_box);
        chfl_trajectory_close(file);
        // get principal axis
        calculate_refpos_principal_components(
            refconf_pos, refconf_box, nmols, mols_firstatom, mols_natoms,
            mols_moltypenr, moltypes_atommasses, mols_mass, moltypes_rot_treat,
            arguments.no_pbc,
            atom_refpos_principal_components); // output
    }

    // open trajectory for calculations
    file = chfl_trajectory_open(trajectory_file, 'r');
    // skip frames
    verbPrintf(verbosity, "skipping %llu frames\n", arguments.skip_frames);
    for (size_t t = 0; t < arguments.skip_frames; t++) {
        chfl_trajectory_read(file, frame);
    }

    // output arrays
    const size_t ndos = 15;
    const char *dos_names[15] = {"trn_x",  "trn_y",  "trn_z",  "rot_x",
                                 "rot_y",  "rot_z",  "vib_x",  "vib_y",
                                 "vib_z",  "roto_a", "roto_b", "roto_c",
                                 "vibc_x", "vibc_y", "vibc_z"};
    // order is: trn_xyz, rot_xyz, vib_xyz, rot_omega_abc
    // this one has all the dos in it
    float *moltypes_dos_samples =
        calloc(nmoltypes * ndos * nsamples * nfrequencies, sizeof(float));
    // this one all cross spectra
    float *cross_spectra_samples =
        calloc(ncross_spectra * nsamples * nfrequencies, sizeof(float));
    // moments of inertia and its std. deviation
    float *moltypes_samples_moments_of_inertia =
        calloc(nmoltypes * nsamples * 3, sizeof(float));
    float *moltypes_samples_moments_of_inertia_squared =
        calloc(nmoltypes * nsamples * 3, sizeof(float));
    float *moltypes_samples_moments_of_inertia_std =
        calloc(nmoltypes * nsamples * 3, sizeof(float));
    // coriolis energy term
    float *moltypes_samples_coriolis =
        calloc(nmoltypes * nsamples, sizeof(float));

    // TIMING: parse end
    timings[0] += omp_get_wtime() - begin;
    begin = omp_get_wtime();

    // start samples loop
    verbPrintf(verbosity, "going through %zu samples\n", nsamples);
    for (size_t sample = 0; sample < nsamples; sample++) {
        verbPrintf(verbosity, "now doing sample %zu\n", sample);

        // moi/coiolis of mols (this sample)
        float *mol_moments_of_inertia = calloc(nmols * 3, sizeof(float));
        float *mol_moments_of_inertia_squared =
            calloc(nmols * 3, sizeof(float));
        float *mol_coriolis = calloc(nmols, sizeof(float));

        // start block loop
        verbPrintf(verbosity, "going through %zu blocks\n", nblocks);
        for (size_t block = 0; block < nblocks; block++) {
            verbPrintf(verbosity, "now doing block %zu\n", block);

            verbPrintf(verbosity, "start reading trajectory block\n");
            float *block_pos = calloc(natoms * 3 * nblocksteps, sizeof(float));
            float *block_vel = calloc(natoms * 3 * nblocksteps, sizeof(float));
            float *block_box = calloc(3 * nblocksteps, sizeof(float));
            get_traj_pos_vel_box(file, nblocksteps, natoms, block_pos,
                                 block_vel, block_box);

            // TIMING: read_trj end
            timings[1] += omp_get_wtime() - begin;
            begin = omp_get_wtime();

            verbPrintf(verbosity, "start decomposition\n");
            // series to be fourier transformed per block
            float *mol_velocities_sqrt_m_trn =
                calloc(nmols * 3 * nblocksteps, sizeof(float));
            float *mol_omegas_sqrt_i_rot =
                calloc(nmols * 3 * nblocksteps, sizeof(float));
            float *atom_velocities_sqrt_m_vib =
                calloc(natoms * 3 * nblocksteps, sizeof(float));
            float *atom_velocities_sqrt_m_rot =
                calloc(natoms * 3 * nblocksteps, sizeof(float));
            float *atom_velocities_sqrt_m_vibc =
                calloc(natoms * 3 * nblocksteps, sizeof(float));
            // per block vectors
            float *mol_block_moments_of_inertia =
                calloc(nmols * 3 * nblocksteps, sizeof(float));
            float *mol_block_moments_of_inertia_squared =
                calloc(nmols * 3 * nblocksteps, sizeof(float));
            // per block numbers
            float *mol_block_coriolis =
                calloc(nmols * nblocksteps, sizeof(float));
            decompose_velocities(
                block_pos, block_vel, block_box, nblocksteps, natoms, nmols,
                mols_firstatom, mols_natoms, mols_moltypenr,
                moltypes_atommasses, mols_mass, moltypes_rot_treat,
                moltypes_abc_indicators, arguments.no_pbc,
                atom_refpos_principal_components,
                mol_velocities_sqrt_m_trn, // output
                mol_omegas_sqrt_i_rot, atom_velocities_sqrt_m_vib,
                atom_velocities_sqrt_m_rot, atom_velocities_sqrt_m_vibc,
                mol_block_moments_of_inertia,
                mol_block_moments_of_inertia_squared, mol_block_coriolis);

            // TIMING: vel_decomp end
            timings[2] += omp_get_wtime() - begin;
            begin = omp_get_wtime();

            verbPrintf(verbosity, "start DoS calculation (FFT)\n");
            dos_calculation(
                nmoltypes, nblocksteps, nfrequencies, moltypes_firstmol,
                moltypes_firstatom, moltypes_nmols, moltypes_natomspermol,
                mol_velocities_sqrt_m_trn, mol_omegas_sqrt_i_rot,
                atom_velocities_sqrt_m_vib, atom_velocities_sqrt_m_rot,
                atom_velocities_sqrt_m_vibc, ndos, nsamples, sample,
                ncross_spectra, cross_spectra_def,
                moltypes_dos_samples, // output
                cross_spectra_samples);

            // moi summation over all nblocksteps (this block)
            for (size_t i = 0; i < nmols; i++) {
                for (size_t t = 0; t < nblocksteps; t++) {
                    for (size_t abc = 0; abc < 3; abc++) {
                        mol_moments_of_inertia[3 * i + abc] +=
                            mol_block_moments_of_inertia[3 * nblocksteps * i +
                                                         nblocksteps * abc + t];
                        mol_moments_of_inertia_squared[3 * i + abc] +=
                            mol_block_moments_of_inertia_squared
                                [3 * nblocksteps * i + nblocksteps * abc + t];
                    }
                }
            }

            // coriolis summation over all nblocksteps (this block)
            for (size_t i = 0; i < nmols; i++) {
                for (size_t t = 0; t < nblocksteps; t++) {
                    mol_coriolis[i] += mol_block_coriolis[nblocksteps * i + t];
                }
            }

            // TIMING: fft end
            timings[3] += omp_get_wtime() - begin;
            begin = omp_get_wtime();

            // free block stuff
            free(block_pos);
            free(block_vel);
            free(block_box);
            free(mol_velocities_sqrt_m_trn);
            free(mol_omegas_sqrt_i_rot);
            free(atom_velocities_sqrt_m_vib);
            free(atom_velocities_sqrt_m_rot);
            free(atom_velocities_sqrt_m_vibc);
            free(mol_block_moments_of_inertia);
            free(mol_block_moments_of_inertia_squared);
            free(mol_block_coriolis);
        }
        verbPrintf(verbosity, "finished all blocks\n");

        // divide moi by number of blocks and number of blocksteps
        cblas_sscal(nmols * 3, 1.0 / (float)nblocks / (float)nblocksteps,
                    mol_moments_of_inertia, 1);
        cblas_sscal(nmols * 3, 1.0 / (float)nblocks / (float)nblocksteps,
                    mol_moments_of_inertia_squared, 1);
        // divide coriolis by number of blocks and number of blocksteps
        cblas_sscal(nmols, 1.0 / (float)nblocks / (float)nblocksteps,
                    mol_coriolis, 1);

        // moi summation over all molecules (this sample)
        for (size_t i = 0; i < nmols; i++) {
            for (size_t abc = 0; abc < 3; abc++) {
                size_t moi_index =
                    mols_moltypenr[i] * nsamples * 3 + sample * 3 + abc;
                moltypes_samples_moments_of_inertia[moi_index] +=
                    mol_moments_of_inertia[3 * i + abc];
                moltypes_samples_moments_of_inertia_squared[moi_index] +=
                    mol_moments_of_inertia_squared[3 * i + abc];
            }
        }
        // coriolis summation over all molecules (this sample)
        for (size_t i = 0; i < nmols; i++) {
            size_t coriolis_index = mols_moltypenr[i] * nsamples + sample;
            moltypes_samples_coriolis[coriolis_index] += mol_coriolis[i];
        }
        free(mol_moments_of_inertia);
        free(mol_moments_of_inertia_squared);
        free(mol_coriolis);
    }
    verbPrintf(verbosity, "finished all samples\n");
    chfl_trajectory_close(file);

    // divide moi by nmols
    for (size_t h = 0; h < nmoltypes; h++) {
        cblas_sscal(3 * nsamples, 1.0 / (float)moltypes_nmols[h],
                    &moltypes_samples_moments_of_inertia[h * nsamples * 3], 1);
    }
    // divide moi squares by nmols
    for (size_t h = 0; h < nmoltypes; h++) {
        cblas_sscal(
            3 * nsamples, 1.0 / (float)moltypes_nmols[h],
            &moltypes_samples_moments_of_inertia_squared[h * nsamples * 3], 1);
    }
    // calculate std. deviation of moi
    for (size_t q = 0; q < nmoltypes * nsamples * 3; q++) {
        moltypes_samples_moments_of_inertia_std[q] =
            sqrtf(moltypes_samples_moments_of_inertia_squared[q] -
                  powf(moltypes_samples_moments_of_inertia[q], 2.0));
    }

    // average Coriolis energy term
    for (size_t h = 0; h < nmoltypes; h++) {
        cblas_sscal(nsamples, 1.0 / (float)moltypes_nmols[h],
                    &moltypes_samples_coriolis[h * nsamples], 1);
    }

    // normalize dos
    for (size_t h = 0; h < nmoltypes; h++) {
        float norm_factor = 1.0 / (float)nblocks;        // normalize for blocks
        norm_factor *= framelength / (float)nblocksteps; // normalize DFT
        norm_factor /= (float)moltypes_nmols[h];         // normalize nmols
        size_t dos_index = h * ndos * nsamples * nfrequencies;
        cblas_sscal(ndos * nsamples * nfrequencies, norm_factor,
                    &moltypes_dos_samples[dos_index], 1);
    }

    // normalize cross spectra
    float norm_factor = 1.0 / (float)nblocks;
    norm_factor *= framelength / (float)nblocksteps;
    cblas_sscal(ncross_spectra * nsamples * nfrequencies, norm_factor,
                &cross_spectra_samples[0], 1);

    // write dos.json
    result = write_dos(arguments.outfile, nsamples, nblocksteps, nfrequencies,
                       framelength, ndos, ncross_spectra, dos_names, nmoltypes,
                       moltypes_dos_samples, cross_spectra_samples,
                       moltypes_samples_moments_of_inertia,
                       moltypes_samples_moments_of_inertia_std,
                       cross_spectra_def, moltypes_samples_coriolis);
    if (result != 0) {
        fprintf(stderr, "ERROR: Could not write json to file.\n");
        exit(1);
    }

    // free output
    free(moltypes_dos_samples);
    free(cross_spectra_samples);
    free(moltypes_samples_moments_of_inertia);
    free(moltypes_samples_moments_of_inertia_squared);
    free(moltypes_samples_moments_of_inertia_std);
    free(moltypes_samples_coriolis);

    // free input arrays
    free_dosparams_arrays(nmoltypes, &moltypes_nmols, &moltypes_natomspermol,
                          &moltypes_atommasses, &moltypes_rot_treat,
                          &moltypes_abc_indicators, &ncross_spectra,
                          &cross_spectra_def);

    // free convenience arrays
    free(moltypes_firstmol);
    free(moltypes_firstatom);
    free(mols_moltypenr);
    free(mols_natoms);
    free(mols_mass);
    free(mols_firstatom);

    // free other
    chfl_free(frame);
    free(refconf_pos);
    free(refconf_box);
    free(atom_refpos_principal_components);

    // TIMING: write end
    timings[4] += omp_get_wtime() - begin;

    verbPrintf(arguments.verbosity, "timings in seconds:\n", timings[0]);
    verbPrintf(arguments.verbosity, "parsing input: %g\n", timings[0]);
    verbPrintf(arguments.verbosity, "reading trajectory: %g\n", timings[1]);
    verbPrintf(arguments.verbosity, "velocity decomposition: %g\n", timings[2]);
    verbPrintf(arguments.verbosity, "fast Fourier transform: %g\n", timings[3]);
    verbPrintf(arguments.verbosity, "writing output: %g\n", timings[4]);

    return 0;
}
