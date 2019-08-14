#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdbool.h>
#include <argp.h>
#include <cblas.h>
#include "xdrfile.h"
#include "xdrfile_trr.h"
#include "velocity-decomposition.c"
#include "fft.c"
#include "verbPrintf.c"

// constants
#define K_GRO 0.00831445986144858  // k_B in gromacs units

// CLI stuff
const char* argp_program_version = VERSION;
const char* argp_program_bug_address = "<bernhardt@cpc.tu-darmstadt.de>";
static char doc[] = "dos-calc -- a programm to calculate densities of states from trajectories";
static char args_doc[] = "";

static struct argp_option options[] = {
    {"dump",        'd', 0,      0, "Dump velocities of first block"},
    {"verbose",     'v', 0,      0, "Produce verbose output"},
    {"components",  'c', 0,      0, "Compute xyz components of trn and vib and abc of rot"},
    {"cross",       'x', 0,      0, "Compute cross-spectrum of translational with rotationalas DoF"},
    {"rot-alt",     'a', 0,      0, "Additionally output alternative DoS_rot, calculated from v_rot"},
    {"file",        'f', "FILE", 0, "Input .trr trajectory file (default: traj.trr)"},
    {"no-pbc",      'p', 0,      0, "Do not recombine molecules seperated by periodic boundary conditions"},
    {"temperature", 't', "TEMP", 0, "Temperature in trajectory (in K); if given, the DoS normalized"},
    { 0 }
};

struct arguments
{
    bool dump_vel;
    bool verbosity;
    bool calc_components;
    bool calc_cross;
    bool calc_rot_alt;
    char *file;
    bool no_pbc;
    float temperature;
};

static error_t parse_opt (int key, char *arg, struct argp_state *state)
{
    /* Get the input argument from argp_parse, which we
       know is a pointer to our arguments structure. */
    struct arguments *arguments = state->input;

    switch (key)
    {
        case 'd':
            arguments->dump_vel = true;
            break;
        case 'v':
            arguments->verbosity = true;
            break;
        case 'c':
            arguments->calc_components = true;
            break;
        case 'x':
            arguments->calc_cross = true;
            break;
        case 'a':
            arguments->calc_rot_alt = true;
            break;
        case 'f':
            arguments->file = arg;
            break;
        case 'p':
            arguments->no_pbc = true;
            break;
        case 't':
            arguments->temperature = strtof(arg, NULL);
            break;

        case ARGP_KEY_ARG:
            if (state->arg_num >= 0)
                /* Too many arguments. */
                argp_usage (state);

            break;

        case ARGP_KEY_END:
            if (state->arg_num < 0)
                /* Not enough arguments. */
                argp_usage (state);
            break;

        default:
            return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

/* Our argp parser. */
static struct argp argp = {options, parse_opt, args_doc, doc};

int main( int argc, char *argv[] )
{
    // command line arguments
    struct arguments arguments;

    // Default values.
    arguments.dump_vel = false;
    arguments.verbosity = false;
    arguments.calc_components = false;
    arguments.calc_cross = false;
    arguments.calc_rot_alt = false;
    arguments.file = "traj.trr";
    arguments.no_pbc = false;
    arguments.temperature = 0.0;

    // parse command line arguments
    argp_parse(&argp, argc, argv, 0, 0, &arguments);

    bool verbosity = arguments.verbosity;
    bool normalize;
    if (arguments.temperature != 0.0)
    {
        normalize = true;
        verbPrintf(verbosity, "Temperature given. DoS will be normalized\n");
    }
    else
    {
        normalize = false;
        verbPrintf(verbosity, "Temperature missing. DoS will not be normalized\n");
    }

    // input that will be scanned
    int nsamples;
    int nblocks;
    long nblocksteps;
    int nmoltypes;

    // convenience variables
    int natoms;
    int nmols;

    // file pointers for output
    FILE* f;
    FILE* fa;
    FILE* fb;
    FILE* fc;
    FILE* ftr;
    FILE* ftv;
    FILE* frv;

    // start scanning
    verbPrintf(verbosity, "Start reading user input\n");
    verbPrintf(verbosity, "Reading integers\n");

    verbPrintf(verbosity, "nsamples: ");
    if (scanf("%d", &nsamples) != 1)
    {
        fprintf(stderr, "ERROR: Failed to read nsamples.\n");
        return 1;
    }
    verbPrintf(verbosity, "%d\n", nsamples);

    verbPrintf(verbosity, "nblocks: ");
    if (scanf("%d", &nblocks) != 1)
    {
        fprintf(stderr, "ERROR: Failed to read nblocks.\n");
        return 1;
    }
    verbPrintf(verbosity, "%d\n", nblocks);

    verbPrintf(verbosity, "nblocksteps: ");
    if (scanf("%ld", &nblocksteps) != 1)
    {
        fprintf(stderr, "ERROR: Failed to read nblocksteps.\n");
        return 1;
    }
    verbPrintf(verbosity, "%ld\n", nblocksteps);
    int nfftsteps = nblocksteps / 2 + 1;

    verbPrintf(verbosity, "nmoltypes: ");
    if (scanf("%d", &nmoltypes) != 1)
    {
        fprintf(stderr, "ERROR: Failed to read nmoltypes.\n");
        return 1;
    }
    verbPrintf(verbosity, "%d\n", nmoltypes);

    // scanning moltypes
    int* moltypes_nmols = calloc(nmoltypes, sizeof(int));
    int* moltypes_natomspermol = calloc(nmoltypes, sizeof(int));
    float** moltypes_atommasses = (float**) malloc(nmoltypes * sizeof(float *));
    char* moltypes_rot_treat = calloc(nmoltypes, sizeof(char));
    int** moltypes_abc_indicators = (int**) malloc(nmoltypes * sizeof(int *));

    for (int h=0; h<nmoltypes; h++)
    {
        verbPrintf(verbosity, "now reading moltype %d of %d\n", h+1, nmoltypes);

        verbPrintf(verbosity, "moltype %d nmols: ", h+1);
        if (scanf("%d", &moltypes_nmols[h]) != 1)
        {
            fprintf(stderr, "ERROR: Failed to read nmols.\n");
            return 1;
        }
        verbPrintf(verbosity, "%d\n", moltypes_nmols[h]);

        verbPrintf(verbosity, "moltype %d natomspermol: ", h+1);
        if (scanf("%d", &moltypes_natomspermol[h]) != 1)
        {
            fprintf(stderr, "ERROR: Failed to read natomspermol.\n");
            return 1;
        }
        verbPrintf(verbosity, "%d\n", moltypes_natomspermol[h]);

        verbPrintf(verbosity, "moltype %d atommasses: ", h+1);
        moltypes_atommasses[h] = (float*) malloc(moltypes_natomspermol[h] * sizeof(float));
        for (int j=0; j<moltypes_natomspermol[h]; j++)
        {
            if (scanf("%f", &moltypes_atommasses[h][j]) != 1)
            {
                fprintf(stderr, "ERROR: Failed to read atommass.\n");
                return 1;
            }
        }
        for (int j=0; j<moltypes_natomspermol[h]; j++)
        {
            verbPrintf(verbosity, "%f ", moltypes_atommasses[h][j]);
        }
        verbPrintf(verbosity, "\n");

        verbPrintf(verbosity, "moltype %d rot_treat: ", h+1);
        if (scanf(" %c", &moltypes_rot_treat[h]) != 1)
        {
            fprintf(stderr, "ERROR: Failed to read rot_treat.\n");
            return 1;
        }
        verbPrintf(verbosity, "%c\n", moltypes_rot_treat[h]);

        verbPrintf(verbosity, "moltype %d abc_indicators: ", h+1);
        moltypes_abc_indicators[h] = (int*) malloc(4 * sizeof(int));
        for (int j=0; j<4; j++)
        {
            if (scanf("%d", &moltypes_abc_indicators[h][j]) != 1)
            {
                fprintf(stderr, "ERROR: Failed to read abc_indicators.\n");
                return 1;
            }
        }
        verbPrintf(verbosity, "%d %d %d %d\n", moltypes_abc_indicators[h][0], moltypes_abc_indicators[h][1], moltypes_abc_indicators[h][2], moltypes_abc_indicators[h][3]);
    }


    // generate convenience variables and arrays
    nmols = 0;
    natoms = 0;
    int* moltypes_firstmol = calloc(nmoltypes, sizeof(int));
    int* moltypes_firstatom = calloc(nmoltypes, sizeof(int));

    for (int h=0; h<nmoltypes; h++)
    {
        moltypes_firstmol[h] = nmols;
        moltypes_firstatom[h] = natoms;
        nmols += moltypes_nmols[h];
        natoms += moltypes_nmols[h] * moltypes_natomspermol[h];
    }

    int* mols_moltypenr = calloc(nmols, sizeof(int));
    int* mols_natoms = calloc(nmols, sizeof(int));
    float* mols_mass = calloc(nmols, sizeof(float));

    for (int i=0; i<nmols; i++)
    {
        for (int h=0; h<nmoltypes; h++)
        {
            if (i >= moltypes_firstmol[h] && i < moltypes_firstmol[h] + moltypes_nmols[h])
                mols_moltypenr[i] = h;
        }
        mols_natoms[i] = moltypes_natomspermol[mols_moltypenr[i]];
        mols_mass[i] = 0;
        for (int j=0; j<mols_natoms[i]; j++)
            mols_mass[i] += moltypes_atommasses[mols_moltypenr[i]][j];
    }

    int* mols_firstatom = calloc(nmols, sizeof(int));
    for (int h=0; h<nmoltypes; h++)
    {
        for (int ii=0; ii<moltypes_nmols[h]; ii++)
        {
            mols_firstatom[moltypes_firstmol[h]+ii] = moltypes_firstatom[h] + moltypes_natomspermol[h] * ii;
        }
    }

    // read natoms_traj
    verbPrintf(verbosity, "starting with file %s\n", arguments.file);
    int natoms_traj;
    if (read_trr_natoms(arguments.file , &natoms_traj) != 0) {
        fprintf(stderr, "ERROR: could not read number of atoms from %s\n", arguments.file);
        return 1;
    }

    // warning when more or less natoms in traj
    if (natoms > natoms_traj) {
        fprintf(stderr, "ERROR: The topology you give has more atoms than first frame of the trajectory\n");
        return 1;
    }
    else if (natoms < natoms_traj) {
        fprintf(stderr, "WARNING: The topology you give has less atoms than first frame of the trajectory\n");
        fprintf(stderr, "         Some atoms are ignored in every frame\n");
    }

    // open file
    XDRFILE* traj = xdrfile_open(arguments.file , "r");
    if (traj == NULL) {
        fprintf(stderr, "ERROR: could not open file %s\n", arguments.file);
        return 1;
    }

    // read first two frames
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
    traj = xdrfile_open(arguments.file , "r");

    verbPrintf(verbosity, "going through %d samples\n", nsamples);
    for (int sample=0; sample<nsamples; sample++)
    {
        verbPrintf(verbosity, "now doing sample %d\n", sample);

        // output arrays
        float* mol_moments_of_inertia = calloc(nmols*3, sizeof(float));
        float* moltypes_moments_of_inertia = calloc(nmoltypes*3, sizeof(float));
        float* moltypes_dos_raw_trn = calloc(nmoltypes*nfftsteps, sizeof(float));
        float* moltypes_dos_raw_trn_x = calloc(nmoltypes*nfftsteps, sizeof(float));
        float* moltypes_dos_raw_trn_y = calloc(nmoltypes*nfftsteps, sizeof(float));
        float* moltypes_dos_raw_trn_z = calloc(nmoltypes*nfftsteps, sizeof(float));
        float* moltypes_dos_raw_rot = calloc(nmoltypes*nfftsteps, sizeof(float));
        float* moltypes_dos_raw_rot_a = calloc(nmoltypes*nfftsteps, sizeof(float));
        float* moltypes_dos_raw_rot_b = calloc(nmoltypes*nfftsteps, sizeof(float));
        float* moltypes_dos_raw_rot_c = calloc(nmoltypes*nfftsteps, sizeof(float));
        float* moltypes_dos_raw_rot_alt = calloc(nmoltypes*nfftsteps, sizeof(float)); // alternative dos_rot
        float* moltypes_dos_raw_rot_alt_x = calloc(nmoltypes*nfftsteps, sizeof(float));
        float* moltypes_dos_raw_rot_alt_y = calloc(nmoltypes*nfftsteps, sizeof(float));
        float* moltypes_dos_raw_rot_alt_z = calloc(nmoltypes*nfftsteps, sizeof(float));
        float* moltypes_dos_raw_vib = calloc(nmoltypes*nfftsteps, sizeof(float));
        float* moltypes_dos_raw_vib_x = calloc(nmoltypes*nfftsteps, sizeof(float));
        float* moltypes_dos_raw_vib_y = calloc(nmoltypes*nfftsteps, sizeof(float));
        float* moltypes_dos_raw_vib_z = calloc(nmoltypes*nfftsteps, sizeof(float));
        // cross spectra
        float* moltypes_dos_raw_x_trn_rot = calloc(nmoltypes*nfftsteps, sizeof(float));
        float* moltypes_dos_raw_x_trn_vib = calloc(nmoltypes*nfftsteps, sizeof(float));
        float* moltypes_dos_raw_x_rot_vib = calloc(nmoltypes*nfftsteps, sizeof(float));

        verbPrintf(verbosity, "going through %d blocks\n", nblocks);
        for (int block=0; block<nblocks; block++)
        {
            verbPrintf(verbosity, "now doing block %d\n", block);
            verbPrintf(verbosity, "start decomposition\n");

            float* mol_velocities_sqrt_m_trn = calloc(nmols*3*nblocksteps, sizeof(float));
            float* mol_omegas_sqrt_i_rot = calloc(nmols*3*nblocksteps, sizeof(float));
            float* atom_velocities_sqrt_m_vib = calloc(natoms*3*nblocksteps, sizeof(float));
            float* atom_velocities_sqrt_m_rot = calloc(natoms*3*nblocksteps, sizeof(float));

            decomposeVelocities (traj,
                    nblocksteps,
                    natoms,
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

            // dump first block only!
            if (arguments.dump_vel == 1 && block == 0)
            {
                verbPrintf(verbosity, "start dumping (first block only)\n");
                f = fopen("mol_omegas_sqrt_i_rot.txt", "w");
                for (int h=0; h<nmoltypes; h++)
                {
                    int first_dof = moltypes_firstmol[h]*3;
                    int last_dof = moltypes_firstmol[h]*3 + moltypes_nmols[h]*3;
                    for (int i=first_dof; i<last_dof; i++)
                    {
                        for (int t=0; t<nblocksteps; t++)
                        {
                            if (t!=0) fprintf(f, " ");
                            fprintf(f, "%f", mol_omegas_sqrt_i_rot[i*nblocksteps + t]);
                        }
                        fprintf(f, "\n");
                    }
                }
                fclose(f);

                f = fopen("mol_velocities_sqrt_m_trn.txt", "w");
                for (int h=0; h<nmoltypes; h++)
                {
                    int first_dof = moltypes_firstmol[h]*3;
                    int last_dof = moltypes_firstmol[h]*3 + moltypes_nmols[h]*3;
                    for (int i=first_dof; i<last_dof; i++)
                    {
                        for (int t=0; t<nblocksteps; t++)
                        {
                            if (t!=0) fprintf(f, " ");
                            fprintf(f, "%f", mol_velocities_sqrt_m_trn[i*nblocksteps + t]);
                        }
                        fprintf(f, "\n");
                    }
                }
                fclose(f);

                f = fopen("atom_velocities_sqrt_m_vib.txt", "w");
                for (int h=0; h<nmoltypes; h++)
                {
                    int first_dof_vib = moltypes_firstatom[h]*3;
                    int last_dof_vib = moltypes_firstatom[h]*3 + moltypes_nmols[h]*moltypes_natomspermol[h]*3;
                    for (int i=first_dof_vib; i<last_dof_vib; i++)
                    {
                        for (int t=0; t<nblocksteps; t++)
                        {
                            if (t!=0) fprintf(f, " ");
                            fprintf(f, "%f", atom_velocities_sqrt_m_vib[i*nblocksteps + t]);
                        }
                        fprintf(f, "\n");
                    }
                }
                fclose(f);

                f = fopen("atom_velocities_sqrt_m_rot.txt", "w");
                for (int h=0; h<nmoltypes; h++)
                {
                    int first_dof_vib = moltypes_firstatom[h]*3;
                    int last_dof_vib = moltypes_firstatom[h]*3 + moltypes_nmols[h]*moltypes_natomspermol[h]*3;
                    for (int i=first_dof_vib; i<last_dof_vib; i++)
                    {
                        for (int t=0; t<nblocksteps; t++)
                        {
                            if (t!=0) fprintf(f, " ");
                            fprintf(f, "%f", atom_velocities_sqrt_m_rot[i*nblocksteps + t]);
                        }
                        fprintf(f, "\n");
                    }
                }
                fclose(f);
            }

            verbPrintf(verbosity, "start DoS calculation (FFT)\n");
            DOSCalculation (nmoltypes,
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
                    arguments.calc_components,
                    arguments.calc_cross,
                    arguments.calc_rot_alt,
                    moltypes_dos_raw_trn, //output
                    moltypes_dos_raw_trn_x,
                    moltypes_dos_raw_trn_y,
                    moltypes_dos_raw_trn_z,
                    moltypes_dos_raw_rot,
                    moltypes_dos_raw_rot_a,
                    moltypes_dos_raw_rot_b,
                    moltypes_dos_raw_rot_c,
                    moltypes_dos_raw_rot_alt,
                    moltypes_dos_raw_rot_alt_x,
                    moltypes_dos_raw_rot_alt_y,
                    moltypes_dos_raw_rot_alt_z,
                    moltypes_dos_raw_vib,
                    moltypes_dos_raw_vib_x,
                    moltypes_dos_raw_vib_y,
                    moltypes_dos_raw_vib_z,
                    moltypes_dos_raw_x_trn_rot,
                    moltypes_dos_raw_x_trn_vib,
                    moltypes_dos_raw_x_rot_vib);

            //free velocities
            free(mol_velocities_sqrt_m_trn);
            free(mol_omegas_sqrt_i_rot);
            free(atom_velocities_sqrt_m_vib);
            free(atom_velocities_sqrt_m_rot);
        }
        verbPrintf(verbosity, "finished all blocks\n");

        // divide moi by number of blocks and number of blocksteps
        cblas_sscal(nmols*3, 1.0 / (float) nblocks / (float) nblocksteps, mol_moments_of_inertia, 1);

        // normalize
        float norm_factor;
        for (int h=0; h<nmoltypes; h++)
        {
            if (normalize) {
               norm_factor = 1.0 / (float)nblocks;
               norm_factor *= 2.0 * framelength / nblocksteps / moltypes_nmols[h] / K_GRO / arguments.temperature;
            }
            else
            {
               norm_factor = 1.0 / (float)nblocks;
            }
            cblas_sscal(nfftsteps, norm_factor, &moltypes_dos_raw_trn[h*nfftsteps], 1);
            cblas_sscal(nfftsteps, norm_factor, &moltypes_dos_raw_trn_x[h*nfftsteps], 1);
            cblas_sscal(nfftsteps, norm_factor, &moltypes_dos_raw_trn_y[h*nfftsteps], 1);
            cblas_sscal(nfftsteps, norm_factor, &moltypes_dos_raw_trn_z[h*nfftsteps], 1);
            cblas_sscal(nfftsteps, norm_factor, &moltypes_dos_raw_rot[h*nfftsteps], 1);
            cblas_sscal(nfftsteps, norm_factor, &moltypes_dos_raw_rot_a[h*nfftsteps], 1);
            cblas_sscal(nfftsteps, norm_factor, &moltypes_dos_raw_rot_b[h*nfftsteps], 1);
            cblas_sscal(nfftsteps, norm_factor, &moltypes_dos_raw_rot_c[h*nfftsteps], 1);
            cblas_sscal(nfftsteps, norm_factor, &moltypes_dos_raw_vib[h*nfftsteps], 1);
            cblas_sscal(nfftsteps, norm_factor, &moltypes_dos_raw_vib_x[h*nfftsteps], 1);
            cblas_sscal(nfftsteps, norm_factor, &moltypes_dos_raw_vib_y[h*nfftsteps], 1);
            cblas_sscal(nfftsteps, norm_factor, &moltypes_dos_raw_vib_z[h*nfftsteps], 1);
            cblas_sscal(nfftsteps, norm_factor, &moltypes_dos_raw_rot_alt[h*nfftsteps], 1);
            cblas_sscal(nfftsteps, norm_factor, &moltypes_dos_raw_rot_alt_x[h*nfftsteps], 1);
            cblas_sscal(nfftsteps, norm_factor, &moltypes_dos_raw_rot_alt_y[h*nfftsteps], 1);
            cblas_sscal(nfftsteps, norm_factor, &moltypes_dos_raw_rot_alt_z[h*nfftsteps], 1);
            cblas_sscal(nfftsteps, norm_factor, &moltypes_dos_raw_x_trn_rot[h*nfftsteps], 1);
            cblas_sscal(nfftsteps, norm_factor, &moltypes_dos_raw_x_trn_vib[h*nfftsteps], 1);
            cblas_sscal(nfftsteps, norm_factor, &moltypes_dos_raw_x_rot_vib[h*nfftsteps], 1);
        }

        // moments of inertia of moltype
        for (int i=0; i<nmols; i++)
        {
            moltypes_moments_of_inertia[3*mols_moltypenr[i]+0] += mol_moments_of_inertia[3*i+0];
            moltypes_moments_of_inertia[3*mols_moltypenr[i]+1] += mol_moments_of_inertia[3*i+1];
            moltypes_moments_of_inertia[3*mols_moltypenr[i]+2] += mol_moments_of_inertia[3*i+2];
        }
        for (int h=0; h<nmoltypes; h++)
        {
            cblas_sscal(3, 1.0 / (float) moltypes_nmols[h], &moltypes_moments_of_inertia[3*h], 1);
        }

        // translational spectrum
        char filename[255] = {0};
        sprintf(filename, "sample%d_dos_trn.txt", sample);
        f = fopen(filename, "w");

        for (int h=0; h<nmoltypes; h++)
        {
            for (int t=0; t<nfftsteps; t++)
            {
                if (t!=0) fprintf(f, " ");
                fprintf(f, "%f", moltypes_dos_raw_trn[h*nfftsteps + t]);
            }
            fprintf(f, "\n");
        }
        fclose(f);

        // individual axis (xyz) translational spectrum
        if (arguments.calc_components)
        {
            char filenamea[255] = {0};
            char filenameb[255] = {0};
            char filenamec[255] = {0};
            sprintf(filenamea, "sample%d_dos_trn_x.txt", sample);
            sprintf(filenameb, "sample%d_dos_trn_y.txt", sample);
            sprintf(filenamec, "sample%d_dos_trn_z.txt", sample);
            fa = fopen(filenamea, "w");
            fb = fopen(filenameb, "w");
            fc = fopen(filenamec, "w");
            for (int h=0; h<nmoltypes; h++)
            {
                for (int t=0; t<nfftsteps; t++)
                {
                    if (t!=0) fprintf(fa, " ");
                    if (t!=0) fprintf(fb, " ");
                    if (t!=0) fprintf(fc, " ");
                    fprintf(fa, "%f", moltypes_dos_raw_trn_x[h*nfftsteps + t]);
                    fprintf(fb, "%f", moltypes_dos_raw_trn_y[h*nfftsteps + t]);
                    fprintf(fc, "%f", moltypes_dos_raw_trn_z[h*nfftsteps + t]);
                }
                fprintf(fa, "\n");
                fprintf(fb, "\n");
                fprintf(fc, "\n");
            }
            fclose(fa);
            fclose(fb);
            fclose(fc);
        }

        // rotational spectrum
        sprintf(filename, "sample%d_dos_rot.txt", sample);
        f = fopen(filename, "w");

        for (int h=0; h<nmoltypes; h++)
        {
            for (int t=0; t<nfftsteps; t++)
            {
                if (t!=0) fprintf(f, " ");
                fprintf(f, "%f", moltypes_dos_raw_rot[h*nfftsteps + t]);
            }
            fprintf(f, "\n");
        }
        fclose(f);

        // individual axis (abc) rotational spectrum
        if (arguments.calc_components)
        {
            char filenamea[255] = {0};
            char filenameb[255] = {0};
            char filenamec[255] = {0};
            sprintf(filenamea, "sample%d_dos_rot_a.txt", sample);
            sprintf(filenameb, "sample%d_dos_rot_b.txt", sample);
            sprintf(filenamec, "sample%d_dos_rot_c.txt", sample);
            fa = fopen(filenamea, "w");
            fb = fopen(filenameb, "w");
            fc = fopen(filenamec, "w");
            for (int h=0; h<nmoltypes; h++)
            {
                for (int t=0; t<nfftsteps; t++)
                {
                    if (t!=0) fprintf(fa, " ");
                    if (t!=0) fprintf(fb, " ");
                    if (t!=0) fprintf(fc, " ");
                    fprintf(fa, "%f", moltypes_dos_raw_rot_a[h*nfftsteps + t]);
                    fprintf(fb, "%f", moltypes_dos_raw_rot_b[h*nfftsteps + t]);
                    fprintf(fc, "%f", moltypes_dos_raw_rot_c[h*nfftsteps + t]);
                }
                fprintf(fa, "\n");
                fprintf(fb, "\n");
                fprintf(fc, "\n");
            }
            fclose(fa);
            fclose(fb);
            fclose(fc);
        }

        // vibrational spectrum
        sprintf(filename, "sample%d_dos_vib.txt", sample);
        f = fopen(filename, "w");
        for (int h=0; h<nmoltypes; h++)
        {
            for (int t=0; t<nfftsteps; t++)
            {
                if (t!=0) fprintf(f, " ");
                fprintf(f, "%f", moltypes_dos_raw_vib[h*nfftsteps + t]);
            }
            fprintf(f, "\n");
        }
        fclose(f);

        // individual axis (xyz) vibrational spectrum
        if (arguments.calc_components)
        {
            char filenamea[255] = {0};
            char filenameb[255] = {0};
            char filenamec[255] = {0};
            sprintf(filenamea, "sample%d_dos_vib_x.txt", sample);
            sprintf(filenameb, "sample%d_dos_vib_y.txt", sample);
            sprintf(filenamec, "sample%d_dos_vib_z.txt", sample);
            fa = fopen(filenamea, "w");
            fb = fopen(filenameb, "w");
            fc = fopen(filenamec, "w");
            for (int h=0; h<nmoltypes; h++)
            {
                for (int t=0; t<nfftsteps; t++)
                {
                    if (t!=0) fprintf(fa, " ");
                    if (t!=0) fprintf(fb, " ");
                    if (t!=0) fprintf(fc, " ");
                    fprintf(fa, "%f", moltypes_dos_raw_vib_x[h*nfftsteps + t]);
                    fprintf(fb, "%f", moltypes_dos_raw_vib_y[h*nfftsteps + t]);
                    fprintf(fc, "%f", moltypes_dos_raw_vib_z[h*nfftsteps + t]);
                }
                fprintf(fa, "\n");
                fprintf(fb, "\n");
                fprintf(fc, "\n");
            }
            fclose(fa);
            fclose(fb);
            fclose(fc);
        }

        // rotational alternative spectrum
        if (arguments.calc_rot_alt)
        {
            sprintf(filename, "sample%d_dos_rotalt.txt", sample);
            f = fopen(filename, "w");
            for (int h=0; h<nmoltypes; h++)
            {
                for (int t=0; t<nfftsteps; t++)
                {
                    if (t!=0) fprintf(f, " ");
                    fprintf(f, "%f", moltypes_dos_raw_rot_alt[h*nfftsteps + t]);
                }
                fprintf(f, "\n");
            }
            fclose(f);

            // individual axis (xyz) rotational alternative spectrum
            if (arguments.calc_components)
            {
                char filenamea[255] = {0};
                char filenameb[255] = {0};
                char filenamec[255] = {0};
                sprintf(filenamea, "sample%d_dos_rotalt_x.txt", sample);
                sprintf(filenameb, "sample%d_dos_rotalt_y.txt", sample);
                sprintf(filenamec, "sample%d_dos_rotalt_z.txt", sample);
                fa = fopen(filenamea, "w");
                fb = fopen(filenameb, "w");
                fc = fopen(filenamec, "w");
                for (int h=0; h<nmoltypes; h++)
                {
                    for (int t=0; t<nfftsteps; t++)
                    {
                        if (t!=0) fprintf(fa, " ");
                        if (t!=0) fprintf(fb, " ");
                        if (t!=0) fprintf(fc, " ");
                        fprintf(fa, "%f", moltypes_dos_raw_rot_alt_x[h*nfftsteps + t]);
                        fprintf(fb, "%f", moltypes_dos_raw_rot_alt_y[h*nfftsteps + t]);
                        fprintf(fc, "%f", moltypes_dos_raw_rot_alt_z[h*nfftsteps + t]);
                    }
                    fprintf(fa, "\n");
                    fprintf(fb, "\n");
                    fprintf(fc, "\n");
                }
                fclose(fa);
                fclose(fb);
                fclose(fc);
            }
        }

        // cross spectra
        if (arguments.calc_cross)
        {
            char filenametrnrot[255] = {0};
            char filenametrnvib[255] = {0};
            char filenamerotvib[255] = {0};
            sprintf(filenametrnrot, "sample%d_dos_x_trn_rot.txt", sample);
            sprintf(filenametrnvib, "sample%d_dos_x_trn_vib.txt", sample);
            sprintf(filenamerotvib, "sample%d_dos_x_rot_vib.txt", sample);
            ftr = fopen(filenametrnrot, "w");
            ftv = fopen(filenametrnvib, "w");
            frv = fopen(filenamerotvib, "w");
            for (int h=0; h<nmoltypes; h++)
            {
                for (int t=0; t<nfftsteps; t++)
                {
                    if (t!=0) fprintf(ftr, " ");
                    if (t!=0) fprintf(ftv, " ");
                    if (t!=0) fprintf(frv, " ");
                    fprintf(ftr, "%f", moltypes_dos_raw_x_trn_rot[h*nfftsteps + t]);
                    fprintf(ftv, "%f", moltypes_dos_raw_x_trn_vib[h*nfftsteps + t]);
                    fprintf(frv, "%f", moltypes_dos_raw_x_rot_vib[h*nfftsteps + t]);
                }
                fprintf(ftr, "\n");
                fprintf(ftv, "\n");
                fprintf(frv, "\n");
            }
            fclose(ftr);
            fclose(ftv);
            fclose(frv);
        }

        // moments of inertia
        sprintf(filename, "sample%d_moments_of_inertia.txt", sample);
        f = fopen(filename, "w");
        for (int h=0; h<nmoltypes; h++)
        {
            fprintf(f, "%f %f %f\n", moltypes_moments_of_inertia[3*h+0],
                    moltypes_moments_of_inertia[3*h+1], moltypes_moments_of_inertia[3*h+2]);
        }
        fclose(f);

        // free DoSes
        free(mol_moments_of_inertia);
        free(moltypes_moments_of_inertia);
        free(moltypes_dos_raw_trn);
        free(moltypes_dos_raw_trn_x);
        free(moltypes_dos_raw_trn_y);
        free(moltypes_dos_raw_trn_z);
        free(moltypes_dos_raw_rot);
        free(moltypes_dos_raw_rot_a);
        free(moltypes_dos_raw_rot_b);
        free(moltypes_dos_raw_rot_c);
        free(moltypes_dos_raw_rot_alt);
        free(moltypes_dos_raw_rot_alt_x);
        free(moltypes_dos_raw_rot_alt_y);
        free(moltypes_dos_raw_rot_alt_z);
        free(moltypes_dos_raw_vib);
        free(moltypes_dos_raw_vib_x);
        free(moltypes_dos_raw_vib_y);
        free(moltypes_dos_raw_vib_z);
        free(moltypes_dos_raw_x_trn_rot);
        free(moltypes_dos_raw_x_trn_vib);
        free(moltypes_dos_raw_x_rot_vib);
    }
    verbPrintf(verbosity, "finished all samples\n");

    verbPrintf(verbosity, "generating frequencies.txt\n");
    f = fopen("frequencies.txt", "w");
    for (int i=0; i<nfftsteps; i++)
    {
        if (i!=0) fprintf(f, " ");
        fprintf(f, "%f", i / (framelength * nblocksteps));
    }
    fclose(f);

    // free arrays
    free(moltypes_nmols);
    free(moltypes_natomspermol);

    for (int h=0; h<nmoltypes; h++)
    {
        free(moltypes_atommasses[h]);
        free(moltypes_abc_indicators[h]);
    }
    free(moltypes_atommasses);
    free(moltypes_abc_indicators);

    free(moltypes_rot_treat);
    free(moltypes_firstmol);
    free(moltypes_firstatom);
    free(mols_moltypenr);
    free(mols_natoms);
    free(mols_mass);
    free(mols_firstatom);

    return 0;
}
