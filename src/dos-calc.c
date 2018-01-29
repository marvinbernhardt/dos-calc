#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdbool.h>
#include <argp.h>
#include <cblas.h>
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/trrio.h"
#include "velocity-decomposition.c"
#include "fft.c"
#include "verbPrintf.c"

// CLI stuff
const char* argp_program_version = VERSION;
const char* argp_program_bug_address = "<bernhardt@cpc.tu-darmstadt.de>";
static char doc[] = "dos-calc -- a programm to calculate densities of states from trajectories";
static char args_doc[] = "";

static struct argp_option options[] = {
    {"dump",     'd', 0,      0,  "Dump velocities of first block" },
    {"verbose",  'v', 0,      0,  "Produce verbose output" },
    {"file",     'f', "FILE", 0,  "Input .trr trajectory file (default: traj.trr)"},
    { 0 }
};

struct arguments
{
    bool verbosity;
    bool dump_vel;
    char *file;
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
        case 'f':
            arguments->file = arg;
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
    arguments.file = "traj.trr";

    // parse command line arguments
    argp_parse(&argp, argc, argv, 0, 0, &arguments);
    bool verbosity = arguments.verbosity;

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

    // results of decomposeVelocities and DOSCalculation
    int result = 0;
    int result2 = 0;

    // start scanning
    verbPrintf(verbosity, "Start reading user input\n");
    verbPrintf(verbosity, "Reading integers\n");

    verbPrintf(verbosity, "nsamples: ");
    scanf("%d", &nsamples);
    verbPrintf(verbosity, "%d\n", nsamples);

    verbPrintf(verbosity, "nblocks: ");
    scanf("%d", &nblocks);
    verbPrintf(verbosity, "%d\n", nblocks);

    verbPrintf(verbosity, "nblocksteps: ");
    scanf("%ld", &nblocksteps);
    verbPrintf(verbosity, "%ld\n", nblocksteps);
    int nfftsteps = nblocksteps / 2 + 1;

    verbPrintf(verbosity, "nmoltypes: ");
    scanf("%d", &nmoltypes);
    verbPrintf(verbosity, "%d\n", nmoltypes);

    // scanning moltypes
    int* moltypes_nmols = calloc(nmoltypes, sizeof(int));
    int* moltypes_natomspermol = calloc(nmoltypes, sizeof(int));
    float** moltypes_atommasses = (float**) malloc(nmoltypes * sizeof(float *));
    char* moltypes_rot_treat = calloc(nmoltypes, sizeof(char));
    int** moltypes_abc_indicators = (int**) malloc(nmoltypes * sizeof(int *));

    for (int h=0; h<nmoltypes; h++)
    {
        verbPrintf(verbosity, "now reading moltype %d of %d.\n", h+1, nmoltypes);

        verbPrintf(verbosity, "moltype %d nmols: ", h+1);
        scanf("%d", &moltypes_nmols[h]);
        verbPrintf(verbosity, "%d\n", moltypes_nmols[h]);

        verbPrintf(verbosity, "moltype %d natomspermol: ", h+1);
        scanf("%d", &moltypes_natomspermol[h]);
        verbPrintf(verbosity, "%d\n", moltypes_natomspermol[h]);

        verbPrintf(verbosity, "moltype %d atommasses: ", h+1);
        moltypes_atommasses[h] = (float*) malloc(moltypes_natomspermol[h] * sizeof(float));
        for (int j=0; j<moltypes_natomspermol[h]; j++)
        {
            scanf("%f", &moltypes_atommasses[h][j]);
        }
        for (int j=0; j<moltypes_natomspermol[h]; j++)
        {
            verbPrintf(verbosity, "%f ", moltypes_atommasses[h][j]);
        }
        verbPrintf(verbosity, "\n");

        verbPrintf(verbosity, "moltype %d rot_treat: ", h+1);
        scanf(" %c", &moltypes_rot_treat[h]);
        verbPrintf(verbosity, "%c\n", moltypes_rot_treat[h]);

        verbPrintf(verbosity, "moltype %d abc_indicators: ", h+1);
        moltypes_abc_indicators[h] = (int*) malloc(4 * sizeof(int));
        for (int j=0; j<4; j++)
        {
            scanf("%d", &moltypes_abc_indicators[h][j]);
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


    // open file
    verbPrintf(verbosity, "starting with file %s\n", arguments.file);
    t_fileio* trj_in = gmx_trr_open(arguments.file, "r");

    // read header
    verbPrintf(verbosity, "start reading header\n");
    gmx_trr_header_t header;
    gmx_bool bOK;
    gmx_trr_read_frame_header(trj_in, &header, &bOK);
    if (!bOK) {
        printf("trajectory header broken\n");
        return 1;
    }

    // go back to start of file
    gmx_fio_rewind(trj_in);

    verbPrintf(verbosity, "going through %d samples\n", nsamples);
    for (int sample=0; sample<nsamples; sample++)
    {
        verbPrintf(verbosity, "now doing sample %d\n", sample);

        // output arrays
        float* mol_moments_of_inertia = calloc(nmols*3, sizeof(float));
        float* moltypes_dos_raw_trn = calloc(nmoltypes*nfftsteps, sizeof(float));
        float* moltypes_dos_raw_rot = calloc(nmoltypes*nfftsteps, sizeof(float));
        float* moltypes_dos_raw_rot_a = calloc(nmoltypes*nfftsteps, sizeof(float));
        float* moltypes_dos_raw_rot_b = calloc(nmoltypes*nfftsteps, sizeof(float));
        float* moltypes_dos_raw_rot_c = calloc(nmoltypes*nfftsteps, sizeof(float));
        float* moltypes_dos_raw_vib = calloc(nmoltypes*nfftsteps, sizeof(float));

        verbPrintf(verbosity, "going through %d blocks\n", nblocks);
        for (int block=0; block<nblocks; block++)
        {
            verbPrintf(verbosity, "now doing block %d\n", block);
            verbPrintf(verbosity, "start decomposition\n");

            float* mol_velocities_sqrt_m_trn = calloc(nmols*3*nblocksteps, sizeof(float));
            float* mol_omegas_sqrt_i_rot = calloc(nmols*3*nblocksteps, sizeof(float));
            float* atom_velocities_sqrt_m_vib = calloc(natoms*3*nblocksteps, sizeof(float));

            result = decomposeVelocities (trj_in,
                    header,
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
                    mol_velocities_sqrt_m_trn,
                    mol_omegas_sqrt_i_rot,
                    atom_velocities_sqrt_m_vib,
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
            }

            verbPrintf(verbosity, "start DoS calculation (FFT)\n");
            result2 = DOSCalculation (nmoltypes,
                    nblocksteps,
                    nfftsteps,
                    moltypes_firstmol,
                    moltypes_firstatom,
                    moltypes_nmols,
                    moltypes_natomspermol,
                    mol_velocities_sqrt_m_trn,
                    mol_omegas_sqrt_i_rot,
                    atom_velocities_sqrt_m_vib,
                    moltypes_dos_raw_trn, //output
                    moltypes_dos_raw_rot,
                    moltypes_dos_raw_rot_a,
                    moltypes_dos_raw_rot_b,
                    moltypes_dos_raw_rot_c,
                    moltypes_dos_raw_vib);

            //free velocities
            free(mol_velocities_sqrt_m_trn);
            free(mol_omegas_sqrt_i_rot);
            free(atom_velocities_sqrt_m_vib);
        }
        verbPrintf(verbosity, "finished all blocks\n");

        // divide results by number of blocks (and number of blocksteps for the moi)
        cblas_sscal(nmols*3, 1.0 / (float) nblocks / (float) nblocksteps, mol_moments_of_inertia, 1);
        cblas_sscal(nmoltypes*nfftsteps, 1.0 / (float)nblocks, moltypes_dos_raw_trn, 1);
        cblas_sscal(nmoltypes*nfftsteps, 1.0 / (float)nblocks, moltypes_dos_raw_rot, 1);
        cblas_sscal(nmoltypes*nfftsteps, 1.0 / (float)nblocks, moltypes_dos_raw_rot_a, 1);
        cblas_sscal(nmoltypes*nfftsteps, 1.0 / (float)nblocks, moltypes_dos_raw_rot_b, 1);
        cblas_sscal(nmoltypes*nfftsteps, 1.0 / (float)nblocks, moltypes_dos_raw_rot_c, 1);
        cblas_sscal(nmoltypes*nfftsteps, 1.0 / (float)nblocks, moltypes_dos_raw_vib, 1);

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

        sprintf(filename, "sample%d_moi.txt", sample);
        f = fopen(filename, "w");
        for (int i=0; i<nmols; i++)
        {
            fprintf(f, "%f %f %f\n", mol_moments_of_inertia[3*i+0],
                    mol_moments_of_inertia[3*i+1], mol_moments_of_inertia[3*i+2]);
        }
        fclose(f);

        // free DoSes
        free(mol_moments_of_inertia);
        free(moltypes_dos_raw_trn);
        free(moltypes_dos_raw_rot);
        free(moltypes_dos_raw_rot_a);
        free(moltypes_dos_raw_rot_b);
        free(moltypes_dos_raw_rot_c);
        free(moltypes_dos_raw_vib);

    }
    verbPrintf(verbosity, "finished all samples\n");

    // free input arrays
    free(moltypes_firstmol);
    free(moltypes_firstatom);
    free(moltypes_nmols);
    free(moltypes_natomspermol);
    free(moltypes_abc_indicators);
    free(mols_firstatom);
    free(mols_natoms);
    free(mols_mass);
    free(mols_moltypenr);

    return result + result2;
}
