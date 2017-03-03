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
const char* argp_program_version = "dos-calc develop";
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
    bool verbosity, dump_vel;
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
    int natoms;
    int nmols;
    int nmoltypes;

    // file pointers for output
    FILE* f;
    FILE* fa;
    FILE* fb;
    FILE* fc;

    // results of decomposeVelocities and DOSCalculation
    int result = 0;
    int result2 = 0;

    // start scanning
    scanf("%d", &nsamples);
    scanf("%d", &nblocks);
    scanf("%ld", &nblocksteps);

    int nfftsteps = nblocksteps / 2 + 1;

    scanf("%d", &natoms);
    scanf("%d", &nmols);
    scanf("%d", &nmoltypes);

    verbPrintf(verbosity, "Integers scanned:\n");
    verbPrintf(verbosity, "%d\n", nsamples);
    verbPrintf(verbosity, "%d\n", nblocks);
    verbPrintf(verbosity, "%ld\n", nblocksteps);
    verbPrintf(verbosity, "%d\n", natoms);
    verbPrintf(verbosity, "%d\n", nmols);
    verbPrintf(verbosity, "%d\n", nmoltypes);

    int* moltype_firstmol = calloc(nmoltypes, sizeof(int));
    for (int h=0; h<nmoltypes; h++)
    {
        scanf("%d", &moltype_firstmol[h]);
    }

    int* moltype_firstatom = calloc(nmoltypes, sizeof(int));
    for (int h=0; h<nmoltypes; h++)
    {
        scanf("%d", &moltype_firstatom[h]);
    }

    int* moltype_nmols = calloc(nmoltypes, sizeof(int));
    for (int h=0; h<nmoltypes; h++)
    {
        scanf("%d", &moltype_nmols[h]);
    }

    int* moltype_natomtypes = calloc(nmoltypes, sizeof(int));
    for (int h=0; h<nmoltypes; h++)
    {
        scanf("%d", &moltype_natomtypes[h]);
    }

    int* moltype_abc_indicators = calloc(nmoltypes * 4, sizeof(int));
    for (int h=0; h<nmoltypes; h++)
    {
        for (int i=0; i<4; i++)
        {
            scanf("%d", &moltype_abc_indicators[4*h + i]);
        }
    }

    int* mol_firstatom = calloc(nmols, sizeof(int));
    for (int i=0; i<nmols; i++)
    {
        scanf("%d", &mol_firstatom[i]);
    }

    int* mol_natoms = calloc(nmols, sizeof(int));
    for (int i=0; i<nmols; i++)
    {
        scanf("%d", &mol_natoms[i]);
    }

    float* mol_mass = calloc(nmols, sizeof(float));
    for (int i=0; i<nmols; i++)
    {
        scanf("%f", &mol_mass[i]);
    }

    int* mol_moltypenr = calloc(nmols, sizeof(int));
    for (int i=0; i<nmols; i++)
    {
        scanf("%d", &mol_moltypenr[i]);
    }

    float* atom_mass = calloc(natoms, sizeof(float));
    for (int j=0; j<natoms; j++)
    {
        scanf("%f", &atom_mass[j]);
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
        float* moltype_dos_raw_trn = calloc(nmoltypes*nfftsteps, sizeof(float));
        float* moltype_dos_raw_rot = calloc(nmoltypes*nfftsteps, sizeof(float));
        float* moltype_dos_raw_rot_a = calloc(nmoltypes*nfftsteps, sizeof(float));
        float* moltype_dos_raw_rot_b = calloc(nmoltypes*nfftsteps, sizeof(float));
        float* moltype_dos_raw_rot_c = calloc(nmoltypes*nfftsteps, sizeof(float));
        float* moltype_dos_raw_vib = calloc(nmoltypes*nfftsteps, sizeof(float));

        verbPrintf(verbosity, "going through %d blocks\n", nblocks);
        for (int block=0; block<nblocks; block++)
        {
            verbPrintf(verbosity, "now doing block %d\n", block);
            verbPrintf(verbosity, "start decomposition\n");

            float* mol_velocities_trn = calloc(nmols*3*nblocksteps, sizeof(float));
            float* omegas_sqrt_i = calloc(nmols*3*nblocksteps, sizeof(float));
            float* velocities_vib = calloc(natoms*3*nblocksteps, sizeof(float));

            result = decomposeVelocities (trj_in,
                    header,
                    nblocksteps,
                    natoms, 
                    nmols, 
                    nmoltypes,
                    mol_firstatom,
                    mol_natoms,
                    mol_moltypenr,
                    atom_mass,
                    mol_mass,
                    moltype_natomtypes,
                    moltype_abc_indicators,
                    mol_velocities_trn,
                    omegas_sqrt_i,
                    velocities_vib,
                    mol_moments_of_inertia);

            // dump first block only!
            if (arguments.dump_vel == 1 && block == 0)
            {
                verbPrintf(verbosity, "start dumping (first block only)\n");
                f = fopen("mol_omega_sqrt_i.txt", "w");
                for (int h=0; h<nmoltypes; h++)
                {
                    int first_dof = moltype_firstmol[h]*3;
                    int last_dof = moltype_firstmol[h]*3 + moltype_nmols[h]*3;
                    for (int i=first_dof; i<last_dof; i++)
                    {
                        for (int t=0; t<nblocksteps; t++)
                        {
                            if (t!=0) fprintf(f, " ");
                            fprintf(f, "%f", omegas_sqrt_i[i*nblocksteps + t]);
                        }
                        fprintf(f, "\n");
                    }
                }
                fclose(f);

                f = fopen("mol_velocities.txt", "w");
                for (int h=0; h<nmoltypes; h++)
                {
                    int first_dof = moltype_firstmol[h]*3;
                    int last_dof = moltype_firstmol[h]*3 + moltype_nmols[h]*3;
                    for (int i=first_dof; i<last_dof; i++)
                    {
                        for (int t=0; t<nblocksteps; t++)
                        {
                            if (t!=0) fprintf(f, " ");
                            fprintf(f, "%f", mol_velocities_trn[i*nblocksteps + t]);
                        }
                        fprintf(f, "\n");
                    }
                }
                fclose(f);

                f = fopen("atom_velocities_vib.txt", "w");
                for (int h=0; h<nmoltypes; h++)
                {
                    int first_dof_vib = moltype_firstatom[h]*3;
                    int last_dof_vib = moltype_firstatom[h]*3 + moltype_nmols[h]*moltype_natomtypes[h]*3;
                    for (int i=first_dof_vib; i<last_dof_vib; i++)
                    {
                        for (int t=0; t<nblocksteps; t++)
                        {
                            if (t!=0) fprintf(f, " ");
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
                    moltype_firstmol,
                    moltype_firstatom,
                    moltype_nmols,
                    moltype_natomtypes,
                    atom_mass,
                    mol_velocities_trn,
                    omegas_sqrt_i,
                    velocities_vib,
                    moltype_dos_raw_trn,
                    moltype_dos_raw_rot,
                    moltype_dos_raw_rot_a,
                    moltype_dos_raw_rot_b,
                    moltype_dos_raw_rot_c,
                    moltype_dos_raw_vib);

            //free velocities
            free(mol_velocities_trn);
            free(omegas_sqrt_i);
            free(velocities_vib);
        } 
        verbPrintf(verbosity, "finished all blocks\n");

        // divide results by number of blocks
        cblas_sscal(nmols*3, 1.0 / (float)nblocks, mol_moments_of_inertia, 1);
        cblas_sscal(nmoltypes*nfftsteps, 1.0 / (float)nblocks, moltype_dos_raw_trn, 1);
        cblas_sscal(nmoltypes*nfftsteps, 1.0 / (float)nblocks, moltype_dos_raw_rot, 1);
        cblas_sscal(nmoltypes*nfftsteps, 1.0 / (float)nblocks, moltype_dos_raw_rot_a, 1);
        cblas_sscal(nmoltypes*nfftsteps, 1.0 / (float)nblocks, moltype_dos_raw_rot_b, 1);
        cblas_sscal(nmoltypes*nfftsteps, 1.0 / (float)nblocks, moltype_dos_raw_rot_c, 1);
        cblas_sscal(nmoltypes*nfftsteps, 1.0 / (float)nblocks, moltype_dos_raw_vib, 1);

        char filename[255] = {0};
        sprintf(filename, "sample%d_dos_trn.txt", sample);
        f = fopen(filename, "w");

        for (int h=0; h<nmoltypes; h++)
        {
            for (int t=0; t<nfftsteps; t++)
            {
                if (t!=0) fprintf(f, " ");
                fprintf(f, "%f", moltype_dos_raw_trn[h*nfftsteps + t]);
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
                fprintf(f, "%f", moltype_dos_raw_rot[h*nfftsteps + t]);
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
                fprintf(fa, "%f", moltype_dos_raw_rot_a[h*nfftsteps + t]);
                fprintf(fb, "%f", moltype_dos_raw_rot_b[h*nfftsteps + t]);
                fprintf(fc, "%f", moltype_dos_raw_rot_c[h*nfftsteps + t]);
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
                fprintf(f, "%f", moltype_dos_raw_vib[h*nfftsteps + t]);
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
        free(moltype_dos_raw_trn);
        free(moltype_dos_raw_rot);
        free(moltype_dos_raw_rot_a);
        free(moltype_dos_raw_rot_b);
        free(moltype_dos_raw_rot_c);
        free(moltype_dos_raw_vib);

    }
    verbPrintf(verbosity, "finished all samples\n");

    // free input arrays
    free(moltype_firstmol);
    free(moltype_firstatom);
    free(moltype_nmols);
    free(moltype_natomtypes);
    free(moltype_abc_indicators);
    free(mol_firstatom);
    free(mol_natoms);
    free(mol_mass);
    free(mol_moltypenr);
    free(atom_mass);

    return result + result2;
}
