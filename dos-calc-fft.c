#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <cblas.h>
#include <fftw3.h>

//#define DEBUG

#ifdef DEBUG
#define DPRINT(...) do{ fprintf( stderr, __VA_ARGS__ ); } while( 0 )
#else
#define DPRINT(...)
#endif

//#define VERBOSE

#ifdef VERBOSE
#define VPRINT(...) do{ fprintf( stderr, __VA_ARGS__ ); } while( 0 )
#else
#define VPRINT(...)
#endif

int DOSCalculation (int nmoltypes,
        long ntrajsteps,
        long nfftsteps,
        int* moltype_firstmol,
        int* moltype_firstatom,
        int* moltype_nmols,
        int* moltype_natomtypes,
        float* atom_mass,
        float* mol_velocities_trn,
        float* omegas_sqrt_i,
        float* velocities_vib,
        float* moltype_dos_raw_trn,
        float* moltype_dos_raw_rot,
        float* moltype_dos_raw_rot_a,
        float* moltype_dos_raw_rot_b,
        float* moltype_dos_raw_rot_c,
        float* moltype_dos_raw_vib) 
{
    VPRINT("starting with FFT of translational and rotational dofs\n");

    for (int h=0; h<nmoltypes; h++)
    {
        DPRINT("moltype nr %d\n", h);
        float* fft_in = calloc(ntrajsteps, sizeof(float));
        fftwf_complex* fft_out = fftwf_malloc(sizeof(fftwf_complex) * nfftsteps);
        float* fft_out_squared = calloc(nfftsteps, sizeof(float));
        DPRINT("creating plan\n");
        fftwf_plan p = fftwf_plan_dft_r2c_1d(ntrajsteps, fft_in, fft_out, FFTW_MEASURE);

        int first_dof = moltype_firstmol[h]*3;
        int last_dof = moltype_firstmol[h]*3 + moltype_nmols[h]*3;
        for (int i=first_dof; i<last_dof; i++)
        {
            DPRINT("translational dof nr %d\n", i);
            memcpy(fft_in, &mol_velocities_trn[i*ntrajsteps], ntrajsteps * sizeof(float) );
            fftwf_execute(p);

            DPRINT("abs and square of fft\n");
            for (int t=0; t<nfftsteps; t++)
            {
                fft_out_squared[t] = cabs(fft_out[t] * fft_out[t]);
            }

            cblas_saxpy(nfftsteps, 1.0, fft_out_squared, 1, 
                    &moltype_dos_raw_trn[h*nfftsteps], 1);

            DPRINT("rotational dof nr %d\n", i);
            memcpy(fft_in, &omegas_sqrt_i[i*ntrajsteps], ntrajsteps * sizeof(float) );
            fftwf_execute(p);

            DPRINT("abs and square of fft\n");
            for (int t=0; t<nfftsteps; t++)
            {
                fft_out_squared[t] = cabs(fft_out[t] * fft_out[t]);
            }

            cblas_saxpy(nfftsteps, 1.0, fft_out_squared, 1, 
                    &moltype_dos_raw_rot[h*nfftsteps], 1);
            if ( i % 3 == 0 ) cblas_saxpy(nfftsteps, 1.0, fft_out_squared, 1, &moltype_dos_raw_rot_a[h*nfftsteps], 1);
            if ( i % 3 == 1 ) cblas_saxpy(nfftsteps, 1.0, fft_out_squared, 1, &moltype_dos_raw_rot_b[h*nfftsteps], 1);
            if ( i % 3 == 2 ) cblas_saxpy(nfftsteps, 1.0, fft_out_squared, 1, &moltype_dos_raw_rot_c[h*nfftsteps], 1);
        }


        // vibration
        int first_dof_vib = moltype_firstatom[h]*3;
        int last_dof_vib = moltype_firstatom[h]*3 + moltype_nmols[h]*moltype_natomtypes[h]*3;
        for (int i=first_dof_vib; i<last_dof_vib; i++)
        {
            DPRINT("vibrational dof nr %d\n", i);
            memcpy(fft_in, &velocities_vib[i*ntrajsteps], ntrajsteps * sizeof(float) );
            fftwf_execute(p);

            DPRINT("abs and square of fft\n");
            for (int t=0; t<nfftsteps; t++)
            {
                fft_out_squared[t] = atom_mass[i/3] * cabs(fft_out[t] * fft_out[t]);
            }

            cblas_saxpy(nfftsteps, 1.0, fft_out_squared, 1, 
                    &moltype_dos_raw_vib[h*nfftsteps], 1);
        }

        DPRINT("moltype done\n");

        fftwf_destroy_plan(p);
        free(fft_in);
        fftwf_free(fft_out);
        free(fft_out_squared);
    }

    return 0;
}
