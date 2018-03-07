#include <stdlib.h>
#include <complex.h>
#include <stdbool.h>
#include <math.h>
#include <cblas.h>
#include <fftw3.h>

//#define DEBUG

#ifdef DEBUG
#define DPRINT(...) do{ fprintf( stdout, __VA_ARGS__ ); } while( 0 )
#else
#define DPRINT(...)
#endif

//#define VERBOSE

#ifdef VERBOSE
#define VPRINT(...) do{ fprintf( stdout, __VA_ARGS__ ); } while( 0 )
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
        float* mol_velocities_sqrt_m_trn,
        float* mol_omegas_sqrt_i_rot,
        float* atom_velocities_sqrt_m_vib,
        bool calc_cross,
        float* moltype_dos_raw_trn,  // output
        float* moltype_dos_raw_rot,
        float* moltype_dos_raw_rot_a,
        float* moltype_dos_raw_rot_b,
        float* moltype_dos_raw_rot_c,
        float* moltype_dos_raw_vib,
        float* moltype_dos_raw_x_trn_rot_a,
        float* moltype_dos_raw_x_trn_rot_b,
        float* moltype_dos_raw_x_trn_rot_c)
{
    VPRINT("starting with FFT of translational and rotational dofs\n");

    for (int h=0; h<nmoltypes; h++)
    {
        DPRINT("moltype nr %d\n", h);
        float* fft_in = calloc(ntrajsteps, sizeof(float));
        fftwf_complex* fft_out_trn = fftwf_malloc(sizeof(fftwf_complex) * nfftsteps);
        fftwf_complex* fft_out_rot = fftwf_malloc(sizeof(fftwf_complex) * nfftsteps);
        float* fft_out_squared = calloc(nfftsteps, sizeof(float));
        DPRINT("creating plan translation\n");
        fftwf_plan plan_trn = fftwf_plan_dft_r2c_1d(ntrajsteps, fft_in, fft_out_trn, FFTW_MEASURE);
        DPRINT("creating plan rotational\n");
        fftwf_plan plan_rot = fftwf_plan_dft_r2c_1d(ntrajsteps, fft_in, fft_out_rot, FFTW_MEASURE);

        int first_dof = moltype_firstmol[h]*3;
        int last_dof = moltype_firstmol[h]*3 + moltype_nmols[h]*3;
        for (int i=first_dof; i<last_dof; i++)
        {
            DPRINT("translational dof nr %d\n", i);
            memcpy(fft_in, &mol_velocities_sqrt_m_trn[i*ntrajsteps], ntrajsteps * sizeof(float) );
            fftwf_execute(plan_trn);

            DPRINT("abs and square of fft\n");
            for (int t=0; t<nfftsteps; t++)
            {
                fft_out_squared[t] = cabs(fft_out_trn[t] * fft_out_trn[t]);
            }

            cblas_saxpy(nfftsteps, 1.0, fft_out_squared, 1,
                    &moltype_dos_raw_trn[h*nfftsteps], 1);

            DPRINT("rotational dof nr %d\n", i);
            memcpy(fft_in, &mol_omegas_sqrt_i_rot[i*ntrajsteps], ntrajsteps * sizeof(float) );
            fftwf_execute(plan_rot);

            DPRINT("abs and square of fft\n");
            for (int t=0; t<nfftsteps; t++)
            {
                fft_out_squared[t] = cabs(fft_out_rot[t] * fft_out_rot[t]);
            }

            cblas_saxpy(nfftsteps, 1.0, fft_out_squared, 1,
                    &moltype_dos_raw_rot[h*nfftsteps], 1);
            if ( i % 3 == 0 ) cblas_saxpy(nfftsteps, 1.0, fft_out_squared, 1, &moltype_dos_raw_rot_a[h*nfftsteps], 1);
            if ( i % 3 == 1 ) cblas_saxpy(nfftsteps, 1.0, fft_out_squared, 1, &moltype_dos_raw_rot_b[h*nfftsteps], 1);
            if ( i % 3 == 2 ) cblas_saxpy(nfftsteps, 1.0, fft_out_squared, 1, &moltype_dos_raw_rot_c[h*nfftsteps], 1);

            if (calc_cross)
            {
                DPRINT("cross terms\n");
                DPRINT("abs and multiply fft\n");
                // NOTE: this currently only correlates x with a, y with b and z with c
                // this is only meaningfull for isotropic systems, not for systems with a surface or electric field
                for (int t=0; t<nfftsteps; t++)
                {
                    fft_out_squared[t] = cabs(fft_out_rot[t] * fft_out_trn[t]);
                }

                if ( i % 3 == 0 ) cblas_saxpy(nfftsteps, 1.0, fft_out_squared, 1, &moltype_dos_raw_x_trn_rot_a[h*nfftsteps], 1);
                if ( i % 3 == 1 ) cblas_saxpy(nfftsteps, 1.0, fft_out_squared, 1, &moltype_dos_raw_x_trn_rot_b[h*nfftsteps], 1);
                if ( i % 3 == 2 ) cblas_saxpy(nfftsteps, 1.0, fft_out_squared, 1, &moltype_dos_raw_x_trn_rot_c[h*nfftsteps], 1);
                }
        }


        // vibration
        fftwf_complex* fft_out_vib = fftwf_malloc(sizeof(fftwf_complex) * nfftsteps);
        DPRINT("creating plan vibration\n");
        fftwf_plan plan_vib = fftwf_plan_dft_r2c_1d(ntrajsteps, fft_in, fft_out_vib, FFTW_MEASURE);

        int first_dof_vib = moltype_firstatom[h]*3;
        int last_dof_vib = moltype_firstatom[h]*3 + moltype_nmols[h]*moltype_natomtypes[h]*3;
        for (int i=first_dof_vib; i<last_dof_vib; i++)
        {
            DPRINT("vibrational dof nr %d\n", i);
            memcpy(fft_in, &atom_velocities_sqrt_m_vib[i*ntrajsteps], ntrajsteps * sizeof(float) );
            fftwf_execute(plan_vib);

            DPRINT("abs and square of fft\n");
            for (int t=0; t<nfftsteps; t++)
            {
                fft_out_squared[t] = cabs(fft_out_vib[t] * fft_out_vib[t]);
            }

            cblas_saxpy(nfftsteps, 1.0, fft_out_squared, 1,
                    &moltype_dos_raw_vib[h*nfftsteps], 1);
        }

        DPRINT("moltype done\n");

        fftwf_destroy_plan(plan_trn);
        fftwf_destroy_plan(plan_rot);
        fftwf_destroy_plan(plan_vib);
        free(fft_in);
        fftwf_free(fft_out_trn);
        fftwf_free(fft_out_rot);
        fftwf_free(fft_out_vib);
        free(fft_out_squared);
    }

    return 0;
}
