#include <stdlib.h>
#include <string.h>
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

int DOSCalculation (int nmoltypes,
        long ntrajsteps,
        long nfftsteps,
        int* moltype_firstmol,
        int* moltype_firstatom,
        int* moltype_nmols,
        int* moltype_natomspermol,
        float* mol_velocities_sqrt_m_trn,
        float* mol_omegas_sqrt_i_rot,
        float* atom_velocities_sqrt_m_vib,
        bool calc_components,
        bool calc_cross,
        float* moltype_dos_raw_trn,  // output
        float* moltype_dos_raw_trn_x,
        float* moltype_dos_raw_trn_y,
        float* moltype_dos_raw_trn_z,
        float* moltype_dos_raw_rot,
        float* moltype_dos_raw_rot_a,
        float* moltype_dos_raw_rot_b,
        float* moltype_dos_raw_rot_c,
        float* moltype_dos_raw_vib,
        float* moltype_dos_raw_vib_x,
        float* moltype_dos_raw_vib_y,
        float* moltype_dos_raw_vib_z,
        float* moltype_dos_raw_x_trn_rot,
        float* moltype_dos_raw_x_trn_vib,
        float* moltype_dos_raw_x_rot_vib)
{
    DPRINT("starting with FFT of translational and rotational dofs\n");

    for (int h=0; h<nmoltypes; h++)
    {
        DPRINT("moltype nr %d\n", h);
        float* fft_in = calloc(ntrajsteps, sizeof(float));
        fftwf_complex* fft_out_trnx = fftwf_malloc(sizeof(fftwf_complex) * nfftsteps);
        fftwf_complex* fft_out_trny = fftwf_malloc(sizeof(fftwf_complex) * nfftsteps);
        fftwf_complex* fft_out_trnz = fftwf_malloc(sizeof(fftwf_complex) * nfftsteps);
        fftwf_complex* fft_out_rota = fftwf_malloc(sizeof(fftwf_complex) * nfftsteps);
        fftwf_complex* fft_out_rotb = fftwf_malloc(sizeof(fftwf_complex) * nfftsteps);
        fftwf_complex* fft_out_rotc = fftwf_malloc(sizeof(fftwf_complex) * nfftsteps);
        fftwf_complex* fft_out_vib = fftwf_malloc(sizeof(fftwf_complex) * nfftsteps);
        float* fft_out_squared1 = calloc(nfftsteps, sizeof(float));
        float* fft_out_squared2 = calloc(nfftsteps, sizeof(float));
        float* fft_out_squared3 = calloc(nfftsteps, sizeof(float));
        DPRINT("creating plans translation\n");
        fftwf_plan plan_trnx = fftwf_plan_dft_r2c_1d(ntrajsteps, fft_in, fft_out_trnx, FFTW_MEASURE);
        fftwf_plan plan_trny = fftwf_plan_dft_r2c_1d(ntrajsteps, fft_in, fft_out_trny, FFTW_MEASURE);
        fftwf_plan plan_trnz = fftwf_plan_dft_r2c_1d(ntrajsteps, fft_in, fft_out_trnz, FFTW_MEASURE);
        DPRINT("creating plans rotational\n");
        fftwf_plan plan_rota = fftwf_plan_dft_r2c_1d(ntrajsteps, fft_in, fft_out_rota, FFTW_MEASURE);
        fftwf_plan plan_rotb = fftwf_plan_dft_r2c_1d(ntrajsteps, fft_in, fft_out_rotb, FFTW_MEASURE);
        fftwf_plan plan_rotc = fftwf_plan_dft_r2c_1d(ntrajsteps, fft_in, fft_out_rotc, FFTW_MEASURE);
        DPRINT("creating plan vibration\n");
        fftwf_plan plan_vib = fftwf_plan_dft_r2c_1d(ntrajsteps, fft_in, fft_out_vib, FFTW_MEASURE);

        int first_mol = moltype_firstmol[h];
        int last_mol = moltype_firstmol[h] + moltype_nmols[h];
        for (int i=first_mol; i<last_mol; i++)
        {
            DPRINT("molecule nr %d\n", i);
            // translation
            DPRINT("translational fft\n");
            memcpy(fft_in, &mol_velocities_sqrt_m_trn[(3*i+0)*ntrajsteps], ntrajsteps * sizeof(float) );
            fftwf_execute(plan_trnx);
            memcpy(fft_in, &mol_velocities_sqrt_m_trn[(3*i+1)*ntrajsteps], ntrajsteps * sizeof(float) );
            fftwf_execute(plan_trny);
            memcpy(fft_in, &mol_velocities_sqrt_m_trn[(3*i+2)*ntrajsteps], ntrajsteps * sizeof(float) );
            fftwf_execute(plan_trnz);

            DPRINT("abs and square of fft\n");
            for (int t=0; t<nfftsteps; t++)
            {
                fft_out_squared1[t] = cabs(fft_out_trnx[t] * fft_out_trnx[t]);
                fft_out_squared2[t] = cabs(fft_out_trny[t] * fft_out_trny[t]);
                fft_out_squared3[t] = cabs(fft_out_trnz[t] * fft_out_trnz[t]);
            }
            DPRINT("add to moltype dos\n");
            cblas_saxpy(nfftsteps, 1.0, fft_out_squared1, 1, &moltype_dos_raw_trn[h*nfftsteps], 1);
            cblas_saxpy(nfftsteps, 1.0, fft_out_squared2, 1, &moltype_dos_raw_trn[h*nfftsteps], 1);
            cblas_saxpy(nfftsteps, 1.0, fft_out_squared3, 1, &moltype_dos_raw_trn[h*nfftsteps], 1);

            if (calc_components)
            {
                DPRINT("add to moltype dos (trn components)\n");
                cblas_saxpy(nfftsteps, 1.0, fft_out_squared1, 1, &moltype_dos_raw_trn_x[h*nfftsteps], 1);
                cblas_saxpy(nfftsteps, 1.0, fft_out_squared2, 1, &moltype_dos_raw_trn_y[h*nfftsteps], 1);
                cblas_saxpy(nfftsteps, 1.0, fft_out_squared3, 1, &moltype_dos_raw_trn_z[h*nfftsteps], 1);
            }

            // rotation
            DPRINT("rotational fft\n");
            memcpy(fft_in, &mol_omegas_sqrt_i_rot[(3*i+0)*ntrajsteps], ntrajsteps * sizeof(float) );
            fftwf_execute(plan_rota);
            memcpy(fft_in, &mol_omegas_sqrt_i_rot[(3*i+1)*ntrajsteps], ntrajsteps * sizeof(float) );
            fftwf_execute(plan_rotb);
            memcpy(fft_in, &mol_omegas_sqrt_i_rot[(3*i+2)*ntrajsteps], ntrajsteps * sizeof(float) );
            fftwf_execute(plan_rotc);

            DPRINT("abs and square of fft\n");
            for (int t=0; t<nfftsteps; t++)
            {
                fft_out_squared1[t] = cabs(fft_out_rota[t] * fft_out_rota[t]);
                fft_out_squared2[t] = cabs(fft_out_rotb[t] * fft_out_rotb[t]);
                fft_out_squared3[t] = cabs(fft_out_rotc[t] * fft_out_rotc[t]);
            }

            DPRINT("add to moltype dos\n");
            cblas_saxpy(nfftsteps, 1.0, fft_out_squared1, 1, &moltype_dos_raw_rot[h*nfftsteps], 1);
            cblas_saxpy(nfftsteps, 1.0, fft_out_squared2, 1, &moltype_dos_raw_rot[h*nfftsteps], 1);
            cblas_saxpy(nfftsteps, 1.0, fft_out_squared3, 1, &moltype_dos_raw_rot[h*nfftsteps], 1);

            if (calc_components)
            {
                DPRINT("add to moltype dos (rot components)\n");
                cblas_saxpy(nfftsteps, 1.0, fft_out_squared1, 1, &moltype_dos_raw_rot_a[h*nfftsteps], 1);
                cblas_saxpy(nfftsteps, 1.0, fft_out_squared2, 1, &moltype_dos_raw_rot_b[h*nfftsteps], 1);
                cblas_saxpy(nfftsteps, 1.0, fft_out_squared3, 1, &moltype_dos_raw_rot_c[h*nfftsteps], 1);
            }

            // cross trn rot
            if (calc_cross)
            {
                DPRINT("cross terms\n");
                DPRINT("abs and multiply fft\n");
                for (int t=0; t<nfftsteps; t++)
                {
                    fft_out_squared1[t] = cabs(fft_out_trnx[t] * fft_out_rota[t]);
                    fft_out_squared2[t] = cabs(fft_out_trnx[t] * fft_out_rotb[t]);
                    fft_out_squared3[t] = cabs(fft_out_trnx[t] * fft_out_rotc[t]);
                }
                DPRINT("add to moltype dos\n");
                cblas_saxpy(nfftsteps, 1.0, fft_out_squared1, 1, &moltype_dos_raw_x_trn_rot[h*nfftsteps], 1);
                cblas_saxpy(nfftsteps, 1.0, fft_out_squared2, 1, &moltype_dos_raw_x_trn_rot[h*nfftsteps], 1);
                cblas_saxpy(nfftsteps, 1.0, fft_out_squared3, 1, &moltype_dos_raw_x_trn_rot[h*nfftsteps], 1);

                for (int t=0; t<nfftsteps; t++)
                {
                    fft_out_squared1[t] = cabs(fft_out_trny[t] * fft_out_rota[t]);
                    fft_out_squared2[t] = cabs(fft_out_trny[t] * fft_out_rotb[t]);
                    fft_out_squared3[t] = cabs(fft_out_trny[t] * fft_out_rotc[t]);
                }
                cblas_saxpy(nfftsteps, 1.0, fft_out_squared1, 1, &moltype_dos_raw_x_trn_rot[h*nfftsteps], 1);
                cblas_saxpy(nfftsteps, 1.0, fft_out_squared2, 1, &moltype_dos_raw_x_trn_rot[h*nfftsteps], 1);
                cblas_saxpy(nfftsteps, 1.0, fft_out_squared3, 1, &moltype_dos_raw_x_trn_rot[h*nfftsteps], 1);

                for (int t=0; t<nfftsteps; t++)
                {
                    fft_out_squared1[t] = cabs(fft_out_trnz[t] * fft_out_rota[t]);
                    fft_out_squared2[t] = cabs(fft_out_trnz[t] * fft_out_rotb[t]);
                    fft_out_squared3[t] = cabs(fft_out_trnz[t] * fft_out_rotc[t]);
                }
                cblas_saxpy(nfftsteps, 1.0, fft_out_squared1, 1, &moltype_dos_raw_x_trn_rot[h*nfftsteps], 1);
                cblas_saxpy(nfftsteps, 1.0, fft_out_squared2, 1, &moltype_dos_raw_x_trn_rot[h*nfftsteps], 1);
                cblas_saxpy(nfftsteps, 1.0, fft_out_squared3, 1, &moltype_dos_raw_x_trn_rot[h*nfftsteps], 1);
            }

            // vibration and other cross terms
            DPRINT("vibrational fft\n");
            int first_dof_vib = 3*moltype_firstatom[h] + 3*(i-moltype_firstmol[h])*moltype_natomspermol[h];
            int last_dof_vib = first_dof_vib + 3*moltype_natomspermol[h];
            // loop over all vibrational dof of the molecule
            for (int df=first_dof_vib; df<last_dof_vib; df++)
            {
                DPRINT("vibrational fft nr %d\n", df);
                memcpy(fft_in, &atom_velocities_sqrt_m_vib[df*ntrajsteps], ntrajsteps * sizeof(float) );
                fftwf_execute(plan_vib);

                DPRINT("abs and square of fft\n");
                for (int t=0; t<nfftsteps; t++)
                {
                    fft_out_squared1[t] = cabs(fft_out_vib[t] * fft_out_vib[t]);
                }

                DPRINT("add to moltype dos\n");
                cblas_saxpy(nfftsteps, 1.0, fft_out_squared1, 1, &moltype_dos_raw_vib[h*nfftsteps], 1);

                if (calc_components)
                {
                    DPRINT("add to moltype dos (vib components)\n");
                    if (df % 3 == 0) {cblas_saxpy(nfftsteps, 1.0, fft_out_squared1, 1, &moltype_dos_raw_vib_x[h*nfftsteps], 1);}
                    if (df % 3 == 1) {cblas_saxpy(nfftsteps, 1.0, fft_out_squared1, 1, &moltype_dos_raw_vib_y[h*nfftsteps], 1);}
                    if (df % 3 == 2) {cblas_saxpy(nfftsteps, 1.0, fft_out_squared1, 1, &moltype_dos_raw_vib_z[h*nfftsteps], 1);}
                }

                DPRINT("cross terms trn vib\n");
                DPRINT("abs and multiply fft\n");
                for (int t=0; t<nfftsteps; t++)
                {
                    fft_out_squared1[t] = cabs(fft_out_trnx[t] * fft_out_vib[t]);
                    fft_out_squared2[t] = cabs(fft_out_trny[t] * fft_out_vib[t]);
                    fft_out_squared3[t] = cabs(fft_out_trnz[t] * fft_out_vib[t]);
                }
                DPRINT("add to moltype dos\n");
                cblas_saxpy(nfftsteps, 1.0, fft_out_squared1, 1, &moltype_dos_raw_x_trn_vib[h*nfftsteps], 1);
                cblas_saxpy(nfftsteps, 1.0, fft_out_squared2, 1, &moltype_dos_raw_x_trn_vib[h*nfftsteps], 1);
                cblas_saxpy(nfftsteps, 1.0, fft_out_squared3, 1, &moltype_dos_raw_x_trn_vib[h*nfftsteps], 1);

                DPRINT("cross terms rot vib\n");
                DPRINT("abs and multiply fft\n");
                for (int t=0; t<nfftsteps; t++)
                {
                    fft_out_squared1[t] = cabs(fft_out_rota[t] * fft_out_vib[t]);
                    fft_out_squared2[t] = cabs(fft_out_rotb[t] * fft_out_vib[t]);
                    fft_out_squared3[t] = cabs(fft_out_rotc[t] * fft_out_vib[t]);
                }
                DPRINT("add to moltype dos\n");
                cblas_saxpy(nfftsteps, 1.0, fft_out_squared1, 1, &moltype_dos_raw_x_rot_vib[h*nfftsteps], 1);
                cblas_saxpy(nfftsteps, 1.0, fft_out_squared2, 1, &moltype_dos_raw_x_rot_vib[h*nfftsteps], 1);
                cblas_saxpy(nfftsteps, 1.0, fft_out_squared3, 1, &moltype_dos_raw_x_rot_vib[h*nfftsteps], 1);
            }
        }
        DPRINT("moltype done\n");

        fftwf_destroy_plan(plan_trnx);
        fftwf_destroy_plan(plan_trny);
        fftwf_destroy_plan(plan_trnz);
        fftwf_destroy_plan(plan_rota);
        fftwf_destroy_plan(plan_rotb);
        fftwf_destroy_plan(plan_rotc);
        fftwf_destroy_plan(plan_vib);
        free(fft_in);
        fftwf_free(fft_out_trnx);
        fftwf_free(fft_out_trny);
        fftwf_free(fft_out_trnz);
        fftwf_free(fft_out_rota);
        fftwf_free(fft_out_rotb);
        fftwf_free(fft_out_rotc);
        fftwf_free(fft_out_vib);
        free(fft_out_squared1);
        free(fft_out_squared2);
        free(fft_out_squared3);
    }

    return 0;
}
