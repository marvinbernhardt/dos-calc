#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <stdbool.h>
#include <math.h>
#include <cblas.h>
#include <fftw3.h>

#ifdef DEBUG
#define DPRINT(...) do{ fprintf( stdout, __VA_ARGS__ ); } while( 0 )
#else
#define DPRINT(...)
#endif

void DOSCalculation (int nmoltypes,
        long nblocksteps,
        long nfftsteps,
        int* moltype_firstmol,
        int* moltype_firstatom,
        int* moltype_nmols,
        int* moltype_natomspermol,
        float* mol_velocities_sqrt_m_trn,
        float* mol_omegas_sqrt_i_rot,
        float* atom_velocities_sqrt_m_vib,
        float* atom_velocities_sqrt_m_rot,
        bool calc_components,
        bool calc_cross,
        bool calc_rot_alt,
        float* moltype_dos_raw_trn,  // output
        float* moltype_dos_raw_trn_x,
        float* moltype_dos_raw_trn_y,
        float* moltype_dos_raw_trn_z,
        float* moltype_dos_raw_rot,
        float* moltype_dos_raw_rot_a,
        float* moltype_dos_raw_rot_b,
        float* moltype_dos_raw_rot_c,
        float* moltype_dos_raw_rot_alt,
        float* moltype_dos_raw_rot_alt_x,
        float* moltype_dos_raw_rot_alt_y,
        float* moltype_dos_raw_rot_alt_z,
        float* moltype_dos_raw_vib,
        float* moltype_dos_raw_vib_x,
        float* moltype_dos_raw_vib_y,
        float* moltype_dos_raw_vib_z,
        float* moltype_dos_raw_x_trn_rot,
        float* moltype_dos_raw_x_trn_vib,
        float* moltype_dos_raw_x_rot_vib)
{
    DPRINT("starting with FFT of translational and rotational dofs\n");

    // no dynamic teams!!!
    omp_set_dynamic(0);

    for (int h=0; h<nmoltypes; h++)
    {
        DPRINT("moltype nr %d\n", h);
        static fftwf_complex* fft_out_trnx;
        static fftwf_complex* fft_out_trny;
        static fftwf_complex* fft_out_trnz;
        static fftwf_complex* fft_out_rota;
        static fftwf_complex* fft_out_rotb;
        static fftwf_complex* fft_out_rotc;
        static fftwf_complex* fft_out_vib;
        static fftwf_complex* fft_out_rotalt;
        #pragma omp threadprivate(fft_out_trnx, fft_out_trny, fft_out_trnz, fft_out_rota, fft_out_rotb, fft_out_rotc, fft_out_vib, fft_out_rotalt)

        #pragma omp parallel
        {
            fft_out_trnx = fftwf_malloc(sizeof(fftwf_complex) * nfftsteps);
            fft_out_trny = fftwf_malloc(sizeof(fftwf_complex) * nfftsteps);
            fft_out_trnz = fftwf_malloc(sizeof(fftwf_complex) * nfftsteps);
            fft_out_rota = fftwf_malloc(sizeof(fftwf_complex) * nfftsteps);
            fft_out_rotb = fftwf_malloc(sizeof(fftwf_complex) * nfftsteps);
            fft_out_rotc = fftwf_malloc(sizeof(fftwf_complex) * nfftsteps);
            fft_out_vib = fftwf_malloc(sizeof(fftwf_complex) * nfftsteps);
            fft_out_rotalt = fftwf_malloc(sizeof(fftwf_complex) * nfftsteps);
        }
        static float* fft_in;
        static float* fft_out_squared1;
        static float* fft_out_squared2;
        static float* fft_out_squared3;
        #pragma omp threadprivate(fft_in, fft_out_squared1, fft_out_squared2, fft_out_squared3)
        #pragma omp parallel
        {
            fft_in = calloc(nblocksteps, sizeof(float));
            fft_out_squared1 = calloc(nfftsteps, sizeof(float));
            fft_out_squared2 = calloc(nfftsteps, sizeof(float));
            fft_out_squared3 = calloc(nfftsteps, sizeof(float));
        }
        DPRINT("creating fftw plans\n");
        static fftwf_plan plan_trnx;
        static fftwf_plan plan_trny;
        static fftwf_plan plan_trnz;
        static fftwf_plan plan_rota;
        static fftwf_plan plan_rotb;
        static fftwf_plan plan_rotc;
        static fftwf_plan plan_vib;
        static fftwf_plan plan_rotalt;
        #pragma omp threadprivate(plan_trnx, plan_trny, plan_trnz, plan_rota, plan_rotb, plan_rotc, plan_vib, plan_rotalt)
        #pragma omp parallel
        {
            plan_trnx = fftwf_plan_dft_r2c_1d(nblocksteps, fft_in, fft_out_trnx, FFTW_MEASURE);
            plan_trny = fftwf_plan_dft_r2c_1d(nblocksteps, fft_in, fft_out_trny, FFTW_MEASURE);
            plan_trnz = fftwf_plan_dft_r2c_1d(nblocksteps, fft_in, fft_out_trnz, FFTW_MEASURE);
            plan_rota = fftwf_plan_dft_r2c_1d(nblocksteps, fft_in, fft_out_rota, FFTW_MEASURE);
            plan_rotb = fftwf_plan_dft_r2c_1d(nblocksteps, fft_in, fft_out_rotb, FFTW_MEASURE);
            plan_rotc = fftwf_plan_dft_r2c_1d(nblocksteps, fft_in, fft_out_rotc, FFTW_MEASURE);
            plan_vib = fftwf_plan_dft_r2c_1d(nblocksteps, fft_in, fft_out_vib, FFTW_MEASURE);
            plan_rotalt = fftwf_plan_dft_r2c_1d(nblocksteps, fft_in, fft_out_rotalt, FFTW_MEASURE);
        }

        int first_mol = moltype_firstmol[h];
        int last_mol = moltype_firstmol[h] + moltype_nmols[h];
        #pragma omp parallel for
        for (int i=first_mol; i<last_mol; i++)
        {
            DPRINT("molecule nr %d\n", i);
            // translation
            DPRINT("translational fft\n");
            memcpy(fft_in, &mol_velocities_sqrt_m_trn[(3*i+0)*nblocksteps], nblocksteps * sizeof(float) );
            fftwf_execute(plan_trnx);
            memcpy(fft_in, &mol_velocities_sqrt_m_trn[(3*i+1)*nblocksteps], nblocksteps * sizeof(float) );
            fftwf_execute(plan_trny);
            memcpy(fft_in, &mol_velocities_sqrt_m_trn[(3*i+2)*nblocksteps], nblocksteps * sizeof(float) );
            fftwf_execute(plan_trnz);

            DPRINT("abs and square of fft\n");
            for (int t=0; t<nfftsteps; t++)
            {
                fft_out_squared1[t] = cabs(fft_out_trnx[t] * fft_out_trnx[t]);
                fft_out_squared2[t] = cabs(fft_out_trny[t] * fft_out_trny[t]);
                fft_out_squared3[t] = cabs(fft_out_trnz[t] * fft_out_trnz[t]);
            }
            DPRINT("add to moltype dos\n");
            #pragma omp critical (add_dos_trn)
            {
                cblas_saxpy(nfftsteps, 1.0, fft_out_squared1, 1, &moltype_dos_raw_trn[h*nfftsteps], 1);
                cblas_saxpy(nfftsteps, 1.0, fft_out_squared2, 1, &moltype_dos_raw_trn[h*nfftsteps], 1);
                cblas_saxpy(nfftsteps, 1.0, fft_out_squared3, 1, &moltype_dos_raw_trn[h*nfftsteps], 1);

                if (calc_components)
                {
                    DPRINT("add to moltype dos (trn components)\n");
                    {
                        cblas_saxpy(nfftsteps, 1.0, fft_out_squared1, 1, &moltype_dos_raw_trn_x[h*nfftsteps], 1);
                        cblas_saxpy(nfftsteps, 1.0, fft_out_squared2, 1, &moltype_dos_raw_trn_y[h*nfftsteps], 1);
                        cblas_saxpy(nfftsteps, 1.0, fft_out_squared3, 1, &moltype_dos_raw_trn_z[h*nfftsteps], 1);
                    }
                }
            }

            // rotation
            DPRINT("rotational fft\n");
            memcpy(fft_in, &mol_omegas_sqrt_i_rot[(3*i+0)*nblocksteps], nblocksteps * sizeof(float) );
            fftwf_execute(plan_rota);
            memcpy(fft_in, &mol_omegas_sqrt_i_rot[(3*i+1)*nblocksteps], nblocksteps * sizeof(float) );
            fftwf_execute(plan_rotb);
            memcpy(fft_in, &mol_omegas_sqrt_i_rot[(3*i+2)*nblocksteps], nblocksteps * sizeof(float) );
            fftwf_execute(plan_rotc);

            DPRINT("abs and square of fft\n");
            for (int t=0; t<nfftsteps; t++)
            {
                fft_out_squared1[t] = cabs(fft_out_rota[t] * fft_out_rota[t]);
                fft_out_squared2[t] = cabs(fft_out_rotb[t] * fft_out_rotb[t]);
                fft_out_squared3[t] = cabs(fft_out_rotc[t] * fft_out_rotc[t]);
            }

            DPRINT("add to moltype dos\n");
            #pragma omp critical (add_dos_rot)
            {
                cblas_saxpy(nfftsteps, 1.0, fft_out_squared1, 1, &moltype_dos_raw_rot[h*nfftsteps], 1);
                cblas_saxpy(nfftsteps, 1.0, fft_out_squared2, 1, &moltype_dos_raw_rot[h*nfftsteps], 1);
                cblas_saxpy(nfftsteps, 1.0, fft_out_squared3, 1, &moltype_dos_raw_rot[h*nfftsteps], 1);

                if (calc_components)
                {
                    DPRINT("add to moltype dos (rot components)\n");
                    {
                        cblas_saxpy(nfftsteps, 1.0, fft_out_squared1, 1, &moltype_dos_raw_rot_a[h*nfftsteps], 1);
                        cblas_saxpy(nfftsteps, 1.0, fft_out_squared2, 1, &moltype_dos_raw_rot_b[h*nfftsteps], 1);
                        cblas_saxpy(nfftsteps, 1.0, fft_out_squared3, 1, &moltype_dos_raw_rot_c[h*nfftsteps], 1);
                    }
                }
            }

            // cross trn rot
            if (calc_cross)
            {
                #pragma omp critical (add_dos_cross)
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
            }
        }

        for (int i=first_mol; i<last_mol; i++)
        {
            // vibration and other cross terms
            // and alternative rotation
            DPRINT("vibrational fft\n");
            int first_dof_vib = 3*moltype_firstatom[h] + 3*(i-moltype_firstmol[h])*moltype_natomspermol[h];
            int last_dof_vib = first_dof_vib + 3*moltype_natomspermol[h];
            // loop over all vibrational dof of the molecule
            for (int df=first_dof_vib; df<last_dof_vib; df++)
            {
                DPRINT("vibrational fft nr %d\n", df);
                memcpy(fft_in, &atom_velocities_sqrt_m_vib[df*nblocksteps], nblocksteps * sizeof(float) );
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

                if (calc_rot_alt)
                {
                    DPRINT("rotational alternative fft nr %d\n", df);
                    memcpy(fft_in, &atom_velocities_sqrt_m_rot[df*nblocksteps], nblocksteps * sizeof(float) );
                    fftwf_execute(plan_rotalt);

                    DPRINT("abs and square of fft\n");
                    for (int t=0; t<nfftsteps; t++)
                    {
                        fft_out_squared1[t] = cabs(fft_out_rotalt[t] * fft_out_rotalt[t]);
                    }

                    DPRINT("add to moltype dos\n");
                    cblas_saxpy(nfftsteps, 1.0, fft_out_squared1, 1, &moltype_dos_raw_rot_alt[h*nfftsteps], 1);

                    if (calc_components)
                    {
                        DPRINT("add to moltype dos (rotalt components)\n");
                        if (df % 3 == 0) {cblas_saxpy(nfftsteps, 1.0, fft_out_squared1, 1, &moltype_dos_raw_rot_alt_x[h*nfftsteps], 1);}
                        if (df % 3 == 1) {cblas_saxpy(nfftsteps, 1.0, fft_out_squared1, 1, &moltype_dos_raw_rot_alt_y[h*nfftsteps], 1);}
                        if (df % 3 == 2) {cblas_saxpy(nfftsteps, 1.0, fft_out_squared1, 1, &moltype_dos_raw_rot_alt_z[h*nfftsteps], 1);}
                    }
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
        fftwf_destroy_plan(plan_rotalt);
        free(fft_in);
        fftwf_free(fft_out_trnx);
        fftwf_free(fft_out_trny);
        fftwf_free(fft_out_trnz);
        fftwf_free(fft_out_rota);
        fftwf_free(fft_out_rotb);
        fftwf_free(fft_out_rotc);
        fftwf_free(fft_out_vib);
        fftwf_free(fft_out_rotalt);
        free(fft_out_squared1);
        free(fft_out_squared2);
        free(fft_out_squared3);
    }
}
