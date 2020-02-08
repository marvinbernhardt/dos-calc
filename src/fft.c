#include "structs.h"
#include <cblas.h>
#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#ifdef DEBUG
#define DPRINT(...)                                                            \
  do {                                                                         \
    fprintf(stderr, __VA_ARGS__);                                              \
  } while (0)
#else
#define DPRINT(...)
#endif

size_t gen_dof_fourier_index(unsigned long nfrequencies,
                             size_t *moltype_firstmol,
                             size_t *moltype_firstatom,
                             size_t *moltype_natomspermol,
                             size_t moltype, // query variables
                             char type, size_t dof, size_t i0) {
  size_t mol_natoms = moltype_natomspermol[moltype];
  size_t trn_start = 0;
  size_t rot_xyz_start = 3;
  size_t vib_start = rot_xyz_start + 3 * mol_natoms;
  size_t rot_omega_start = vib_start + 3 * mol_natoms;

  // add specified type
  size_t dof_index = 0;
  if (type == 't')
    dof_index += trn_start;
  else if (type == 'r')
    dof_index += rot_xyz_start;
  else if (type == 'v')
    dof_index += vib_start;
  else if (type == 'o')
    dof_index += rot_omega_start;
  else {
    fprintf(stderr, "ERROR: unknown cross spectrum dof_type: '%c'\n", type);
    exit(1);
  }

  // add specified dof
  dof_index += dof;

  size_t dof_fourier_index =
      nfrequencies *
      (6 * moltype_firstmol[moltype] + 6 * moltype_firstatom[moltype] + 6 * i0 +
       6 * mol_natoms * i0 + dof_index);
  return dof_fourier_index;
}

void dos_calculation(
    size_t nmoltypes, unsigned long nblocksteps, unsigned long nfrequencies,
    size_t *moltype_firstmol, size_t *moltype_firstatom, size_t *moltype_nmols,
    size_t *moltype_natomspermol, float *mol_velocities_sqrt_m_trn,
    float *mol_omegas_sqrt_i_rot, float *atom_velocities_sqrt_m_vib,
    float *atom_velocities_sqrt_m_rot, size_t ndos, size_t nsamples,
    size_t sample, size_t ncross_spectra, cross_spectrum_def *cross_spectra_def,
    float *moltypes_dos_samples, // output
    float *cross_spectra_samples) {
  // finding the number of dof
  size_t ndof = 0;
  for (size_t h = 0; h < nmoltypes; h++) {
    // trn and rot_omega
    ndof += 6 * moltype_nmols[h];
    // vib and rot_xyz
    ndof += 6 * moltype_nmols[h] * moltype_natomspermol[h];
  }

  // array that will hold all FT of the time series
  // this is for cross spectra calculation later
  // order of dof is: trn rot_xyz vib rot_omega
  fftwf_complex *dof_fourier =
      calloc(ndof * nfrequencies, sizeof(fftwf_complex));

  // stuff for fftw
  float *fft_in = calloc(nblocksteps, sizeof(float));
  fftwf_complex *fft_out = fftwf_malloc(sizeof(fftwf_complex) * nfrequencies);
  float *fft_out_squared = calloc(nfrequencies, sizeof(float));
  DPRINT("creating plan\n");
  fftwf_plan plan =
      fftwf_plan_dft_r2c_1d(nblocksteps, fft_in, fft_out, FFTW_MEASURE);

  // fourier all dof
  for (size_t h = 0; h < nmoltypes; h++) {
    DPRINT("moltype nr %zu\n", h);

    size_t first_mol = moltype_firstmol[h];
    size_t last_mol = moltype_firstmol[h] + moltype_nmols[h];
    for (size_t i = first_mol; i < last_mol; i++) {
      DPRINT("molecule nr %zu\n", i);
      size_t i0 = i - moltype_firstmol[h];

      // convenience stuff
      size_t mol_natoms = moltype_natomspermol[h];
      size_t mol_ndof = 6 + 6 * mol_natoms;
      size_t trn_start = 0;
      size_t trn_end = 3;
      size_t rot_xyz_start = trn_end;
      size_t rot_xyz_end = trn_end + 3 * mol_natoms;
      size_t vib_start = rot_xyz_end;
      size_t vib_end = rot_xyz_end + 3 * mol_natoms;
      size_t rot_omega_start = vib_end;
      size_t rot_omega_end = vib_end + 3;
      size_t mol_first_dof_atomic =
          3 * moltype_firstatom[h] + 3 * i0 * mol_natoms;
      size_t dos;

      for (size_t dof = 0; dof < mol_ndof; dof++) {
        size_t xyz = dof % 3;
        if ((dof >= trn_start) && (dof < trn_end)) {
          DPRINT("translational fft\n");
          memcpy(fft_in,
                 &mol_velocities_sqrt_m_trn[(3 * i + xyz) * nblocksteps],
                 nblocksteps * sizeof(float));
          dos = 0 + xyz;
        }
        if ((dof >= rot_xyz_start) && (dof < rot_xyz_end)) {
          DPRINT("rotational_xyz fft\n");
          size_t dof_rot = dof - rot_xyz_start;
          memcpy(fft_in,
                 &atom_velocities_sqrt_m_rot[(mol_first_dof_atomic + dof_rot) *
                                             nblocksteps],
                 nblocksteps * sizeof(float));
          dos = 3 + xyz;
        }
        if ((dof >= vib_start) && (dof < vib_end)) {
          DPRINT("vibrational fft\n");
          size_t d_vib = dof - vib_start;
          memcpy(fft_in,
                 &atom_velocities_sqrt_m_vib[(mol_first_dof_atomic + d_vib) *
                                             nblocksteps],
                 nblocksteps * sizeof(float));
          dos = 6 + xyz;
        }
        if ((dof >= rot_omega_start) && (dof < rot_omega_end)) {
          DPRINT("rotational_omega fft\n");
          memcpy(fft_in, &mol_omegas_sqrt_i_rot[(3 * i + xyz) * nblocksteps],
                 nblocksteps * sizeof(float));
          dos = 9 + xyz;
        }

        // execute fftw
        fftwf_execute(plan);

        // save all fourier transforms
        size_t dof_index =
            nfrequencies * (6 * moltype_firstmol[h] + 6 * moltype_firstatom[h] +
                            6 * i0 + 6 * mol_natoms * i0 + dof);
        memcpy(&dof_fourier[dof_index], fft_out,
               nfrequencies * sizeof(fftwf_complex));

        // square and add to dos
        for (unsigned long t = 0; t < nfrequencies; t++) {
          fft_out_squared[t] = cabs(fft_out[t] * fft_out[t]);
        }
        size_t dos_index = h * ndos * nsamples * nfrequencies +
                           dos * nsamples * nfrequencies +
                           sample * nfrequencies;
        cblas_saxpy(nfrequencies, 1.0, fft_out_squared, 1,
                    &moltypes_dos_samples[dos_index], 1);
      }
    }
    DPRINT("moltype done\n");
  }
  DPRINT("all moltypes done\n");

  // fftwf_complex *dof_fourier = calloc(ndof*nfrequencies,
  // sizeof(fftwf_complex));
  for (size_t d = 0; d < ndof; d++) {
    DPRINT("dof_fourier of dof %zu:\n", d);
    for (size_t t = 0; t < nfrequencies; t++) {
      DPRINT("%f + i%f, ", creal(dof_fourier[d * nfrequencies + t]),
             cimag(dof_fourier[d * nfrequencies + t]));
    }
    DPRINT("\n");
  }

  DPRINT("begin with cross spectra\n");
  for (size_t d = 0; d < ncross_spectra; d++) {
    DPRINT("Cross spectrum name: %s\n", cross_spectra_def[d].name);
    char cross_spectrum_type = cross_spectra_def[d].type;
    size_t ncross_contribs = 0;
    float *cross_spectrum = calloc(nfrequencies, sizeof(float));
    for (size_t p = 0; p < cross_spectra_def[d].ndof_pair_defs; p++) {
      // convenience
      size_t moltypeA = cross_spectra_def[d].dof_pair_defs[p].dofA_moltype;
      size_t moltypeB = cross_spectra_def[d].dof_pair_defs[p].dofB_moltype;
      char typeA = cross_spectra_def[d].dof_pair_defs[p].dofA_type;
      char typeB = cross_spectra_def[d].dof_pair_defs[p].dofB_type;

      for (size_t dA = 0; dA < cross_spectra_def[d].dof_pair_defs[p].ndofA;
           dA++) {
        for (size_t dB = 0; dB < cross_spectra_def[d].dof_pair_defs[p].ndofB;
             dB++) {
          // convenience
          size_t dofA = cross_spectra_def[d].dof_pair_defs[p].dofA_list[dA];
          size_t dofB = cross_spectra_def[d].dof_pair_defs[p].dofB_list[dB];

          DPRINT("Cross spectrum contribution: %zu %c %zu - %zu %c %zu\n",
                 moltypeA, typeA, dofA, moltypeB, typeB, dofB);

          if (cross_spectrum_type == 'e') {
            size_t nmolsA = moltype_nmols[moltypeA];
            size_t nmolsB = moltype_nmols[moltypeB];
            for (size_t iA = 0; iA < nmolsA; iA++) {
              for (size_t iB = 0; iB < nmolsB; iB++) {
                ncross_contribs++;
                // find index in
                size_t dof_fourier_indexA = gen_dof_fourier_index(
                    nfrequencies, moltype_firstmol, moltype_firstatom,
                    moltype_natomspermol, moltypeA, typeA, dofA, iA);
                size_t dof_fourier_indexB = gen_dof_fourier_index(
                    nfrequencies, moltype_firstmol, moltype_firstatom,
                    moltype_natomspermol, moltypeB, typeB, dofB, iB);

                for (unsigned long t = 0; t < nfrequencies; t++) {
                  fft_out_squared[t] =
                      cabs(dof_fourier[dof_fourier_indexA + t] *
                           dof_fourier[dof_fourier_indexB + t]);
                }

                cblas_saxpy(nfrequencies, 1.0, fft_out_squared, 1,
                            cross_spectrum, 1);
              }
            }
          } else if (cross_spectrum_type == 'i') {
            size_t nmolsA = moltype_nmols[moltypeA];
            for (size_t iA = 0; iA < nmolsA; iA++) {
              ncross_contribs++;
              // find index in
              size_t dof_fourier_indexA = gen_dof_fourier_index(
                  nfrequencies, moltype_firstmol, moltype_firstatom,
                  moltype_natomspermol, moltypeA, typeA, dofA, iA);
              size_t dof_fourier_indexB = gen_dof_fourier_index(
                  nfrequencies, moltype_firstmol, moltype_firstatom,
                  moltype_natomspermol, moltypeB, typeB, dofB, iA);

              for (unsigned long t = 0; t < nfrequencies; t++) {
                fft_out_squared[t] = cabs(dof_fourier[dof_fourier_indexA + t] *
                                          dof_fourier[dof_fourier_indexB + t]);
              }

              cblas_saxpy(nfrequencies, 1.0, fft_out_squared, 1, cross_spectrum,
                          1);
            }
          }
        }
      }
    }

    // scale by 1/ncross_contribs and add to array
    size_t cross_index = d * nsamples * nfrequencies + sample * nfrequencies;
    DPRINT("ncross_contribs: %zu\n", ncross_contribs);
    cblas_saxpy(nfrequencies, 1.0 / (float)ncross_contribs, cross_spectrum, 1,
                &cross_spectra_samples[cross_index], 1);
    free(cross_spectrum);
  }

  fftwf_destroy_plan(plan);
  free(fft_in);
  fftwf_free(fft_out);
  free(fft_out_squared);
  free(dof_fourier);
}
