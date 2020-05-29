#include "linear-algebra.c"
#include <cblas.h>
#include <lapacke.h>
#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

float sqrt_neg_zero(float number) {
    if (number < -0.001) {
        printf("Trying to take square-root of negative number!\n");
        printf("%f\n", number);
        exit(0);
    } else if (number < 0)
        return 0;
    else
        return sqrt(number);
}

void recombine_molecule(float *box, size_t m_natoms, float *positions) {
    for (size_t dim = 0; dim < 3; dim++) {
        for (size_t j = 1; j < m_natoms; j++) {
            float dist_to_firstatom =
                positions[3 * j + dim] - positions[3 * 0 + dim];
            if (dist_to_firstatom > 0.5 * box[dim]) {
                positions[3 * j + dim] -= box[dim];
            }
            if (dist_to_firstatom < -0.5 * box[dim]) {
                positions[3 * j + dim] += box[dim];
            }
        }
    }
}

void calculate_refpos_principal_components(
    float *refpos, float *box, size_t nmols, size_t *mol_firstatom,
    size_t *mol_natoms, size_t *mol_moltypenr, float **moltypes_atommasses,
    float *mol_mass, char *moltype_rot_treat, bool no_pbc,
    float *atom_refpos_principal_components) // from here output
{
    // iterate molecules
    for (size_t i = 0; i < nmols; i++) {
        size_t m_firstatom = mol_firstatom[i];
        size_t m_natoms = mol_natoms[i];
        size_t m_moltype = mol_moltypenr[i];
        float m_mass = mol_mass[i];
        float *m_atommasses = moltypes_atommasses[m_moltype];
        char m_rot_treat = moltype_rot_treat[m_moltype];
        // large arrays
        float *positions = calloc(3 * m_natoms, sizeof(float));
        float *positions_rel = calloc(3 * m_natoms, sizeof(float));
        // next mol if not eckart frame decomposition
        if (!(m_rot_treat == 'e' || m_rot_treat == 'p')) {
            continue;
        }
        // reading molecule positions
        for (size_t j = 0; j < m_natoms; j++) {
            size_t jj = m_firstatom + j;
            positions[3 * j + 0] = refpos[3 * jj + 0];
            positions[3 * j + 1] = refpos[3 * jj + 1];
            positions[3 * j + 2] = refpos[3 * jj + 2];
        }
        // recombination
        if (no_pbc == false) {
            recombine_molecule(box, m_natoms, positions);
        }
        // calc molecule com
        float center_of_mass[3] = {0.0, 0.0, 0.0};
        for (size_t dim = 0; dim < 3; dim++) {
            center_of_mass[dim] =
                cblas_sdot(m_natoms, m_atommasses, 1, &positions[0 + dim], 3);
            center_of_mass[dim] /= m_mass;
        }
        // calc molecule atoms relative positions
        for (size_t j = 0; j < m_natoms; j++) {
            for (size_t dim = 0; dim < 3; dim++) {
                positions_rel[3 * j + dim] =
                    positions[3 * j + dim] - center_of_mass[dim];
            }
        }
        // calc moi tensor
        float moi_tensor[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        moiTensor(m_natoms, positions_rel, m_atommasses, moi_tensor);
        // calc moments of inertia and eigenvectors
        float eigenvectors[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        cblas_scopy(9, moi_tensor, 1, eigenvectors, 1);
        float moments_of_inertia[3] = {0.0, 0.0, 0.0};
        if (LAPACKE_ssyev(LAPACK_ROW_MAJOR, 'V', 'U', 3, eigenvectors, 3,
                          moments_of_inertia) > 0) {
            fprintf(stderr,
                    "ERROR: LAPACKE_ssyev failed to compute eigenvalues of\n");
            fprintf(
                stderr,
                "       the refpos moment of inertia tensor of molecule %zu\n",
                i);
            exit(1);
        }
        // first dimension: atom, second dimension: dim
        float *refpos_principal_components =
            calloc(3 * m_natoms, sizeof(float));
        float eigenvectors_inv[9] = {0.0, 0.0, 0.0, 0.0, 0.0,
                                     0.0, 0.0, 0.0, 0.0};
        cblas_scopy(9, eigenvectors, 1, eigenvectors_inv, 1);
        if (invert_matrix(eigenvectors_inv, 3) != 0) {
            fprintf(stderr,
                    "ERROR: non-invertible moi-eigenvector matrix of molecule "
                    "%zu\n",
                    i);
            exit(1);
        }
        // c = inv(v) @ posref_rel
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 3, m_natoms, 3,
                    1.0, eigenvectors_inv, 3, positions_rel, 3, 0.0,
                    refpos_principal_components, m_natoms);
        // save in array
        cblas_scopy(m_natoms * 3, refpos_principal_components, 1,
                    &atom_refpos_principal_components[m_firstatom * 3], 1);
        // free arrays
        free(refpos_principal_components);
        free(positions);
        free(positions_rel);
        // TODO: check first eckart condition
    }
}

void decompose_velocities(
    float *block_pos, float *block_vel, float *block_box,
    unsigned long nblocksteps, size_t natoms, size_t nmols,
    size_t *mol_firstatom, size_t *mol_natoms, size_t *mol_moltypenr,
    float **moltypes_atommasses, float *mol_mass, char *moltype_rot_treat,
    int **moltype_abc_indicators, bool no_pbc,
    float *atom_refpos_principal_components,
    float *mol_velocities_sqrt_m_trn, // from here output
    float *mol_omegas_sqrt_i_rot, float *atom_velocities_sqrt_m_vib,
    float *atom_velocities_sqrt_m_rot, float *atom_velocities_sqrt_m_vibc,
    float *mol_block_moments_of_inertia,
    float *mol_block_moments_of_inertia_squared, float *mol_block_coriolis) {

    // no dynamic teams since molecules all cause roughly the same work
    omp_set_dynamic(0);

    for (unsigned long t = 0; t < nblocksteps; t++) {
        // arrays for intermediate results to split up molecule loop
        // loop over molecules

#pragma omp parallel for
        for (size_t i = 0; i < nmols; i++) {
            // convenience variables
            size_t m_firstatom = mol_firstatom[i];
            size_t m_natoms = mol_natoms[i];
            size_t m_moltype = mol_moltypenr[i];
            float m_mass = mol_mass[i];
            float *m_atommasses = moltypes_atommasses[m_moltype];
            int *m_abc_indicators = moltype_abc_indicators[m_moltype];
            char m_rot_treat = moltype_rot_treat[m_moltype];

            // large arrays
            // TODO: Better would be to iterate over moltypes and allocate those
            // outside molecule loop
            float *positions = calloc(3 * m_natoms, sizeof(float));
            float *velocities = calloc(3 * m_natoms, sizeof(float));
            float *positions_rel = calloc(3 * m_natoms, sizeof(float));
            float *velocities_rot = calloc(3 * m_natoms, sizeof(float));

            // reading into threadprivate arrays
            for (size_t j = 0; j < m_natoms; j++) {
                size_t jj = m_firstatom + j;
                positions[3 * j + 0] = block_pos[3 * natoms * t + 3 * jj + 0];
                positions[3 * j + 1] = block_pos[3 * natoms * t + 3 * jj + 1];
                positions[3 * j + 2] = block_pos[3 * natoms * t + 3 * jj + 2];
                velocities[3 * j + 0] = block_vel[3 * natoms * t + 3 * jj + 0];
                velocities[3 * j + 1] = block_vel[3 * natoms * t + 3 * jj + 1];
                velocities[3 * j + 2] = block_vel[3 * natoms * t + 3 * jj + 2];
            }

            // recombination
            if (no_pbc == false) {
                recombine_molecule(&block_box[3 * t], m_natoms, positions);
            }

            // single atoms
            if (m_natoms == 1) {
                size_t atom = m_firstatom;
                for (size_t dim = 0; dim < 3; dim++) {
                    mol_velocities_sqrt_m_trn[3 * nblocksteps * i +
                                              nblocksteps * dim + t] =
                        velocities[dim] * sqrt(m_mass);
                    mol_omegas_sqrt_i_rot[3 * nblocksteps * i +
                                          nblocksteps * dim + t] = 0;
                    atom_velocities_sqrt_m_vib[3 * nblocksteps * atom +
                                               nblocksteps * dim + t] = 0;
                    atom_velocities_sqrt_m_rot[3 * nblocksteps * atom +
                                               nblocksteps * dim + t] = 0;
                }
                continue;
            }

            // unseperated dos -> output to dos_vib
            if (m_rot_treat == 'u') {
                for (size_t j = 0; j < m_natoms; j++) {
                    size_t atom = m_firstatom + j;
                    for (size_t dim = 0; dim < 3; dim++) {
                        mol_velocities_sqrt_m_trn[3 * nblocksteps * i +
                                                  nblocksteps * dim + t] = 0;
                        mol_omegas_sqrt_i_rot[3 * nblocksteps * i +
                                              nblocksteps * dim + t] = 0;
                        atom_velocities_sqrt_m_vib[3 * nblocksteps * atom +
                                                   nblocksteps * dim + t] =
                            velocities[3 * j + dim] * sqrt(m_atommasses[j]);
                        atom_velocities_sqrt_m_rot[3 * nblocksteps * atom +
                                                   nblocksteps * dim + t] = 0;
                    }
                }
                continue;
            }

            // calc molecule velocity and molecule com
            float center_of_mass[3] = {0.0, 0.0, 0.0};
            float mol_velocity_trn[3] = {0.0, 0.0, 0.0};
            // this might also be done by a matrix multiplication but maybe
            // overkill
            for (size_t dim = 0; dim < 3; dim++) {
                mol_velocity_trn[dim] = cblas_sdot(m_natoms, m_atommasses, 1,
                                                   &velocities[0 + dim], 3);
                mol_velocity_trn[dim] /= m_mass;

                // output-array mol_velocities_sqrt_m_trn
                mol_velocities_sqrt_m_trn[3 * nblocksteps * i +
                                          nblocksteps * dim + t] =
                    mol_velocity_trn[dim] * sqrt(m_mass);

                center_of_mass[dim] = cblas_sdot(m_natoms, m_atommasses, 1,
                                                 &positions[0 + dim], 3);
                center_of_mass[dim] /= m_mass;
            }

            // calc molecule atoms relative positions
            for (size_t j = 0; j < m_natoms; j++) {
                for (size_t dim = 0; dim < 3; dim++) {
                    positions_rel[3 * j + dim] =
                        positions[3 * j + dim] - center_of_mass[dim];
                }
            }

            // calc angular momentum
            float cross_product[3] = {0.0, 0.0, 0.0};
            float angular_momentum[3] = {0.0, 0.0, 0.0};
            for (size_t j = 0; j < m_natoms; j++) {
                crossProduct(&positions_rel[3 * j], &velocities[3 * j],
                             cross_product);
                for (size_t dim = 0; dim < 3; dim++) {
                    angular_momentum[dim] +=
                        m_atommasses[j] * cross_product[dim];
                }
            }

            // calc moi tensor
            float moi_tensor[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
            moiTensor(m_natoms, positions_rel, m_atommasses, moi_tensor);

            // calc angular velocity
            float angular_velocity[3] = {0.0, 0.0, 0.0};
            float moi_tensor_temp[9] = {0.0, 0.0, 0.0, 0.0, 0.0,
                                        0.0, 0.0, 0.0, 0.0};
            cblas_scopy(9, moi_tensor, 1, moi_tensor_temp, 1);
            cblas_scopy(3, angular_momentum, 1, angular_velocity, 1);
            // linear molecules
            if (m_rot_treat == 'l') {
                // find angular velocity from underdetermined system of linear
                // equations
                float S[3] = {0.0, 0.0, 0.0};
                int rank = 1337;

                LAPACKE_sgelsd(LAPACK_ROW_MAJOR, 3, 3, 1, moi_tensor_temp, 3,
                               angular_velocity, 1, S, 0.001, &rank);
            }
            // non-linear molecules
            else {
                // find angular velocity from system of linear equations
                int ipiv[3] = {0, 0, 0};
                LAPACKE_sgesv(LAPACK_ROW_MAJOR, 3, 1, moi_tensor_temp, 3, ipiv,
                              angular_velocity, 1);
            }

            // calc velocities rot
            for (size_t j = 0; j < m_natoms; j++) {
                crossProduct(angular_velocity, &positions_rel[3 * j],
                             cross_product);
                for (size_t dim = 0; dim < 3; dim++) {
                    velocities_rot[3 * j + dim] = cross_product[dim];
                }
            }

            // calc velocities vib
            float velocity_vib[3] = {0.0, 0.0, 0.0};
            for (size_t j = 0; j < m_natoms; j++) {
                cblas_scopy(3, &velocities[3 * j], 1, velocity_vib, 1);

                cblas_saxpy(3, -1.0, &velocities_rot[3 * j], 1, velocity_vib,
                            1);
                cblas_saxpy(3, -1.0, mol_velocity_trn, 1, velocity_vib, 1);

                // write in output array atom_velocities_sqrt_m_vib
                size_t atom = m_firstatom + j;
                for (size_t dim = 0; dim < 3; dim++) {
                    atom_velocities_sqrt_m_vib[3 * nblocksteps * atom +
                                               nblocksteps * dim + t] =
                        velocity_vib[dim] * sqrt(m_atommasses[j]);
                    // DoS_rot_xyz
                    atom_velocities_sqrt_m_rot[3 * nblocksteps * atom +
                                               nblocksteps * dim + t] =
                        velocities_rot[3 * j + dim] * sqrt(m_atommasses[j]);
                }
            }

            // linear molecules
            if (m_rot_treat == 'l') {
                // calculate series for FFT
                for (size_t dim = 0; dim < 3; dim++) {
                    /*
                       if (copysignf(1.0, angular_velocity[dim]) !=
                       copysignf(1.0, angular_momentum[dim]))
                       {
                       fprintf(stderr, "ERROR: angular velocity and
                       angular_momentum don't point in the same direction!\n");
                       exit(1);
                       }
                       */
                    mol_omegas_sqrt_i_rot[3 * nblocksteps * i +
                                          nblocksteps * dim + t] =
                        copysignf(sqrt_neg_zero(angular_velocity[dim] *
                                                angular_momentum[dim]),
                                  angular_velocity[dim]);
                }
                // save moments of inertia for each molecule
                {
                    for (size_t dim = 0; dim < 3; dim++) {
                        mol_block_moments_of_inertia[3 * nblocksteps * i +
                                                     nblocksteps * dim + t] +=
                            moi_tensor[3 * dim + dim];
                        mol_block_moments_of_inertia_squared
                            [3 * nblocksteps * i + nblocksteps * dim + t] +=
                            powf(moi_tensor[3 * dim + dim], 2.0);
                    }
                }
                continue;
            }

            // calc moments of inertia and eigenvectors
            float eigenvectors[9] = {0.0, 0.0, 0.0, 0.0, 0.0,
                                     0.0, 0.0, 0.0, 0.0};
            cblas_scopy(9, moi_tensor, 1, eigenvectors, 1);
            float moments_of_inertia[3] = {0.0, 0.0, 0.0};
            if (LAPACKE_ssyev(LAPACK_ROW_MAJOR, 'V', 'U', 3, eigenvectors, 3,
                              moments_of_inertia) > 0) {
                fprintf(
                    stderr,
                    "ERROR: LAPACKE_ssyev failed to compute eigenvalues of\n");
                fprintf(stderr,
                        "       the moment of inertia tensor of molecule %zu\n",
                        i);
                exit(1);
            }

            // abc auxilary vectors
            float a[3] = {0.0, 0.0, 0.0};
            float b[3] = {0.0, 0.0, 0.0};
            float c[3] = {0.0, 0.0, 0.0};

            // about the abc_indicators
            // two numbers define a, two define b', c is cross product of a and
            // b', b is cross product of c and a vector from atom(first number)
            // to atom(second number) second number can be -1 for com c is
            // always the cross product

            for (size_t dim = 0; dim < 3; dim++) {
                a[dim] = positions_rel[3 * m_abc_indicators[0] + dim];
                if (m_abc_indicators[1] != -1)
                    a[dim] -= positions_rel[3 * m_abc_indicators[1] + dim];

                b[dim] = positions_rel[3 * m_abc_indicators[2] + dim];
                if (m_abc_indicators[3] != -1)
                    b[dim] -= positions_rel[3 * m_abc_indicators[3] + dim];
            }
            // normalize a
            cblas_sscal(3, 1 / cblas_snrm2(3, &a[0], 1), &a[0], 1);

            // calculate c and b
            crossProduct(a, b, c);
            crossProduct(c, a, b);

            // normalize b and c
            cblas_sscal(3, 1 / cblas_snrm2(3, &b[0], 1), &b[0], 1);
            cblas_sscal(3, 1 / cblas_snrm2(3, &c[0], 1), &c[0], 1);

            // check if eigenvector points in same general direction as abc
            // if not flip eigenvector
            if (cblas_sdot(3, &eigenvectors[0], 3, a, 1) < 0.0)
                cblas_sscal(3, -1.0, &eigenvectors[0], 3);
            if (cblas_sdot(3, &eigenvectors[1], 3, b, 1) < 0.0)
                cblas_sscal(3, -1.0, &eigenvectors[1], 3);
            if (cblas_sdot(3, &eigenvectors[2], 3, c, 1) < 0.0)
                cblas_sscal(3, -1.0, &eigenvectors[2], 3);

            // angular velocity for output later
            // can be full ω or Eckart Ω
            float output_angular_velocity[3];
            // cases where it is ω
            if (m_rot_treat == 'f' || m_rot_treat == 'a') {
                cblas_scopy(3, angular_velocity, 1, output_angular_velocity, 1);
            }

            // eckart frame decomposition
            if (m_rot_treat == 'e' || m_rot_treat == 'p' ||
                m_rot_treat == 'E' || m_rot_treat == 'P') {
                // calculate Eckart vectors F (F_i in rows)
                float F[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                for (size_t j = 0; j < m_natoms; j++) {
                    for (size_t dim1 = 0; dim1 < 3; dim1++) {
                        for (size_t dim2 = 0; dim2 < 3; dim2++) {
                            F[3 * dim1 + dim2] +=
                                m_atommasses[j] *
                                atom_refpos_principal_components
                                    [m_firstatom * 3 + m_natoms * dim1 + j] *
                                positions[3 * j + dim2];
                        }
                    }
                }
                // calculate Eckart frame f
                // f has f1,f2,f3 in columns
                float f[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                // planar molecule
                if (m_rot_treat == 'p' || m_rot_treat == 'P') {
                    // gram_matrix = F12 @ F12.T
                    float gram_matrix[4] = {0.0, 0.0, 0.0, 0.0};
                    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 2, 2,
                                3, 1.0, F, 3, F, 3, 0.0, gram_matrix, 2);
                    // invert gram matrix
                    float gram_matrix_inv[4] = {0.0, 0.0, 0.0, 0.0};
                    cblas_scopy(4, gram_matrix, 1, gram_matrix_inv, 1);
                    if (invert_matrix(gram_matrix_inv, 2) != 0) {
                        fprintf(stderr,
                                "ERROR: non-invertible Gram matrix of molecule "
                                "%zu\n",
                                i);
                        exit(1);
                    }
                    // squareroot of inverse gram matrix
                    float gram_matrix_inv_sqrt[4] = {0.0, 0.0, 0.0, 0.0};
                    cblas_scopy(4, gram_matrix_inv, 1, gram_matrix_inv_sqrt, 1);
                    if (squareroot_of_matrix(gram_matrix_inv_sqrt, 2)) {
                        fprintf(stderr,
                                "ERROR: could not take square root of the\n");
                        fprintf(stderr,
                                "       inverse gram matrix of molecule %zu\n",
                                i);
                        exit(1);
                    }
                    // Eckart frame f1 f2
                    cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 3, 2,
                                2, 1.0, F, 3, gram_matrix_inv_sqrt, 2, 0.0, f,
                                3);
                    // f3 = f1 x f2
                    float f1[3] = {f[0], f[3], f[6]};
                    float f2[3] = {f[1], f[4], f[7]};
                    float f3[3] = {0.0, 0.0, 0.0};
                    crossProduct(f1, f2, f3);
                    cblas_scopy(3, f3, 1, &f[2], 3);
                }
                // non-planar molecule
                else if (m_rot_treat == 'e' || m_rot_treat == 'E') {
                    // gram_matrix = F12 @ F12.T
                    float gram_matrix[9] = {0.0, 0.0, 0.0, 0.0, 0.0,
                                            0.0, 0.0, 0.0, 0.0};
                    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 3, 3,
                                3, 1.0, F, 3, F, 3, 0.0, gram_matrix, 3);
                    // invert gram matrix
                    float gram_matrix_inv[9] = {0.0, 0.0, 0.0, 0.0, 0.0,
                                                0.0, 0.0, 0.0, 0.0};
                    cblas_scopy(9, gram_matrix, 1, gram_matrix_inv, 1);
                    if (invert_matrix(gram_matrix_inv, 3) != 0) {
                        fprintf(stderr,
                                "ERROR: non-invertible Gram matrix of molecule "
                                "%zu\n",
                                i);
                        exit(1);
                    }
                    // squareroot of inverse gram matrix
                    float gram_matrix_inv_sqrt[9] = {0.0, 0.0, 0.0, 0.0, 0.0,
                                                     0.0, 0.0, 0.0, 0.0};

                    cblas_scopy(9, gram_matrix_inv, 1, gram_matrix_inv_sqrt, 1);
                    if (squareroot_of_matrix(gram_matrix_inv_sqrt, 3)) {
                        fprintf(stderr,
                                "ERROR: could not take square root of the\n");
                        fprintf(stderr,
                                "       inverse gram matrix of molecule %zu\n",
                                i);
                        exit(1);
                    }
                    // Eckart frame
                    cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 3, 3,
                                3, 1.0, F, 3, gram_matrix_inv_sqrt, 3, 0.0, f,
                                3);
                }
                // check F, f from Louck et al. Σ_i F_i x f_i = 0
                float check_louck[3] = {0.0, 0.0, 0.0};
                float cross_product[3] = {0.0, 0.0, 0.0};
                for (size_t dim = 0; dim < 3; dim++) {
                    crossProductInc(&F[3 * dim], 1, &f[dim], 3, cross_product,
                                    1);
                    cblas_saxpy(3, 1.0, cross_product, 1, check_louck, 1);
                }
                // Not sure what number is reasonable here
                if (fabs(check_louck[0]) > 1.0e-3 ||
                    fabs(check_louck[2]) > 1.0e-3 ||
                    fabs(check_louck[2]) > 1.0e-3) {
                    fprintf(
                        stderr,
                        "ERROR: Check on Eckart frame and Eckart vectors\n");
                    fprintf(stderr, "       failed for molecule %zu\n", i);
                    exit(1);
                }
                // positions of reference in lab frame, one atom per row
                // sablic (8)
                float *c_alpha = calloc(m_natoms * 3, sizeof(float));
                for (size_t j = 0; j < m_natoms; j++) {
                    for (size_t dim1 = 0; dim1 < 3; dim1++) {
                        for (size_t dim2 = 0; dim2 < 3; dim2++) {
                            c_alpha[3 * j + dim2] +=
                                atom_refpos_principal_components
                                    [m_firstatom * 3 + m_natoms * dim1 + j] *
                                f[3 * dim2 + dim1];
                        }
                    }
                }
                // possible TODO: calculate rho and check second Eckart
                // condition

                // calculate J' tensor
                float J_prime[9] = {0.0, 0.0, 0.0, 0.0, 0.0,
                                    0.0, 0.0, 0.0, 0.0};
                moiTensorMixed(m_natoms, positions_rel, c_alpha, m_atommasses,
                               J_prime);
                // Eckart angular momentum
                float eckart_angular_momentum[3] = {0.0, 0.0, 0.0};
                for (size_t j = 0; j < m_natoms; j++) {
                    crossProduct(&c_alpha[3 * j], &velocities[3 * j],
                                 cross_product);
                    for (size_t dim = 0; dim < 3; dim++) {
                        eckart_angular_momentum[dim] +=
                            m_atommasses[j] * cross_product[dim];
                    }
                }
                // Eckart angular velocity Ω
                float eckart_angular_velocity[3] = {0.0, 0.0, 0.0};
                // J'
                float J_prime_temp[9] = {0.0, 0.0, 0.0, 0.0, 0.0,
                                         0.0, 0.0, 0.0, 0.0};
                cblas_scopy(9, J_prime, 1, J_prime_temp, 1);
                cblas_scopy(3, eckart_angular_momentum, 1,
                            eckart_angular_velocity, 1);
                // find Ω from system of linear equations
                int ipiv[3] = {0, 0, 0};
                LAPACKE_sgesv(LAPACK_ROW_MAJOR, 3, 1, J_prime_temp, 3, ipiv,
                              eckart_angular_velocity, 1);
                // save Ω for later
                cblas_scopy(3, eckart_angular_velocity, 1,
                            output_angular_velocity, 1);
                // total angular velocity ω - Eckart angular velocity Ω
                float omega_minus_Omega[3] = {0.0, 0.0, 0.0};
                cblas_scopy(3, angular_velocity, 1, omega_minus_Omega, 1);
                cblas_saxpy(3, -1.0, eckart_angular_velocity, 1,
                            omega_minus_Omega, 1);
                // vibrational motion coupled with rotation u
                // u_j = (ω - Ω) x δr_j
                // and Coriolis energy term
                // Σ_j u_j · (Ω x δr_j)
                float velocity_vibc[3] = {0.0, 0.0, 0.0};
                for (size_t j = 0; j < m_natoms; j++) {
                    size_t atom = m_firstatom + j;
                    // vibc
                    crossProduct(omega_minus_Omega, &positions_rel[3 * j],
                                 velocity_vibc);
                    // cross product for coriolis
                    crossProduct(eckart_angular_velocity, &positions_rel[3 * j],
                                 cross_product);
                    // add coriolis energy
                    mol_block_coriolis[nblocksteps * i + t] +=
                        m_atommasses[j] *
                        cblas_sdot(3, velocity_vibc, 1, cross_product, 1);
                    // write in output array atom_velocities_sqrt_m_vibc
                    for (size_t dim = 0; dim < 3; dim++) {
                        atom_velocities_sqrt_m_vibc[3 * nblocksteps * atom +
                                                    nblocksteps * dim + t] =
                            velocity_vibc[dim] * sqrt(m_atommasses[j]);
                    }
                }
                // writing in output arrays
                /*
                for (size_t dim = 0; dim < 3; dim++) {
                    // the rotational velocities
                    mol_omegas_sqrt_i_rot[3 * nblocksteps * i +
                                          nblocksteps * dim + t] =
                        eckart_angular_velocity_pa[dim] *
                        sqrt(moments_of_inertia[dim]);
                    // save moments of inertia for each molecule
                    mol_block_moments_of_inertia[3 * nblocksteps * i +
                                                 nblocksteps * dim + t] +=
                        moments_of_inertia[dim];
                    mol_block_moments_of_inertia_squared[3 * nblocksteps * i +
                                                         nblocksteps * dim +
                                                         t] +=
                        powf(moments_of_inertia[dim], 2.0);
                }
                */

                free(c_alpha);
            }

            // using auxillary frame decompostion
            if (m_rot_treat == 'a' || m_rot_treat == 'E' ||
                m_rot_treat == 'P') {
                // matrix abc which is used for transfroming in aux frame
                float abc[9];
                cblas_scopy(3, a, 1, &abc[0], 3);
                cblas_scopy(3, b, 1, &abc[1], 3);
                cblas_scopy(3, c, 1, &abc[2], 3);

                // calculate angular velocity in abc coords
                float angular_velocity_abc[3] = {0.0, 0.0, 0.0};
                // angular_velocity_abc = angular_velocity @ abc
                for (size_t dim = 0; dim < 3; dim++) {
                    angular_velocity_abc[dim] =
                        cblas_sdot(3, output_angular_velocity, 1, &abc[dim], 3);
                }

                // calc moi_tensor in abc coords
                float temp_matrix[9] = {0.0, 0.0, 0.0, 0.0, 0.0,
                                        0.0, 0.0, 0.0, 0.0};
                float abc_inv[9] = {0.0, 0.0, 0.0, 0.0, 0.0,
                                    0.0, 0.0, 0.0, 0.0};
                cblas_scopy(9, abc, 1, abc_inv, 1);
                if (invert_matrix(abc_inv, 3) != 0) {
                    fprintf(stderr,
                            "ERROR: non-invertible abc matrix of molecule "
                            "%zu\n",
                            i);
                    exit(1);
                }
                // temp_matrix = abc_inv @ moi_tensor
                cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3,
                            1.0, abc_inv, 3, moi_tensor, 3, 0.0, temp_matrix,
                            3);

                // moi_abc = temp_matrix @ abc
                float moi_tensor_abc[9] = {0.0, 0.0, 0.0, 0.0, 0.0,
                                           0.0, 0.0, 0.0, 0.0};
                cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3,
                            1.0, temp_matrix, 3, abc, 3, 0.0, moi_tensor_abc,
                            3);

                // writing in output arrays the rotational velocities
                for (size_t dim = 0; dim < 3; dim++) {
                    // 'a'bc as rotational axis
                    mol_omegas_sqrt_i_rot[3 * nblocksteps * i +
                                          nblocksteps * dim + t] =
                        angular_velocity_abc[dim] *
                        sqrt(moi_tensor_abc[3 * dim + dim]);
                    // save moments of inertia for each molecule
                    mol_block_moments_of_inertia[3 * nblocksteps * i +
                                                 nblocksteps * dim + t] +=
                        moi_tensor_abc[3 * dim + dim];
                    mol_block_moments_of_inertia_squared[3 * nblocksteps * i +
                                                         nblocksteps * dim +
                                                         t] +=
                        powf(moi_tensor_abc[3 * dim + dim], 2.0);
                }
                continue;
            }

            // using principal axi frame decompostion
            if (m_rot_treat == 'f' || m_rot_treat == 'e' ||
                m_rot_treat == 'p') {
                // calculate angular velocity in pa coords
                float angular_velocity_pa[3] = {0.0, 0.0, 0.0};
                for (size_t dim = 0; dim < 3; dim++) {
                    angular_velocity_pa[dim] = cblas_sdot(
                        3, output_angular_velocity, 1, &eigenvectors[dim], 3);
                }
                // writing in output arrays
                for (size_t dim = 0; dim < 3; dim++) {
                    // the rotational velocities
                    mol_omegas_sqrt_i_rot[3 * nblocksteps * i +
                                          nblocksteps * dim + t] =
                        angular_velocity_pa[dim] *
                        sqrt(moments_of_inertia[dim]);
                    // save moments of inertia for each molecule
                    mol_block_moments_of_inertia[3 * nblocksteps * i +
                                                 nblocksteps * dim + t] +=
                        moments_of_inertia[dim];
                    mol_block_moments_of_inertia_squared[3 * nblocksteps * i +
                                                         nblocksteps * dim +
                                                         t] +=
                        powf(moments_of_inertia[dim], 2.0);
                }
            }

            free(velocities);
            free(positions);
            free(positions_rel);
            free(velocities_rot);
        }
    }
}
