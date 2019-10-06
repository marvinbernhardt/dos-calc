#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <cblas.h>
#include <lapacke.h>
#include "xdrfile.h"
#include "xdrfile_trr.h"
#include "linear-algebra.c"

#ifdef DEBUG
#define DPRINT(...) do{ fprintf( stdout, __VA_ARGS__ ); } while( 0 )
#else
#define DPRINT(...)
#endif


float sqrt_neg_zero(float number)
{
    if (number < -0.001)
    {
        printf("angular_momentum and angular_velocity too different!\n");
        exit(0);
    }
    else if (number < 0)
        return 0;
    else
        return sqrt(number);
}

void decomposeVelocities (XDRFILE *traj,
        unsigned long nblocksteps,
        size_t natoms_traj,
        size_t nmols,
        size_t nmoltypes,
        size_t *mol_firstatom,
        size_t *mol_natoms,
        size_t *mol_moltypenr,
        float **moltypes_atommasses,
        float *mol_mass,
        size_t *moltype_natomspermol,
        char *moltype_rot_treat,
        int **moltype_abc_indicators,
        bool no_pbc,
        float *mol_velocities_sqrt_m_trn,  // from here output
        float *mol_omegas_sqrt_i_rot,
        float *atom_velocities_sqrt_m_vib,
        float *atom_velocities_sqrt_m_rot,
        float *mol_moments_of_inertia)
{
    // for reading of frame
    int step;
    float time;
    float lambda;
    matrix box;
    rvec *r;
    rvec *v;
    int has_prop;
    r = calloc(natoms_traj, sizeof(*r));
    v = calloc(natoms_traj, sizeof(*v));

    // no dynamic teams!!!
    omp_set_dynamic(0);

    // per-molecule arrays
    static float *positions;
    static float *velocities;
    static float *positions_rel;
    static float *velocities_rot;
    #pragma omp threadprivate(positions, velocities, positions_rel, velocities_rot)

    // largest molecule n_atoms
    size_t mol_natoms_max = 0;
    for (size_t h=0; h<nmoltypes; h++)
    {
        if (moltype_natomspermol[h] > mol_natoms_max)
        {
            mol_natoms_max = moltype_natomspermol[h];
        }
    }

    #pragma omp parallel
    {
        // allocate positions, velocities (work arrays private to each thread)
        positions = calloc(3*mol_natoms_max, sizeof(float));
        velocities = calloc(3*mol_natoms_max, sizeof(float));
        positions_rel = calloc(3*mol_natoms_max, sizeof(float));
        velocities_rot = calloc(3*mol_natoms_max, sizeof(float));
    }

    DPRINT("start reading frame\n");
    for (unsigned long t=0; t<nblocksteps; t++)
    {
        if (read_trr(traj, natoms_traj, &step, &time, &lambda, box, r, v, NULL, &has_prop) != 0)
        {
            fprintf(stderr, "ERROR: Reading frame %lu failed\n", t);
            exit(1);
        }

        DPRINT("There are %zu atoms at step %lu (time %f). My box is: %f %f %f \n",
                natoms_traj, t, time, box[0][0], box[1][1], box[2][2]);

        // loop over molecules
        #pragma omp parallel for
        for (size_t i=0; i<nmols; i++)
        {
            DPRINT("\ndoing molecule %zu\n", i);
            size_t m_firstatom = mol_firstatom[i];
            size_t m_natoms = mol_natoms[i];
            size_t m_moltype = mol_moltypenr[i];
            float m_mass = mol_mass[i];
            float *m_atommasses = moltypes_atommasses[m_moltype];
            int *m_abc_indicators = moltype_abc_indicators[m_moltype];
            char m_rot_treat = moltype_rot_treat[m_moltype];

            //recombination
            if (no_pbc==false)
            {
                for(size_t dim=0; dim < 3; dim++)
                {
                    for(size_t j=1; j<m_natoms; j++)
                    {
                        size_t jj = m_firstatom + j;
                        float dist_to_firstatom = r[jj][dim] - r[m_firstatom][dim];
                        if(dist_to_firstatom > 0.5 * box[dim][dim])
                        {
                            r[jj][dim] -= box[dim][dim];
                        }
                        if(dist_to_firstatom < -0.5 * box[dim][dim])
                        {
                            r[jj][dim] += box[dim][dim];
                        }
                    }
                }
            }

            // read atoms of one molecule
            for (size_t j=0; j<m_natoms; j++)
            {
                size_t jj = m_firstatom + j;
                DPRINT("scanning atom %zu (nr. %zu in molecule)\n", jj, j);

                positions[3*j+0] = r[jj][0];
                positions[3*j+1] = r[jj][1];
                positions[3*j+2] = r[jj][2];
                velocities[3*j+0] = v[jj][0];
                velocities[3*j+1] = v[jj][1];
                velocities[3*j+2] = v[jj][2];

                DPRINT("%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n",
                        positions[3*j+0], positions[3*j+1], positions[3*j+2],
                        velocities[3*j+0], velocities[3*j+1], velocities[3*j+2]);
            }

            // single atoms
            // ============
            if (m_natoms == 1)
            {
                size_t atom = m_firstatom;
                for (size_t dim=0; dim<3; dim++)
                {
                    mol_velocities_sqrt_m_trn[3*nblocksteps*i + nblocksteps*dim + t] = velocities[dim] * sqrt(m_mass);
                    mol_omegas_sqrt_i_rot[3*nblocksteps*i + nblocksteps*dim + t] = 0;
                    atom_velocities_sqrt_m_vib[3*nblocksteps*atom + nblocksteps*dim + t] = 0;
                    atom_velocities_sqrt_m_rot[3*nblocksteps*atom + nblocksteps*dim + t] = 0;
                }
                continue;
            }

            // single atoms
            // unseperated dos -> output to dos_vib
            if (m_rot_treat == 'u')
            {
                for (size_t j=0; j<m_natoms; j++)
                {
                    size_t atom = m_firstatom + j;
                    for (size_t dim=0; dim<3; dim++)
                    {
                        mol_velocities_sqrt_m_trn[3*nblocksteps*i + nblocksteps*dim + t] = 0;
                        mol_omegas_sqrt_i_rot[3*nblocksteps*i + nblocksteps*dim + t] = 0;
                        atom_velocities_sqrt_m_vib[3*nblocksteps*atom + nblocksteps*dim + t] = velocities[3*j+dim] * sqrt(m_atommasses[j]);
                        atom_velocities_sqrt_m_rot[3*nblocksteps*atom + nblocksteps*dim + t] = 0;
                    }
                }
                continue;
            }

            // calc molecule velocity and molecule com
            float center_of_mass[3] = {0.0, 0.0, 0.0};
            float mol_velocity_trn[3] = {0.0, 0.0, 0.0};
            // this might also be done by a matrix multiplication but maybe overkill
            for (size_t dim=0; dim<3; dim++)
            {
                mol_velocity_trn[dim] = cblas_sdot(m_natoms,
                        m_atommasses, 1,
                        &velocities[0+dim], 3);
                mol_velocity_trn[dim] /= m_mass;

                // output-array mol_velocities_sqrt_m_trn
                mol_velocities_sqrt_m_trn[3*nblocksteps*i + nblocksteps*dim + t] = mol_velocity_trn[dim] * sqrt(m_mass);

                center_of_mass[dim] = cblas_sdot(m_natoms,
                        m_atommasses, 1,
                        &positions[0+dim], 3);
                center_of_mass[dim] /= m_mass;
            }

            DPRINT("mol_velocity_trn: %8.4f%8.4f%8.4f\n",
                    mol_velocity_trn[0],
                    mol_velocity_trn[1],
                    mol_velocity_trn[2]);
            DPRINT("com: %8.4f%8.4f%8.4f\n", center_of_mass[0], center_of_mass[1], center_of_mass[2]);


            // calc molecule atoms relative positions
            for (size_t j=0; j<m_natoms; j++)
            {
                for (size_t dim=0; dim<3; dim++)
                {
                    positions_rel[3*j+dim] = positions[3*j+dim] - center_of_mass[dim];
                }

                DPRINT("pos_rel atom %zu: %8.4f%8.4f%8.4f\n", j,
                        positions_rel[3*j+0], positions_rel[3*j+1], positions_rel[3*j+2]);
            }

            // calc angular momentum
            float cross_product[3] = {0.0, 0.0, 0.0};
            float angular_momentum[3] = {0.0, 0.0, 0.0};
            for (size_t j=0; j<m_natoms; j++)
            {
                crossProduct(&positions_rel[3*j], &velocities[3*j], cross_product);
                for (size_t dim=0; dim<3; dim++)
                {
                    angular_momentum[dim] += m_atommasses[j] * cross_product[dim];
                }
            }

            DPRINT("angular momentum: %8.4f%8.4f%8.4f\n",
                    angular_momentum[0], angular_momentum[1], angular_momentum[2]);


            // calc moi tensor
            float moi_tensor[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
            moiTensor(m_natoms, positions_rel, m_atommasses, moi_tensor);

            DPRINT("moi Tensor:\n");
            DPRINT("%f %f %f\n", moi_tensor[0], moi_tensor[1], moi_tensor[2]);
            DPRINT("%f %f %f\n", moi_tensor[3], moi_tensor[4], moi_tensor[5]);
            DPRINT("%f %f %f\n", moi_tensor[6], moi_tensor[7], moi_tensor[8]);


            // calc angular velocity
            float angular_velocity[3] = {0.0, 0.0, 0.0};
            float moi_tensor_temp[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
            cblas_scopy(9, moi_tensor, 1, moi_tensor_temp, 1);
            cblas_scopy(3, angular_momentum, 1, angular_velocity, 1);

            // linear molecules
            if (m_rot_treat == 'l')
            {
                // find angular velocity from underdetermined system of linear equations
                float S[3] = {0.0, 0.0, 0.0};
                int rank = 1337;

                LAPACKE_sgelsd(LAPACK_ROW_MAJOR, 3, 3, 1, moi_tensor_temp, 3, angular_velocity, 1, S, 0.001, &rank);
            }
            // non-linear molecules
            else
            {
                // find angular velocity from system of linear equations
                int ipiv[3] = {0, 0, 0};
                LAPACKE_sgesv(LAPACK_ROW_MAJOR, 3, 1, moi_tensor_temp, 3, ipiv, angular_velocity, 1);
            }

            DPRINT("angular velocity: %f %f %f\n",
                    angular_velocity[0], angular_velocity[1], angular_velocity[2]);


            // calc velocities rot
            for (size_t j=0; j<m_natoms; j++)
            {
                crossProduct(angular_velocity, &positions_rel[3*j], cross_product);
                for (size_t dim=0; dim<3; dim++)
                {
                    velocities_rot[3*j+dim] = cross_product[dim];
                }

                DPRINT("velocity_rot atom %zu: %8.4f%8.4f%8.4f\n", j,
                        velocities_rot[3*j+0], velocities_rot[3*j+1], velocities_rot[3*j+2]);
            }

            // calc velocities vib
            float velocity_vib[3] = {0.0, 0.0, 0.0};
            for (size_t j=0; j<m_natoms; j++)
            {
                cblas_scopy(3, &velocities[3*j], 1, velocity_vib, 1);

                cblas_saxpy(3, -1.0, &velocities_rot[3*j], 1,
                        velocity_vib, 1);
                cblas_saxpy(3, -1.0, mol_velocity_trn, 1,
                        velocity_vib, 1);

                // output-array atom_velocities_sqrt_m_vib
                size_t atom = m_firstatom + j;
                for (size_t dim=0; dim<3; dim++)
                {
                    atom_velocities_sqrt_m_vib[3*nblocksteps*atom + nblocksteps*dim + t] = velocity_vib[dim] * sqrt(m_atommasses[j]);
                    // alternative DoS_rot
                    atom_velocities_sqrt_m_rot[3*nblocksteps*atom + nblocksteps*dim + t] = velocities_rot[3*j+dim] * sqrt(m_atommasses[j]);
                }

                DPRINT("velocity_vib atom %zu: %8.4f%8.4f%8.4f\n", j,
                        velocity_vib[0],
                        velocity_vib[1],
                        velocity_vib[2]);
            }

            // linear molecules
            if (m_rot_treat == 'l')
            {
                for (size_t dim=0; dim<3; dim++)
                {
                    /*
                    if (copysignf(1.0, angular_velocity[dim]) != copysignf(1.0, angular_momentum[dim]))
                    {
                        fprintf(stderr, "ERROR: angular velocity and angular_momentum don't point in the same direction!\n");
                        exit(1);
                    }
                    */

                    mol_omegas_sqrt_i_rot[3*nblocksteps*i + nblocksteps*dim + t] = copysignf(sqrt_neg_zero(angular_velocity[dim] * angular_momentum[dim]), angular_velocity[dim]);
                }

                // save moments of inertia for each molecule
                for (size_t dim=0; dim<3; dim++)
                {
                    mol_moments_of_inertia[3*i+dim] += moi_tensor[3*dim+dim];
                }

                continue;
            }

            // calc moments of inertia and eigenvectors
            float eigenvectors[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
            cblas_scopy(9, moi_tensor, 1, eigenvectors, 1);
            float moments_of_inertia[3] = {0.0, 0.0, 0.0};
            if (LAPACKE_ssyev(LAPACK_ROW_MAJOR, 'V', 'U',
                        3, eigenvectors, 3, moments_of_inertia) > 0)
            {
                fprintf(stderr, "ERROR: LAPACKE_ssyev failed to compute eigenvalues of\n");
                fprintf(stderr, "       the moment of inertia tensor of molecule %zu\n", i);
                exit(1);
            }

            DPRINT("moments of inertia: %8.4f%8.4f%8.4f\n",
                    moments_of_inertia[0], moments_of_inertia[1], moments_of_inertia[2]);

            // save moments of inertia for each molecule
            for (size_t dim=0; dim<3; dim++)
            {
                mol_moments_of_inertia[3*i+dim] += moments_of_inertia[dim];
            }

            DPRINT("eigenvectors (in columns):\n");
            DPRINT("%f %f %f\n", eigenvectors[0], eigenvectors[1], eigenvectors[2]);
            DPRINT("%f %f %f\n", eigenvectors[3], eigenvectors[4], eigenvectors[5]);
            DPRINT("%f %f %f\n", eigenvectors[6], eigenvectors[7], eigenvectors[8]);

            // abc auxilary vectors
            float a[3] = {0.0, 0.0, 0.0};
            float b[3] = {0.0, 0.0, 0.0};
            float c[3] = {0.0, 0.0, 0.0};

            // about the abc_indicators
            // two numbers define a, two define b', c is cross product of a and b', b is cross product of c and a
            // vector from atom(first number) to atom(second number)
            // second number can be -1 for com
            // c is always the cross product

            DPRINT("abc indicators: (%d %d) (%d %d)\n",
                    m_abc_indicators[0], m_abc_indicators[1], m_abc_indicators[2], m_abc_indicators[3]);

            for (size_t dim=0; dim<3; dim++)
            {
                a[dim] = positions_rel[3*m_abc_indicators[0] + dim];
                if (m_abc_indicators[1] != -1)
                    a[dim] -= positions_rel[3*m_abc_indicators[1] + dim];

                b[dim] = positions_rel[3*m_abc_indicators[2] + dim];
                if (m_abc_indicators[3] != -1)
                    b[dim] -= positions_rel[3*m_abc_indicators[3] + dim];
            }
            // normalize a
            cblas_sscal(3, 1/cblas_snrm2(3, &a[0], 1), &a[0], 1);

            // calculate c and b
            crossProduct(a, b, c);
            crossProduct(c, a, b);

            // normalize b and c
            cblas_sscal(3, 1/cblas_snrm2(3, &b[0], 1), &b[0], 1);
            cblas_sscal(3, 1/cblas_snrm2(3, &c[0], 1), &c[0], 1);

            DPRINT("abc:\n");
            DPRINT("%f %f %f\n", a[0], a[1], a[2]);
            DPRINT("%f %f %f\n", b[0], b[1], b[2]);
            DPRINT("%f %f %f\n", c[0], c[1], c[2]);

            // check if eigenvector points in same general direction as abc
            // if not flip eigenvector
            if (m_rot_treat != 'g')
            {
                if (cblas_sdot(3, &eigenvectors[0], 3, a, 1) < 0.0)
                    cblas_sscal(3, -1.0, &eigenvectors[0], 3);
                if (cblas_sdot(3, &eigenvectors[1], 3, b, 1) < 0.0)
                    cblas_sscal(3, -1.0, &eigenvectors[1], 3);
                if (cblas_sdot(3, &eigenvectors[2], 3, c, 1) < 0.0)
                    cblas_sscal(3, -1.0, &eigenvectors[2], 3);
            }

            DPRINT("eigenvectors unflipped:\n");
            DPRINT("%f %f %f\n", eigenvectors[0], eigenvectors[1], eigenvectors[2]);
            DPRINT("%f %f %f\n", eigenvectors[3], eigenvectors[4], eigenvectors[5]);
            DPRINT("%f %f %f\n", eigenvectors[6], eigenvectors[7], eigenvectors[8]);

            // calculate angular velocity in new coords (moi coords)
            float angular_velocity_nc[3] = {0.0, 0.0, 0.0};
            for (size_t dim=0; dim<3; dim++)
            {
                angular_velocity_nc[dim] = cblas_sdot(3, angular_momentum, 1,
                        &eigenvectors[dim], 3) / moments_of_inertia[dim];
            }
            DPRINT("angular velocity new coord: %8.4f%8.4f%8.4f\n",
                    angular_velocity_nc[0], angular_velocity_nc[1], angular_velocity_nc[2]);

            // calculate angular velocity and momentum in abc coords
            float abc[9];
            cblas_scopy(3, a, 1, &abc[0], 3);
            cblas_scopy(3, b, 1, &abc[1], 3);
            cblas_scopy(3, c, 1, &abc[2], 3);
            float angular_velocity_abc[3] = {0.0, 0.0, 0.0};
            float angular_momentum_abc[3] = {0.0, 0.0, 0.0};
            // angular_velocity_abc = angular_velocity @ abc
            for (size_t dim=0; dim<3; dim++)
            {
                angular_velocity_abc[dim] = cblas_sdot(3, angular_velocity, 1,
                        &abc[dim], 3);
                angular_momentum_abc[dim] = cblas_sdot(3, angular_momentum, 1,
                        &abc[dim], 3);
            }
            DPRINT("angular velocity abc coord: %8.4f%8.4f%8.4f\n",
                    angular_velocity_abc[0], angular_velocity_abc[1], angular_velocity_abc[2]);
            DPRINT("angular momentum abc coord: %8.4f%8.4f%8.4f\n",
                    angular_momentum_abc[0], angular_momentum_abc[1], angular_momentum_abc[2]);

            // calc moi_tensor in abc coords
            float temp_matrix[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
            float abc_inv[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
            float moi_tensor_abc[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
            inverseMatrix3x3(abc, abc_inv);
            // temp_matrix = abc_inv @ moi_tensor
            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0,
                    abc_inv, 3,
                    moi_tensor, 3,
                    0.0, temp_matrix, 3);
            DPRINT("temp_matrix:\n");
            DPRINT("%f %f %f\n", temp_matrix[0], temp_matrix[1], temp_matrix[2]);
            DPRINT("%f %f %f\n", temp_matrix[3], temp_matrix[4], temp_matrix[5]);
            DPRINT("%f %f %f\n", temp_matrix[6], temp_matrix[7], temp_matrix[8]);
            // moi_abc = temp_matrix @ abc
            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0,
                    temp_matrix, 3,
                    abc, 3,
                    0.0, moi_tensor_abc, 3);
            DPRINT("moi_tensor_abc:\n");
            DPRINT("%f %f %f\n", moi_tensor_abc[0], moi_tensor_abc[1], moi_tensor_abc[2]);
            DPRINT("%f %f %f\n", moi_tensor_abc[3], moi_tensor_abc[4], moi_tensor_abc[5]);
            DPRINT("%f %f %f\n", moi_tensor_abc[6], moi_tensor_abc[7], moi_tensor_abc[8]);

            // writing in output arrays the rotational velocities
            for (size_t dim=0; dim<3; dim++)
            {
                // 'f'ollowing the principal or abc axis; original method
                if (m_rot_treat == 'f')
                {
                    mol_omegas_sqrt_i_rot[3*nblocksteps*i + nblocksteps*dim + t] = angular_velocity_nc[dim] * sqrt(moments_of_inertia[dim]);
                }
                // 'a'bc as rotational axis
                else if (m_rot_treat == 'a')
                {
                    mol_omegas_sqrt_i_rot[3*nblocksteps*i + nblocksteps*dim + t] = angular_velocity_abc[dim] * sqrt(moi_tensor_abc[3*dim+dim]);
                }
                // a'b'c as rotational axis, but J/sqrt(I)
                else if (m_rot_treat == 'b')
                {
                    mol_omegas_sqrt_i_rot[3*nblocksteps*i + nblocksteps*dim + t] = angular_momentum_abc[dim] / sqrt(moi_tensor_abc[3*dim+dim]);
                }
                /*
                // angular velocity in 'x'yz (omega times sqrt(I)); does not give total rotational energy
                else if (m_rot_treat == 'x')
                {
                    mol_omegas_sqrt_i_rot[3*nblocksteps*i + nblocksteps*dim + t] = angular_velocity[dim] * sqrt(moi_tensor[3*dim+dim]);
                }
                // angular momentum in x'y'z (L over sqrt(I)); does not give total rotational energy
                else if (m_rot_treat == 'y')
                {
                    mol_omegas_sqrt_i_rot[3*nblocksteps*i + nblocksteps*dim + t] = angular_momentum[dim] / sqrt(moi_tensor[3*dim+dim]);
                }
                */

                // TODO catch other input
            }
            DPRINT("mol_omegas_sqrt_i_rot: %8.4f%8.4f%8.4f\n",
                    mol_omegas_sqrt_i_rot[3*nblocksteps*i + nblocksteps*0 + t],
                    mol_omegas_sqrt_i_rot[3*nblocksteps*i + nblocksteps*1 + t],
                    mol_omegas_sqrt_i_rot[3*nblocksteps*i + nblocksteps*2 + t]);
        }

        // print data for one frame
        DPRINT("molecules:\n");
        for (size_t i=0; i<nmols; i++)
        {
            DPRINT("%5zu%5zu%5zu%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f\n",
                    i, mol_natoms[i], mol_firstatom[i], mol_mass[i],
                    mol_velocities_sqrt_m_trn[3*nblocksteps*i + nblocksteps*0 + t], mol_velocities_sqrt_m_trn[3*nblocksteps*i + nblocksteps*1 + t],
                    mol_velocities_sqrt_m_trn[3*nblocksteps*i + nblocksteps*2 + t], mol_omegas_sqrt_i_rot[3*nblocksteps*i + nblocksteps*0 + t],
                    mol_omegas_sqrt_i_rot[3*nblocksteps*i + nblocksteps*1 + t], mol_omegas_sqrt_i_rot[3*nblocksteps*i + nblocksteps*2 + t]);
        }

        DPRINT("Step %lu (time %f) finished\n", t, time);
    }

    #pragma omp parallel
    {
        free(velocities);
        free(positions);
        free(positions_rel);
        free(velocities_rot);
    }

    // free help arrays
    free(r);
    free(v);
}
