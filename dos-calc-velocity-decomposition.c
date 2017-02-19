#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <cblas.h>
#include <lapacke.h>
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/trrio.h"
#include "linear-algebra.c"

//#define DEBUG
#ifdef DEBUG
#define DPRINT(...) do{ fprintf( stderr, __VA_ARGS__ ); } while( 0 )
#else
#define DPRINT(...)
#endif

int decomposeVelocities (t_fileio* trj_in,
        gmx_trr_header_t header,
        long ntrajsteps,
        int natoms, 
        int nmols, 
        int nmoltypes,
        int* mol_firstatom,
        int* mol_natoms,
        int* mol_moltypenr,
        float* atom_mass,
        float* mol_mass,
        int* moltype_natomtypes,
        int* moltype_abc_indicators,
        float* mol_velocities_trn,
        float* omegas_sqrt_i,
        float* velocities_vib,
        float* mol_moments_of_inertia) 
{

    // allocate positions, velocities, ...
    float* positions = calloc(3*natoms, sizeof(float)); // this one could be smaller, but big molecules might arrive
    float* velocities = calloc(3*natoms, sizeof(float)); // this one could be smaller, but big molecules might arrive
    float* positions_rel = calloc(3*natoms, sizeof(float));
    float* velocities_rot = calloc(3*natoms, sizeof(float)); // this one could be smaller, but big molecules

    // for reading of frame
    long step;
    int t = 0;
    real time;
    real lambda;
    rvec box[3];
    rvec* x = malloc(header.x_size * ntrajsteps);
    rvec* v = malloc(header.v_size * ntrajsteps);

    DPRINT("start reading frame\n");

    while(gmx_trr_read_frame(trj_in, &step, &time, &lambda, box, &header.natoms, x, v, NULL))
    {
        if ( t >= ntrajsteps )
            break;

        DPRINT("There are %i atoms at step %i (time %f). My box is: %f %f %f \n",
                header.natoms, t, time, box[0][0], box[1][1], box[2][2]);

        // loop over molecules
        for (int i=0; i<nmols; i++)
        {
            DPRINT("\ndoing molecule %d\n", i);
            int m_firstatom = mol_firstatom[i];
            int m_natoms = mol_natoms[i];
            int m_moltype = mol_moltypenr[i];
            float m_mass = mol_mass[i];
            int* m_abc_indicators = &moltype_abc_indicators[4 * m_moltype];

            // read atoms of one molecule
            for (int j=0; j<m_natoms; j++)
            {
                int jj = m_firstatom + j;
                DPRINT("scanning atom %d (nr. %d in molecule)\n", jj, j);

                positions[3*j+0] = x[jj][0];
                positions[3*j+1] = x[jj][1];
                positions[3*j+2] = x[jj][2];
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
                for (int dim=0; dim<3; dim++)
                {
                    mol_velocities_trn[3*ntrajsteps*i + ntrajsteps*dim + t] = velocities[dim];
                    omegas_sqrt_i[3*ntrajsteps*i + ntrajsteps*dim + t] = 0;
                }

                for (int j=0; j<m_natoms; j++)
                {
                    int atom = m_firstatom + j;
                    for (int dim=0; dim<3; dim++)
                    {
                        velocities_vib[3*ntrajsteps*atom + ntrajsteps*dim + t] = 0;
                    }
                }

                continue;
            }


            // calc molecule velocity and molecule com
            float center_of_mass[3] = {0.0, 0.0, 0.0};
            float mol_velocity_trn[3] = {0.0, 0.0, 0.0};
            // this might also be done by a matrix multiplication but maybee overkill
            for (int dim=0; dim<3; dim++)
            {
                mol_velocity_trn[dim] = cblas_sdot(m_natoms,
                        &atom_mass[m_firstatom], 1,
                        &velocities[0+dim], 3);
                mol_velocity_trn[dim] /= m_mass;

                // output array
                DPRINT("going to output mol_vel_trn\n");
                DPRINT("i=%d t=%d\n", i, t);
                mol_velocities_trn[3*ntrajsteps*i + ntrajsteps*dim + t] = mol_velocity_trn[dim];
                DPRINT("finished output mol_vel_trn\n");

                center_of_mass[dim] = cblas_sdot(m_natoms,
                        &atom_mass[m_firstatom], 1,
                        &positions[0+dim], 3);
                center_of_mass[dim] /= m_mass;
            }

            DPRINT("mol_velocities_trn: %8.4f%8.4f%8.4f\n", 
                    mol_velocity_trn[0],
                    mol_velocity_trn[1],
                    mol_velocity_trn[2]);
            DPRINT("com: %8.4f%8.4f%8.4f\n", center_of_mass[0], center_of_mass[1], center_of_mass[2]);


            // calc molecule atoms relative positions
            for (int j=0; j<m_natoms; j++)
            {
                for (int dim=0; dim<3; dim++)
                {
                    positions_rel[3*j+dim] = positions[3*j+dim] - center_of_mass[dim];
                }

                DPRINT("pos_rel atom %d: %8.4f%8.4f%8.4f\n", j,
                        positions_rel[3*j+0], positions_rel[3*j+1], positions_rel[3*j+2]);

            }

            // calc angular momentum
            float cross_product[3] = {0.0, 0.0, 0.0};
            float angular_momentum[3] = {0.0, 0.0, 0.0};
            for (int j=0; j<m_natoms; j++)
            {
                crossProduct(&positions_rel[3*j], &velocities[3*j], cross_product);
                for (int dim=0; dim<3; dim++)
                {
                    angular_momentum[dim] += atom_mass[m_firstatom+j] * cross_product[dim];
                }
            }

            DPRINT("angular momentum: %8.4f%8.4f%8.4f\n", 
                    angular_momentum[0], angular_momentum[1], angular_momentum[2]);


            // calc moi tensor
            float moi_tensor[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
            moiTensor(m_natoms, positions_rel, &atom_mass[m_firstatom], moi_tensor);

            DPRINT("moi Tensor:\n");
            DPRINT("%f %f %f\n", moi_tensor[0], moi_tensor[1], moi_tensor[2]);
            DPRINT("%f %f %f\n", moi_tensor[3], moi_tensor[4], moi_tensor[5]);
            DPRINT("%f %f %f\n", moi_tensor[6], moi_tensor[7], moi_tensor[8]);


            // calc angular velocity
            float moi_tensor_inv[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
            float angular_velocity[3] = {0.0, 0.0, 0.0};
            inverseMatrix3x3(moi_tensor, moi_tensor_inv);
            cblas_sgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1.0,
                    moi_tensor_inv, 3,
                    angular_momentum, 1, 
                    0.0, angular_velocity, 1);

            DPRINT("angular velocity: %8.4f%8.4f%8.4f\n", 
                    angular_velocity[0], angular_velocity[1], angular_velocity[2]);


            // calc velocities rot
            for (int j=0; j<m_natoms; j++)
            {
                crossProduct(angular_velocity, &positions_rel[3*j], cross_product);
                for (int dim=0; dim<3; dim++)
                {
                    velocities_rot[3*j+dim] = cross_product[dim];
                }

                DPRINT("velocity_rot atom %d: %8.4f%8.4f%8.4f\n", j,
                        velocities_rot[3*j+0], velocities_rot[3*j+1], velocities_rot[3*j+2]);

            }

            // calc velocities vib
            float velocity_vib[3] = {0.0, 0.0, 0.0};
            for (int j=0; j<m_natoms; j++)
            {
                cblas_scopy(3, &velocities[3*j], 1, velocity_vib, 1);

                cblas_saxpy(3, -1.0, &velocities_rot[3*j], 1, 
                        velocity_vib, 1);
                cblas_saxpy(3, -1.0, mol_velocity_trn, 1, 
                        velocity_vib, 1);

                // output array
                int atom = m_firstatom + j;
                for (int dim=0; dim<3; dim++)
                {
                    velocities_vib[3*ntrajsteps*atom + ntrajsteps*dim + t] = velocity_vib[dim];
                }

                DPRINT("velocity_vib atom %d: %8.4f%8.4f%8.4f\n", j,
                        velocity_vib[0], 
                        velocity_vib[1], 
                        velocity_vib[2]);
            }

            // diatomic molecule
            // =================
            if (m_natoms == 2)
            {
                for (int dim=0; dim<3; dim++)
                {
                    float sqrt_moi = sqrt(moi_tensor[3*dim+dim]);
                    omegas_sqrt_i[3*ntrajsteps*i + ntrajsteps*dim + t] = angular_velocity[dim] * sqrt_moi;
                    mol_moments_of_inertia[3*i+dim] += sqrt_moi;
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
                printf("LAPACKE_ssyev failed to compute eigenvalues of\n");
                printf("the moment of inertia tensor of molecule %d\n", i);
                exit(1);
            }

            DPRINT("moments of inertia: %8.4f%8.4f%8.4f\n", 
                    moments_of_inertia[0], moments_of_inertia[1], moments_of_inertia[2]);

            // save moments of inertia for each molecule
            for (int dim=0; dim<3; dim++)
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
            // two numbers define a, two define b, c is cross product of a and b
            // vector from atom(first number) to atom(second number)
            // second number can be -1 for com 
            // c is always the cross product

            DPRINT("abc indicators: (%d %d) (%d %d)\n", 
                    m_abc_indicators[0], m_abc_indicators[1], m_abc_indicators[2], m_abc_indicators[3]);

            for (int dim=0; dim<3; dim++)
            {
                a[dim] = positions_rel[3*m_abc_indicators[0] + dim];
                if (m_abc_indicators[1] != -1)
                    a[dim] -= positions_rel[3*m_abc_indicators[1] + dim];

                b[dim] = positions_rel[3*m_abc_indicators[2] + dim];
                if (m_abc_indicators[3] != -1)
                    b[dim] -= positions_rel[3*m_abc_indicators[3] + dim];
            }

            crossProduct(a, b, c);

            DPRINT("abc:\n");
            DPRINT("%f %f %f\n", a[0], a[1], a[2]);
            DPRINT("%f %f %f\n", b[0], b[1], b[2]);
            DPRINT("%f %f %f\n", c[0], c[1], c[2]);

            // check if eigenvector points in same general direction as abc
            // if not flip eigenvector
            if (cblas_sdot(3, &eigenvectors[0], 3, a, 1) < 0.0)
                cblas_sscal(3, -1.0, &eigenvectors[0], 3);
            if (cblas_sdot(3, &eigenvectors[1], 3, b, 1) < 0.0)
                cblas_sscal(3, -1.0, &eigenvectors[1], 3);
            if (cblas_sdot(3, &eigenvectors[2], 3, c, 1) < 0.0)
                cblas_sscal(3, -1.0, &eigenvectors[2], 3);

            DPRINT("eigenvectors unflipped:\n");
            DPRINT("%f %f %f\n", eigenvectors[0], eigenvectors[1], eigenvectors[2]);
            DPRINT("%f %f %f\n", eigenvectors[3], eigenvectors[4], eigenvectors[5]);
            DPRINT("%f %f %f\n", eigenvectors[6], eigenvectors[7], eigenvectors[8]);

            // calculate angular velocity in new coords (moi coords)
            float angular_velocity_nc[3] = {0.0, 0.0, 0.0};
            for (int dim=0; dim<3; dim++)
            {
                angular_velocity_nc[dim] = cblas_sdot(3, angular_momentum, 1, 
                        &eigenvectors[dim], 3) / moments_of_inertia[dim];
            }
            DPRINT("angular velocity new coord: %8.4f%8.4f%8.4f\n", 
                    angular_velocity_nc[0], angular_velocity_nc[1], angular_velocity_nc[2]);

            // calculate angular velocity new coord times squareroot of moi
            for (int dim=0; dim<3; dim++)
            {
                omegas_sqrt_i[3*ntrajsteps*i + ntrajsteps*dim + t] = angular_velocity_nc[dim] * sqrt(moments_of_inertia[dim]);
            }
            DPRINT("omegas_sqrt_i: %8.4f%8.4f%8.4f\n", 
                    omegas_sqrt_i[3*ntrajsteps*i + ntrajsteps*0 + t], 
                    omegas_sqrt_i[3*ntrajsteps*i + ntrajsteps*1 + t], 
                    omegas_sqrt_i[3*ntrajsteps*i + ntrajsteps*2 + t]);
        }

        // print data for one frame
        DPRINT("molecules:\n");
        for (int i=0; i<nmols; i++)
        {
            DPRINT("%5d%5d%5d%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f\n",
                    i, mol_natoms[i], mol_firstatom[i], mol_mass[i],
                    mol_velocities_trn[3*ntrajsteps*i + ntrajsteps*0 + t], mol_velocities_trn[3*ntrajsteps*i + ntrajsteps*1 + t], 
                    mol_velocities_trn[3*ntrajsteps*i + ntrajsteps*2 + t], omegas_sqrt_i[3*ntrajsteps*i + ntrajsteps*0 + t],
                    omegas_sqrt_i[3*ntrajsteps*i + ntrajsteps*1 + t], omegas_sqrt_i[3*ntrajsteps*i + ntrajsteps*2 + t]);
        }

        DPRINT("Step %i (time %f) finished\n", t, time);
        t++;
    }    

    // divide by number of steps to get average molecule moment of inertia
    for (int i=0; i<nmols; i++)
    {
        for (int dim=0; dim<3; dim++)
        {
            mol_moments_of_inertia[3*i + dim] /= (float) ntrajsteps;
        }
        DPRINT("molecule %d moments of inertia: %8.4f%8.4f%8.4f\n", i,
                mol_moments_of_inertia[3*i+0], mol_moments_of_inertia[3*i+1], mol_moments_of_inertia[3*i+2]);
    }

    free(velocities);
    free(positions);
    free(positions_rel);
    free(velocities_rot);
    return 0;
}
