#include <cblas.h>
#include <lapacke.h>
#include <math.h>

void crossProduct(float *x, float *y, float *c) {
    c[0] = x[1] * y[2] - x[2] * y[1];
    c[1] = x[2] * y[0] - x[0] * y[2];
    c[2] = x[0] * y[1] - x[1] * y[0];
    return;
}

void crossProductInc(float *x, size_t incx, float *y, size_t incy, float *c,
                     size_t incc) {
    c[0 * incc] = x[1 * incx] * y[2 * incy] - x[2 * incx] * y[1 * incy];
    c[1 * incc] = x[2 * incx] * y[0 * incy] - x[0 * incx] * y[2 * incy];
    c[2 * incc] = x[0 * incx] * y[1 * incy] - x[1 * incx] * y[0 * incy];
    return;
}

void moiTensorMixed(int mol_natoms, float *pos1, float *pos2, float *atom_mass,
                    float *tensor) {
    int a, b, j;
    for (a = 0; a < 9; a++)
        tensor[a] = 0.0;
    for (a = 0; a < 3; a++) {
        for (b = 0; b < 3; b++) {
            for (j = 0; j < mol_natoms; j++) {
                if (a == b) {
                    tensor[3 * a + b] +=
                        atom_mass[j] * (pos1[3 * j + 0] * pos2[3 * j + 0] +
                                        pos1[3 * j + 1] * pos2[3 * j + 1] +
                                        pos1[3 * j + 2] * pos2[3 * j + 2]);
                }
                tensor[3 * a + b] -=
                    atom_mass[j] * (pos1[3 * j + a] * pos2[3 * j + b]);
            }
        }
    }
}

void moiTensor(int mol_natoms, float *pos, float *atom_mass, float *tensor) {
    moiTensorMixed(mol_natoms, pos, pos, atom_mass, tensor);
}

// invert square matrix
lapack_int invert_matrix(float *A, unsigned n) {
    int ipiv[n + 1];
    lapack_int ret;
    ret = LAPACKE_sgetrf(LAPACK_ROW_MAJOR, n, n, A, n, ipiv);
    if (ret != 0)
        return ret;
    ret = LAPACKE_sgetri(LAPACK_ROW_MAJOR, n, A, n, ipiv);
    return ret;
}

// squareroot of matrix
lapack_int squareroot_of_matrix(float *A, unsigned N) {
    lapack_int ret;
    // eigenvectors
    float eigenvalues[N];
    ret = LAPACKE_ssyev(LAPACK_ROW_MAJOR, 'V', 'U', N, A, N, eigenvalues);
    if (ret != 0)
        return ret;
    // calculate B = eigenvectors @ diag(sqrt(eigenvalues))
    float B[N * N];
    for (unsigned a = 0; a < N; a++) {
        for (unsigned b = 0; b < N; b++) {
            B[a * N + b] = A[a * N + b] * sqrt(eigenvalues[b]);
        }
    }
    // calculate C = B @ eigv(A).T
    float C[N * N];
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, N, N, N, 1.0, B, N, A,
                N, 0.0, C, N);
    cblas_scopy(N * N, C, 1, A, 1);
    return ret;
}
