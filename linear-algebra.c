void crossProduct(float* x, float* y, float* c) 
{
    c[0] = x[1] * y[2] - x[2] * y[1];
    c[1] = x[2] * y[0] - x[0] * y[2];
    c[2] = x[0] * y[1] - x[1] * y[0];
    return;
}

void moiTensor(int mol_natoms, float* pos, float* atom_mass, float* tensor) 
{
    int a, b, j;
    for (a=0; a<3; a++)
    {
        for (b=0; b<3; b++)
        {
            for (j=0; j<mol_natoms; j++)
            {
                if (a == b)
                    tensor[3*a+b] += atom_mass[j] * ((pos[3*j+0]*pos[3*j+0] 
                                                      + pos[3*j+1]*pos[3*j+1] 
                                                      + pos[3*j+2]*pos[3*j+2]) 
                                                     - pos[3*j+a] * pos[3*j+b]);
                else
                    tensor[3*a+b] += atom_mass[j] * (- pos[3*j+a] * pos[3*j+b]);
            }
        }
    }
    return;
}

void inverseMatrix3x3(float* m, float* minv)
{
    float det = (m[0*3 + 0] * (m[1*3 + 1] * m[2*3 + 2] - m[2*3 + 1] * m[1*3 + 2])
                 - m[0*3 + 1] * (m[1*3 + 0] * m[2*3 + 2] - m[1*3 + 2] * m[2*3 + 0]) 
                 + m[0*3 + 2] * (m[1*3 + 0] * m[2*3 + 1] - m[1*3 + 1] * m[2*3 + 0]));
    float invdet = 1 / det;

    minv[0*3 + 0] = (m[1*3 + 1] * m[2*3 + 2] - m[2*3 + 1] * m[1*3 + 2]) * invdet;
    minv[0*3 + 1] = (m[0*3 + 2] * m[2*3 + 1] - m[0*3 + 1] * m[2*3 + 2]) * invdet;
    minv[0*3 + 2] = (m[0*3 + 1] * m[1*3 + 2] - m[0*3 + 2] * m[1*3 + 1]) * invdet;
    minv[1*3 + 0] = (m[1*3 + 2] * m[2*3 + 0] - m[1*3 + 0] * m[2*3 + 2]) * invdet;
    minv[1*3 + 1] = (m[0*3 + 0] * m[2*3 + 2] - m[0*3 + 2] * m[2*3 + 0]) * invdet;
    minv[1*3 + 2] = (m[1*3 + 0] * m[0*3 + 2] - m[0*3 + 0] * m[1*3 + 2]) * invdet;
    minv[2*3 + 0] = (m[1*3 + 0] * m[2*3 + 1] - m[2*3 + 0] * m[1*3 + 1]) * invdet;
    minv[2*3 + 1] = (m[2*3 + 0] * m[0*3 + 1] - m[0*3 + 0] * m[2*3 + 1]) * invdet;
    minv[2*3 + 2] = (m[0*3 + 0] * m[1*3 + 1] - m[1*3 + 0] * m[0*3 + 1]) * invdet;

    return;
}

void transposeMatrix3x3(float* m, float* mtrans)
{
    mtrans[0*3 + 0] = m[0*3 + 0];
    mtrans[0*3 + 1] = m[1*3 + 0];
    mtrans[0*3 + 2] = m[2*3 + 0];
    mtrans[1*3 + 0] = m[0*3 + 1];
    mtrans[1*3 + 1] = m[1*3 + 1];
    mtrans[1*3 + 2] = m[2*3 + 1];
    mtrans[2*3 + 0] = m[0*3 + 2];
    mtrans[2*3 + 1] = m[1*3 + 2];
    mtrans[2*3 + 2] = m[2*3 + 2];
    return;
}
