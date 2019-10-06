#include <stdio.h>
#include <stdlib.h>

#ifndef DOSCALC_STRUCTS_DEF
#define DOSCALC_STRUCTS_DEF

typedef struct {
    size_t dofA_moltype;
    size_t dofB_moltype;
    char dofA_type;
    char dofB_type;
    size_t ndofA;
    size_t ndofB;
    size_t *dofA_list;
    size_t *dofB_list;
} dof_pair_def;

typedef struct {
    char name[80];
    char type;
    size_t ndof_pair_defs;
    dof_pair_def *dof_pair_defs;
} cross_spectrum_def;

#endif
