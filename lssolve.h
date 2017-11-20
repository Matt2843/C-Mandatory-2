#ifndef LSSOLVE_H
#define LSSOLVE_H

#include "matrix_io.h"

/* C prototype for LAPACK routine DGELS */
void dgels_(
    const char * trans,  /* 'N' or 'T'             */
    const int * m,       /* rows in A              */
    const int * n,       /* cols in A              */
    const int * nrhs,    /* cols in B              */
    double * A,          /* array A                */
    const int * lda,     /* leading dimension of A */
    double * B,          /* array B                */
    const int * ldb,     /* leading dimension of B */
    double * work,       /* workspace array        */
    int * lwork,         /* workspace size         */
    int * info           /* status code            */
);

#endif
