#ifndef LSSOLVE_H
#define LSSOLVE_H
#include <math.h>
#include "matrix_io.h"

#define max(x,y) (((x) > (y)) ? (x) : (y))
#define min(x,y) (((x) < (y)) ? (x) : (y))
#define square(x) ((x) * (x))

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

// Function prototypes
int tests();
double norm(vector_t*);
int compute_dgels(matrix_t*, vector_t*);

#endif
