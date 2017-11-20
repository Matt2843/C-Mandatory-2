#include "matrix_io.h"
#include "lssolve.h"
#define max(x,y) (((x) > (y)) ? (x) : (y))
#define min(x,y) (((x) < (y)) ? (x) : (y))

int compute_dgels(matrix_t *A, vector_t *b) {
	const char trans = 'T';
	const int m = (int)(A->n);
	const int n = (int)(A->m);
	const int nrhs = 1;
	double * A_lel = *(A->A);
	const int lda = m;
	double * B_lel = b->v;
	const int ldb = n;
	int info = 0;
	int lwork_a = max(1, min(A->m, A->n) + max(min(A->m, A->n), b->n));
	vector_t * work_v = malloc_vector(lwork_a);
	double * work = work_v->v;
	int lwork = (int)(work_v->n);

	dgels_(&trans, &m, &n, &nrhs, A_lel, &lda, B_lel, &ldb, work, &lwork, &info);
	
	free_vector(work_v);
	return info;
}

int main(int argc, char * argv[]) {
	if(argc < 4) return EXIT_FAILURE;
	
	matrix_t * A = read_matrix(argv[1]);
	vector_t * b = read_vector(argv[2]);
	
	int STATE = compute_dgels(A, b);
	
	print_matrix(A);
	
	vector_t res;
	res.n = 3;
	res.v = b->v;

	print_vector(&res);
	write_vector(argv[3], &res);

	free_vector(b);
	free_matrix(A);
	return STATE;
}
