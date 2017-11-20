#include "lssolve.h"

double norm(vector_t * in) {
	double norm = 0.0;
	for(unsigned long i = 0; i < in->n; i++) {
		norm += square(in->v[i]);
	}
	return sqrt(norm); 
}

int compute_dgels(matrix_t *A, vector_t *b) {
	const char trans = 'T';
	const int m = (int)(A->n);
	const int n = (int)(A->m);
	const int nrhs = 1;
	double * A_ = *(A->A);
	const int lda = m;
	double * B_ = b->v;
	const int ldb = n;
	int info = 0;
	int lwork_a = max(1, min(A->m, A->n) + max(min(A->m, A->n), b->n));
	vector_t * work_v = malloc_vector(lwork_a);
	double * work = work_v->v;
	int lwork = (int)(work_v->n);

	dgels_(&trans, &m, &n, &nrhs, A_, &lda, B_, &ldb, work, &lwork, &info);
	
	free_vector(work_v);
	return info;
}

int main(int argc, char * argv[]) {
	if(argc < 4) return EXIT_FAILURE;
	
	matrix_t * A = read_matrix(argv[1]);
	vector_t * b = read_vector(argv[2]);
	
	double b_norm = norm(b);
	
	int STATE = compute_dgels(A, b);
	
	vector_t x;
	x.n = 3;
	x.v = b->v;
	
	vector_t r;
	r.n = b->n - 3;
	r.v = b->v + 3;
	
	double r_norm = norm(&r);

	printf("%lf\n", r_norm/b_norm);

	write_vector(argv[3], &x);
	free_vector(b);
	free_matrix(A);
	return STATE;
}
