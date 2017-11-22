#include "lssolve.h"

/** @brief a function for calculating the euclidean norm of a vector_t.
 *
 * @param in a vector_t
 * @return a double
 */
double norm(vector_t * in) {
	double norm = 0.0;
	for(unsigned long i = 0; i < in->n; i++) {
		norm += square(in->v[i]);
	}
	return sqrt(norm); 
}

/** @brief a support-function for the dgels routine of the LAPACK library.
 *
 *	A computation of the least-squares problem solved by the means of a
 *	QR factorization of the input matrix_t A. The solution is stored in the
 *	vector_t b. 
 *
 * @param A a pointer to a matrix_t
 * @param b a pointer to a vector_t
 * @return info an int
 */
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
	
	if(A->m < A->n) {
		printf("The matrix has more columns than rows, this would provide a solution for the LQ-factorization and not the desired QR-factorization. Terminating ..\n");
		return EXIT_FAILURE;
	}

	double b_norm = norm(b);
	
	int STATE = compute_dgels(A, b);

	if(STATE != 0) {
		if(STATE < 0) printf("The %d'th paramater provided to the dgels routine is illegal\n", abs(STATE));
		else printf("The provided matrix A does not have full rank, the %d'th diagonal element is zero.\n", STATE);
	}

	vector_t x;
	x.n = A->n;
	x.v = b->v;
	
	vector_t r;
	r.n = abs(A->n - A->m);
	r.v = b->v + A->n;
	
	double r_norm = norm(&r);

	printf("%lf\n", r_norm/b_norm);

	write_vector(argv[3], &x);
	free_vector(b);
	free_matrix(A);
	return STATE;
}
