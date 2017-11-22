#include <string.h>
#include <assert.h>
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
 *	vector_t b. The input matrix A always satisfies m >= n, but is not stored
 *	in the desired way for the dgels routine. The transposed solution does however,
 *	and therefore the trans parameter is thus set to 'T'.
 *	And by the documentation of the dgels routine, if TRANS = 'T' and m < n, rows 1 to M of B contain the
 *  least squares solution vectors; the residual sum of squares
 *  for the solution in each column is given by the sum of
 *  squares of elements M+1 to N in that column. Which is the desired solution and why
 *  we initially check that the condition is satisfied in the main function.
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

/** @brief a function containing tests for a fabricated example.
 * @return void
 */
int tests() {
	// Load in the predefined matrix and vector.
	matrix_t *A = read_matrix("A_test.txt"); 
	vector_t *b = read_vector("b_test.txt");
	
	assert(A->m == 10);
	assert(A->n == 3);
	assert(b->n == A->m);
	
	double actual_norm = 1.95336573636377664;
	const double double_machine_epsilon = 2.2204460492503131e-16;
	assert(abs(norm(b) - actual_norm) <= double_machine_epsilon);
	
	int STATE = compute_dgels(A, b);
	assert(STATE == 0);
	
	vector_t x;
	x.n = A->n;
	x.v = b->v;
	
	vector_t r;
	r.n = abs(A->n - A->m);
	r.v = b->v + A->n;
	
	// The solution for the least-squares problem calculated in maple.
	double ls1 = 1.05161973276760;
	double ls2 = -2.04158999101888;
	double ls3 = 1.04294431491269;
	
	assert(abs(x.v[0] - ls1) <= double_machine_epsilon);
	assert(abs(x.v[1] - ls2) <= double_machine_epsilon);
	assert(abs(x.v[2] - ls3) <= double_machine_epsilon);
	
	// The relative residual norm is calculated in matlab.
	double actual_relative_residual = 0.179810934931646;
	assert(abs(norm(&r)/actual_norm - actual_relative_residual) <= double_machine_epsilon);

	free_vector(b);
	free_matrix(A);
	return EXIT_SUCCESS;
}

/** @brief the main engine to test the implemented functions
 * 		
 * 	The function prints the residual norm
 * 	and saves the least-squares solution in the solution.txt file
 *
 * @param argv[1] a path to text-file containing the matrix A OR the string "test" to run the implemented tests.
 * @param argv[2] a path to text-file containing the vector b
 * @param argv[3] a path to text-file for writing solution
 * @return STATE the INFO param of dgels routine 
 */
int main(int argc, char * argv[]) {
	if(argc > 4) return EXIT_FAILURE;
	if(!strcmp(argv[1],"test")) {
		
		return tests();
	}
	
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
