#include <assert.h>
#include "matrix_io.h"
#include "lssolve.h"

int main(int argc, char const *argv[]) {
	matrix_t *A = read_matrix("A.txt");
	vector_t *b = read_vector("b.txt");
	
	assert(A->n == 3);
	assert(A->m == 10);
	assert(b->n == A->n);
	
	
	double actual_norm = 0.0;
	double machine_epsilon = 0.0;
	assert(abs(norm(b) - actual_norm) <= machine_epsilon);
	
	int STATE = dgels_(A,b);
	assert(STATE == 0);
	
	vector_t r;
	double actual_residual = 0.0;
	assert(abs(norm(r)/norm(b) - actual_residual) <= machine_epsilon);

	free_vector(x);
	free_vector(b);
	free_matrix(A);
	return 0;
}
