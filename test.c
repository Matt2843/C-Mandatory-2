#include "matrix_io.h"

int main(int argc, char const *argv[]) {

  matrix_t *A = read_matrix("Al.txt");
  vector_t *b = read_vector("bl.txt");
  vector_t *x = read_vector("xl.txt");

  printf("A: %lu x %lu\n",A->m, A->n);
  printf("x: %lu\n",x->n);
  printf("b: %lu\n",b->n);

  if (b==NULL) printf("b is NULL\n");
  if (x==NULL) printf("x is NULL\n");
  if (A==NULL) printf("A is NULL\n");

  printf("freeing x\n");
  free_vector(x);
  printf("freeing b\n");
  free_vector(b);
  printf("freeing A\n");
  free_matrix(A);

  return 0;
}
