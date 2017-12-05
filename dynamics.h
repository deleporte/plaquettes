typedef struct lkmat lkmat;
struct lkmat{
  double* data;
  lkmat* next;
  lkmat* prev;
};



lkmat* dbMarkovPowers(double*, int, int);

double meanval(double*, double*, int, double*, int);

double* force(double*, double*, int, double*, int);

void cholInv(double*, int);
/* Computes the inverse of a definite positive matrix via Cholesky decomposition

Inputs:
   - double pointer towards the symmetric matrix
   - int of the row dimension

Outputs:
   - No return
   - Changes the values in the pointer to store the inverted matrix
   - No data allocation

Caveats:
   - Symmetry of the input matrix is not checked. Only upper diagnoal
     values are used
   - Nonpositivity will trigger a warning to stdout but not stop execution
   - Only real-valued matrices are supported

Code status: proofed on 2017-08-14
 */

void oneStep(double*, double*, int, double*, int, double);
