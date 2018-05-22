double* generate_C(int,double,double);
/* Generates random values
Inputs:
   - int n concerning the size of the array to generate
   - two doubles containing the minimal and maximal possible values

Outputs:
   - a pointer to an array of size 2^n of random values
   - allocates memory for this pointer

Code status: proofed on 2017-08-14
*/

double* Markov(double*, int);

double* diag(double*, int, double*, double*, double*, double*);

double* curvature(double*, double*, double*, double*, int, double);

double* fullBasis(int);
/* Computes a change of basis with first vectors the null space of curvature
Input:
   - int concerning the shape of the matrix to generate

Output:
   - a pointer to an array of size 2^{2n} containing an orthogonal matrix in 
     colrow format
   - allocates memory for this pointer

Caveats:
   - the basis is expressed in terms of a rotated basis (hence the need for the
     function rotation()

Code status: proofed on 2017-08-14
*/

double* rotation(int);
/* Computes a rotation between the base along x and the base along z
Input:
   - int concerning the shape of the matrix to generate

Output:
   - a pointer to an array of size 2^{2n} concerning an orthogonal matrix in
     colrow format
   - allocates memory for this pointer

Code status: proofed on 2017-08-14
*/

double* reducedCurvature(double*, double*, double*, double*, int, double);
