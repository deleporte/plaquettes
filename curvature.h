#include <armadillo>

using namespace arma;

void generate_C(int,double,double,vec&);
/* Generates random values
Inputs:
   - int n concerning the size of the array to generate
   - two doubles m and M containing the minimal and maximal possible values
   - pointer vC on colvec to be filled

Outputs:
   - resizes vC to 2^n and sets random values between m and M

Code status: needs check
*/

void Markov(vec&, mat&);
/* Generates Markov chain associated with plaquette coefficients
Inputs:
   - pointer C on colvec containing the coefficients
   - pointer T on mat to be filled

Outputs:
   - resizes T to shape (size(C),size(C))
   - fills the non-empty values of T

Code status: needs check
*/

void diag(vec&, cx_vec&, cx_mat&, cx_mat&, mat&);
/* Computes and diagonalizes the Markov chain matrix
Inputs:
   - pointer C on colvec containing the coefficients
   - pointer w on cx_vec which will store the eigenvalues
   - pointers vel,ver,T which will respectively contain left and right eigvecs and Markov matrix

Outputs:
   - changes w, vel, ver and T to the corresponding values

Code status: needs check
*/

void curvature(cx_vec&, cx_mat&, cx_mat&, mat&);
/* Computes curvature in the standard basis, from the eigendecomposition
Inputs:
   - eigenvalues, left and right eigenvectors, empty matrix
Outputs:
   - changes G

Code status: needs check
*/

void fullBasis(int,mat&);
/* Computes a change of basis with first vectors the null space of curvature
Input:
   - int concerning the shape of the matrix to generate
   - empty matrix
Output:
   - Changes the matrix

Caveats:
   - the basis is expressed in terms of a rotated basis (hence the need for the
     function rotation()

Code status: needs check
*/

void rotation(int,mat&);
/* Computes a rotation between the base along x and the base along z
Input:
   - int concerning the shape of the matrix to generate
   - empty matrix V

Output:
   - Changes the matrix

Code status: needs check
*/

void reducedCurvature(cx_vec&, cx_mat&, cx_mat&, mat&);
/* Computes curvature in a convenient basis, from the eigendecomposition
Inputs:
   - eigenvalues, left and right eigenvectors, empty matrix
Outputs:
   - changes Gred

Code status: runs but can give negative matrices
*/
