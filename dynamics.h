void force(cx_vec&, cx_vec&, cx_mat&, cx_mat&, cx_mat&, cx_vec&);
/*

Inputs
   - complex vector containing the coefficients
   - complex vector containing the eigenvalues
   - complex matrices with left and right eigenvectors
   - complex hamiltonian matrix
   - complex vector which will contain the force

Outputs:
   - sets the force

Code status: nothing
*/

void meanval(cx_vec&, cx_vec&, cx_mat&, cx_mat&, cx_mat&);
/*
Computes the mean value of an observable

Inputs
   - complex vector containing the coefficients
   - complex vector containing the eigenvalues
   - complex matrices with left and right eigenvectors
   - complex hamiltonian matrix

Outputs
   - returns a double containing the mean value
*/

void oneStep(cx_vec&, cx_mat&, double);
