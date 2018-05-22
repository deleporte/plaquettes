#include <iostream>
#include <armadillo>
#include <cmath>
#include "diag_arma.h"

using namespace std;
using namespace arma;

extern "C" void sparsediag(double* C, int dim, double* wr, double* wi,
		double* eigvl, double* eigvr)
{
  /* I don't trust this code
I expected armadillo to be able to compute both right and left
eigenvectors at the same time,this is not the case
This presumably leads to errors */
  //constructing the sparse matrix
  vector<double> data;
  vector<int> row_ind, col_ind;
  data.reserve(dim*2);
  row_ind.reserve(dim*2);
  col_ind.reserve(dim*2);
  double a,b;
  umat locations, tlocations;

  for(i=0; i<dim; i++){
    C[(2*i)%dim]=a;
    C[(2*i)%dim+1]=b;
    row_ind.push_back(i);
    row_ind.push_back(i);
    col_ind.push_back((2*i)%dim);
    col_ind.push_back((2*i)%dim+1);
    data.push_back(a/(a+b));
    data.push_back(b/(a+b));
  }
  locations.insert_rows(0,conv_to<urowvec>::from(row_ind));
  locations.insert_cols(1,conv_to<urowvec>::from(col_ind));
  tlocations.insert_rows(0,conv_to<urowvec>::from(col_ind));
  tlocations.insert_cols(1,conv_to<urowvec>::from(row_ind));
  
  sp_mat H=sp_mat(locations, conv_to<colvec>::from(data), dim,dim);
  sp_mat tH=sp_mat(tlocations, conv_to<colvec>::from(data), dim,dim);
  

  cx_vec eigval;
  cx_mat eigvecl;
  cx_mat eigvecr;

  //Actual diagonalisation
  eigs_gen(eigval,eigvecr, H, dim/2); //By default, large magnitude
  eigs_gen(eigval,eigvecl, tH, dim/2);
  //we expect half of the eigenvalues to be zero

  //convert eigval and eigvec
  for(i=0;i<dim/2; i++){
    wr[i]=real(eigval(i));
    wi[i]=imag(eigval(i));
    if(abs(wi[i]<0.0000001)){
      wi[i]=0;
      for(j=0; j<dim; j++){
	eigvl[i*dim+j]=eigvecl(i,j);
	eigvr[i*dim+j]=eigvecr(i,j);
      }
    }
    else{ //hopefully conjugate vectors are stored consecutively
      for(j=0; j<dim; j++){
	eigvl[i*dim+j]=real(eigvecl(i,j));
	eigvr[i*dim+j]=real(eigvecr(i,j));
      }
      for(j=0; j<dim; j++){
	eigvl[(i+1)*dim+j]=imag(eigvecl(i,j));
	eigvr[(i+1)*dim+j]=imag(eigvecr(i,j));
      }
      i++;
    }
  }
}

  
