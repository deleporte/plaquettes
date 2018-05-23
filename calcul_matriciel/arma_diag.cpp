/*
arma_diag.cpp

Contiendra des tests de diagonalisation avec armadillo

compilation : g++ arma_diag.cpp -o arma_diag -O2 -larmadillo -std=c++11

*/



#include <iostream>
#include <armadillo>
#include <iomanip>

using namespace std;
using namespace arma;

int main()
{
  mat A={ {2,1},
	  {1,1} };
  vec spec;
  mat U;
  eig_sym(spec,U,A);
  
  cout<<A*U.col(0)-spec(0)*U.col(0)<<endl;

  return 0;
}
