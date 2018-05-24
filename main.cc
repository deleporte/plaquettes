// -*- compile-command: "g++ main.cc curvature.cc -o myprog.exe -O2 -larmadillo -std=c++11" -*-


// ce qu'on aimerait comme bouton
// choisir SIZE(>=2)
// retourner G ou Gred (attention elles ne sont pas de la meme taille)

#define SIZE 2
#define J -0.9


#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include "curvature.h"
//#include "dynamics.h"

using namespace std;
using namespace arma;


int main(int argc, char* argv[]){

  vec C;
  cx_vec w;
  //double* Cr = NULL;
  //double* Ci = NULL;
  //double* ham = NULL;
  mat Gred,T,G;
  cx_mat vel,ver;
  //double theta;
  //lkmat* T= NULL;
  //lkmat* curr= NULL;
  int i,j;

  srand(time(NULL));
  
  //ham=malloc(16*sizeof(double));
  //for(i=0; i<16; i++){
  //ham[i]=0;
  //}
  //Ising with transverse field
  /* ham[0]=J; */
  /* ham[5]=-J; */
  /* ham[10]=-J; */
  /* ham[15]=J; */
  /* ham[1]=(1-J)/2; */
  /* ham[4]=(1-J)/2; */
  /* ham[2]=(1-J)/2; */
  /* ham[8]=(1-J)/2; */
  /* ham[7]=(1-J)/2; */
  /* ham[13]=(1-J)/2; */
  /* ham[11]=(1-J)/2; */
  /* ham[14]=(1-J)/2; */
  
  //XXZ
  // ham[0]=J;
  //ham[5]=-J;
  //ham[6]=2;
  // ham[9]=2;
  //ham[10]=-J;
  //ham[15]=J;

  //Curvature testing
  generate_C(SIZE,0.1,1.1,C);
  for(i=0; i<pow(2,SIZE); i++){
    printf("C[%d]=%f\n",i,C[i]);
  }
  cout<<"C initialized."<<endl;
  diag(C,w,vel,ver,T);
  cout<<"Markov matrix diagonalized."<<endl;
  reducedCurvature(w,vel,ver,Gred);
  curvature(w,vel,ver,G);
  cout<<"G="<<G<<endl;
  cout<<"spectrum="<<eig_sym(G)<<endl;
  cout<<"Gred="<<Gred<<endl;
  cout<<"spectrum="<<eig_sym(Gred)<<endl;
  
  // C=generate_C(SIZE,0.1,1.1);
  // Cr=malloc((int)pow(2,SIZE)*sizeof(double));
  // Ci=malloc((int)pow(2,SIZE)*sizeof(double));
  // printf("Initial condition.\n");
  // for(i=0; i<pow(2,SIZE); i++){
  //   theta=2*M_PI*rand()/RAND_MAX;
  //   Cr[i]=C[i]*cos(theta);
  //   Ci[i]=C[i]*sin(theta);
  //   printf("C[%d]=%f+%f i\n",i,Cr[i],Ci[i]);
  // }
  // free(C);
  // printf("Press any key to continue.\n");
  // while(!getchar());
  // for(j=0; j<1000; j++){
  //   oneStep(Cr,Ci,(int)pow(2,SIZE),ham,16,0.01/(1+abs(J*J))/pow(2,SIZE));
  //   printf("At step %d: \n",j+1);
  //   for(i=0; i<pow(2,SIZE); i++){
  //     printf("C[%d]=%f+%f i\n",i,Cr[i],Ci[i]);
  //   }
  //   //printf("%lf\n", meanval(Cr,Ci,(int)pow(2,SIZE),ham,16));
  //   //plot
  //   //gnuplot_resetplot(hr);
  //   //gnuplot_cmd(hr, "set yrange [-3:3]");
  //   //gnuplot_plot_x(hr,Cr,(int)pow(2,SIZE),"Real part");
  //   //gnuplot_plot_x(hr,Ci,(int)pow(2,SIZE),"Imag part");
  //   while(!getchar());
  // }
  // //printf("%lf\n", meanval(Cr,Ci,(int)pow(2,SIZE),ham,16));
  // while(!getchar());
  // free(Cr);
  // free(Ci);
  // free(ham);
  // Py_Finalize();
  return 0;
}
