// -*- compile-command: "gcc main.c curvature.c dynamics.c gnuplot_i.c -lm -o myprog" -*-

#define SIZE 2
#define J 0.5

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include "curvature.h"
#include "dynamics.h"
#include "gnuplot_i.h"

int main(int argc, char* argv[]){

  double* C = NULL;
  double* Cr = NULL;
  double* Ci = NULL;
  double* G = NULL;
  double* ham = NULL;
  double* Gred = NULL;
  double theta;
  lkmat* T= NULL;
  lkmat* curr= NULL;
  int i,j;
  gnuplot_ctrl* hr;
  
  srand(time(NULL));
  hr=gnuplot_init();
  gnuplot_setstyle(hr, "points");
  gnuplot_set_xlabel(hr, "Coordinates");
  gnuplot_set_ylabel(hr, "Values");
  
  ham=malloc(16*sizeof(double));
  ham[0]=J;
  ham[5]=-J;
  ham[10]=-J;
  ham[15]=J;
  ham[1]=(1-J)/2;
  ham[4]=(1-J)/2;
  ham[2]=(1-J)/2;
  ham[8]=(1-J)/2;
  ham[7]=(1-J)/2;
  ham[13]=(1-J)/2;
  ham[11]=(1-J)/2;
  ham[14]=(1-J)/2;
  
    
  C=generate_C(SIZE,0.1,1.1);
  Cr=malloc((int)pow(2,SIZE)*sizeof(double));
  Ci=malloc((int)pow(2,SIZE)*sizeof(double));
  printf("Initial condition.\n");
  for(i=0; i<pow(2,SIZE); i++){
    theta=2*M_PI*rand()/RAND_MAX;
    Cr[i]=C[i]*cos(theta);
    Ci[i]=C[i]*sin(theta);
    printf("C[%d]=%f+%f i\n",i,Cr[i],Ci[i]);
  }
  free(C);
  printf("Press any key to continue.\n");
  while(!getchar());
  for(j=0; j<1000; j++){
    oneStep(Cr,Ci,(int)pow(2,SIZE),ham,16,0.01);
    printf("At step %d: \n",j+1);
    for(i=0; i<pow(2,SIZE); i++){
      printf("C[%d]=%f+%f i\n",i,Cr[i],Ci[i]);
    }
    //plot
    gnuplot_resetplot(hr);
    gnuplot_cmd(hr, "set yrange [-3:3]");
    gnuplot_plot_x(hr,Cr,(int)pow(2,SIZE),"Real part");
    gnuplot_plot_x(hr,Ci,(int)pow(2,SIZE),"Imag part");
  }
  while(!getchar());
  free(Cr);
  free(Ci);
  free(ham);
  gnuplot_close(hr);
  return 0;
}
