// -*- compile-command: "gcc main.c curvature.c dynamics.c python.c gnuplot_i.c $(python-config --cflags) $(python-config --libs) -lm -o myprog" -*-

#define SIZE 2
#define J -0.9

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include "laroutines.h"
#include "curvature.h"
#include "dynamics.h"
#include "gnuplot_i.h"
#include <Python.h>
#include "python.h"

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
  PyObject* main_module;
  PyObject* global_dict;
  FILE* pyfile=NULL;
  PyObject* Keys;

  Py_Initialize();
  
  main_module=PyImport_AddModule("__main__");
  pyfile=fopen("eigvals.py", "r");
  PyRun_SimpleFile(pyfile,"eigvals.py");
  global_dict = PyModule_GetDict(main_module);
  if(global_dict==NULL){
    printf("Could not load dictionary.\n");
  }
  Keys=PyDict_Keys(global_dict);
  for(i=0; i<PyList_Size(Keys); i++){
    printf("%s\n",PyString_AsString(PyList_GetItem(Keys,i)));
  }
  srand(time(NULL));
  hr=gnuplot_init();
  gnuplot_setstyle(hr, "points");
  gnuplot_set_xlabel(hr, "Coordinates");
  gnuplot_set_ylabel(hr, "Values");
  
  ham=malloc(16*sizeof(double));
  for(i=0; i<16; i++){
    ham[i]=0;
  }
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
  ham[0]=J;
  ham[5]=-J;
  ham[6]=2;
  ham[9]=2;
  ham[10]=-J;
  ham[15]=J;

  //Curvature testing
  /* C=generate_C(SIZE,0.1,1.1); */
  /* for(i=0; i<pow(2,SIZE); i++){ */
  /*   printf("C[%d]=%f\n",i,C[i]); */
  /* } */
  /* Gred=reducedCurvature(C,pow(2,SIZE),0.0000000000001); */
  /* for(i=0; i<pow(2,SIZE-1); i++){ */
  /*   for(j=0; j<pow(2,SIZE-1); j++){ */
  /*     printf("G[%d,%d]=%f\n",i,j,Gred[i*(int)pow(2,SIZE-1)+j]); */
  /*   } */
  /* } */
  /* free(C); */
  /* free(Gred); */
  
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
    oneStep(Cr,Ci,(int)pow(2,SIZE),ham,16,0.01/(1+abs(J*J))/pow(2,SIZE));
    printf("At step %d: \n",j+1);
    for(i=0; i<pow(2,SIZE); i++){
      printf("C[%d]=%f+%f i\n",i,Cr[i],Ci[i]);
    }
    //printf("%lf\n", meanval(Cr,Ci,(int)pow(2,SIZE),ham,16));
    //plot
    //gnuplot_resetplot(hr);
    //gnuplot_cmd(hr, "set yrange [-3:3]");
    //gnuplot_plot_x(hr,Cr,(int)pow(2,SIZE),"Real part");
    //gnuplot_plot_x(hr,Ci,(int)pow(2,SIZE),"Imag part");
    while(!getchar());
  }
  //printf("%lf\n", meanval(Cr,Ci,(int)pow(2,SIZE),ham,16));
  while(!getchar());
  free(Cr);
  free(Ci);
  free(ham);
  gnuplot_close(hr);
  Py_Finalize();
  return 0;
}
