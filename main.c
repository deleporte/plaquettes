// -*- compile-command: "gcc main.c curvature.c dynamics.c -lm -o myprog" -*-

#define SIZE 2

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "curvature.h"
#include "dynamics.h"

int main(int argc, char* argv[]){

  double* C = NULL;
  double* G = NULL;
  lkmat* T= NULL;
  lkmat* curr= NULL;
  int i,j;
  srand(time(NULL));
  C=generate_C(SIZE,0.1,1.1);
  G=reducedCurvature(C, (int)pow(2,SIZE), 0.00000000000001);
  
  free(G);
  T=dbMarkovPowers(C,(int)pow(2,SIZE), 50);

  //freeing T
  curr=T;
  while(curr->next != NULL){
    curr=curr->next;
  }
  while(curr->prev != NULL){
    curr=curr->prev;
    free(curr->next->data);
    free(curr->next);
  }
  free(curr->data);
  free(curr);
  free(C);
  return 0;
}
