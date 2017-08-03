#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "dynamics.h"

/* typedef struct lkmat lkmat; */
/* struct lkmat{ */
/*   double* data=NULL; */
/*   lkmat* next=NULL; */
/*   lkmat* prev=NULL; */
/* }; */

lkmat* dbMarkovPowers(double* C, int dim, int kmax){
  //As annoying as it might be we have to look at this array
  //It returns arrays of size len(C)*(len(C)**2)
  //ret[i][a,b] is the proba, if we are in a at 0, to be in b at i
  //The data is huge so it is arranged as a linked list
  //It returns a pointer to i=0
  double** Ttemp=NULL;
  double* Cs=NULL;
  lkmat* T=NULL;
  lkmat* curr=NULL;
  double* dtemp;
  int i,j,k,l,m,n;
  int d=(int)(log(dim+0.9)/log(2));
  Ttemp=malloc(d*sizeof(double*));
  
  Ttemp[0]=malloc(dim*sizeof(double));
  for(i=0; i<dim; i++){
    Ttemp[0][i]=C[i]/(C[(i/2)*2]+C[(i/2)*2+1]);
  }
  for(i=1; i<d; i++){
    Ttemp[i]=malloc(dim*(int)pow(2,i)*sizeof(double));
    for(j=0; j<dim*(int)pow(2,i); j++){
      Ttemp[i][j]=Ttemp[i-1][j%2]*Ttemp[0][j%dim];
    }
  }
  T=malloc(sizeof(lkmat));
  T->data=malloc(dim*dim*dim*sizeof(double));
  for(i=0; i<dim*dim*dim; i++){
    if(i/dim/dim==(i/dim)%dim){
      T->data[i]=Ttemp[d-1][i%(dim*dim/2)];
    }
    else{
      T->data[i]=0;
    }
  }
  curr=T;
  for(k=1; k<kmax; k++){
    curr->next=malloc(sizeof(lkmat));
    curr->next->prev=curr;
    curr->next->data=malloc(dim*dim*dim*sizeof(double));
    dtemp=curr->data+dim*dim/2;
    for(i=0; i<dim*dim*dim; i++){
      j=(i/dim/dim)*dim*dim+(i%(dim*dim))/2;
      curr->next->data[i]=Ttemp[0][i%(dim)]*(curr->data[j]+dtemp[j]);
    }
    curr=curr->next;
  }
  curr->next=NULL;
  dtemp=C+dim/2;
  for(i=0; i<dim; i++){
    Ttemp[0][i]=C[i]/(C[i/2]+dtemp[i/2]);
  }
  curr=T;
  for(k=0; k<kmax; k++){
    curr->prev=malloc(sizeof(lkmat));
    curr->prev->next=curr;
    curr->prev->data=malloc(dim*dim*dim*sizeof(double));
    for(i=0; i<dim*dim*dim; i++){
      j=(i/dim/dim)*dim*dim+(i*2)%(dim*dim);
      curr->next->data[i]=Ttemp[0][i/dim/dim]*(curr->data[j]+curr->data[j+1]);
    }
    curr=curr->prev;
  }
  curr->prev=NULL;
  for(i=0; i<d; i++){
    free(Ttemp[i]);
  }
  free(Ttemp);
  return T;
}
