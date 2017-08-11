#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "curvature.h"
#include "dynamics.h"

#define max(a,b) \
  ({ __typeof__ (a) _a = (a); \
  __typeof__ (b) _b = (b); \
  _a > _b ? _a : _b; })

/* typedef struct lkmat lkmat; */
/* struct lkmat{ */
/*   double* data=NULL; */
/*   lkmat* next=NULL; */
/*   lkmat* prev=NULL; */
/* }; */

lkmat* dbMarkovPowers(double* C, int dim, int kmax){
  //As annoying as it might be we have to look at this array
  //It returns arrays of size len(C)*(len(C)**2)
  //T[i][a,b] is the proba, if we are in a at 0, to be in b at i
  //The data is huge so it is arranged as a linked list
  //It returns a pointer to i=0
  double** Ttemp=NULL;
  double* Cs=NULL;
  lkmat* T=NULL;
  lkmat* curr=NULL;
  double dtemp;
  int i,j,k;
  

  int d=(int)(log(dim+0.9)/log(2));
  Ttemp=malloc(d*sizeof(double*));
  
  Ttemp[0]=malloc(dim*sizeof(double));
  for(i=0; i<dim; i++){
    Ttemp[0][i]=C[i]/(C[(i/2)*2]+C[(i/2)*2+1]);
  }
    
  for(i=1; i<d; i++){
    Ttemp[i]=malloc(dim*(int)pow(2,i)*sizeof(double));
    for(j=0; j<dim*(int)pow(2,i); j++){
      Ttemp[i][j]=Ttemp[i-1][j/2]*Ttemp[0][j%dim];
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
  //check integrity
  for(i=0; i<dim; i++){
    dtemp=0;
    for(j=0; j<dim*dim; j++){
      dtemp += T->data[i*dim*dim+j];
    }
    if(fabs(dtemp-1)>0.01){
      printf("Error: T[0] is not stochastic\n");
      i=dim;
    }
  }
  
  curr=T;
  for(k=1; k<kmax; k++){
    curr->next=malloc(sizeof(lkmat));
    curr->next->prev=curr;
    curr->next->data=malloc(dim*dim*dim*sizeof(double));
    for(i=0; i<dim*dim*dim; i++){
      j=(i/dim/dim)*dim*dim+(i%(dim*dim))/2;
      curr->next->data[i]=Ttemp[0][i%(dim)]*(curr->data[j]+curr->data[j+dim*dim/2]);
    }
    curr=curr->next;
    //check integrity
    for(i=0; i<dim; i++){
      dtemp=0;
      for(j=0; j<dim*dim; j++){
	dtemp += curr->data[i*dim*dim+j];
      }
      if(fabs(dtemp-1)>0.01){
	printf("Error: T[%d] is not stochastic\n",k);
	i=dim;
      }
    }
  }

  curr->next=NULL;
  for(i=0; i<dim; i++){
    Ttemp[0][i]=C[i]/(C[i]+C[(i+dim/2)%dim]);
  }
  curr=T;
  for(k=0; k<kmax; k++){
    curr->prev=malloc(sizeof(lkmat));
    curr->prev->next=curr;
    curr->prev->data=malloc(dim*dim*dim*sizeof(double));
    for(i=0; i<dim*dim*dim; i++){
      j=(i/dim/dim)*dim*dim+(i*2)%(dim*dim);
      curr->prev->data[i]=Ttemp[0][(i/dim)%dim]*(curr->data[j]+curr->data[j+1]);
    }
    curr=curr->prev;
    //check integrity
    for(i=0; i<dim; i++){
      dtemp=0;
      for(j=0; j<dim*dim; j++){
	dtemp += curr->data[i*dim*dim+j];
      }
      if(fabs(dtemp-1)>0.01){
	printf("Error: T[%d] is not stochastic\n",-k);
	i=dim;
      }
    }
  }

  curr->prev=NULL;
  for(i=0; i<d; i++){
    free(Ttemp[i]);
  }
  free(Ttemp);
  return T;
}

double* force(double* Cr, double* Ci, int Cdim, double* ham, int hamdim){
  lkmat* T=NULL;
  lkmat* curr;
  double* gr=NULL;
  double* Cnorm=NULL;
  double* gi=NULL;
  double* Fr=NULL;
  double* Fi=NULL;
  int i,j,k,oj,js,ojs;
  int dimg=Cdim*Cdim*hamdim/4;
  int d=log(Cdim+0.9)/log(2);
  double tr,ti;

  Cnorm=malloc(Cdim*sizeof(double));
  for(i=0; i<Cdim; i++){
    Cnorm[i]=Cr[i]*Cr[i]+Ci[i]*Ci[i];
  }

  gr=malloc(dimg*sizeof(double));
  gi=malloc(dimg*sizeof(double));
  for(i=0; i<dimg; i++){
    gr[i]=ham[i/(dimg/(int)sqrt(hamdim))*(int)sqrt(hamdim)+(i/(Cdim/2))%(int)sqrt(hamdim)];
    gi[i]=0;
  }
  for(i=0; i<Cdim+sqrt(hamdim)-1; i++){
    for(j=0; j<dimg; j++){
      oj=j-((j/(Cdim/2))%(int)sqrt(hamdim))*Cdim/2;
      oj+=(j/(dimg/(int)sqrt(hamdim)))*Cdim/2;
      js=(j/(int)pow(2,i))%Cdim;
      ojs=(oj/(int)pow(2,i))%Cdim;
      tr=(Cr[ojs]*Cr[js]+Ci[ojs]*Ci[js])*gr[j];
      tr-=(Ci[ojs]*Cr[js]-Cr[ojs]*Ci[js])*gi[j];
      tr/=Cnorm[js];
      ti=(Cr[ojs]*Cr[js]+Ci[ojs]*Ci[js])*gi[j];
      ti+=(Ci[ojs]*Cr[js]-Cr[ojs]*Ci[js])*gr[j];
      ti/=Cnorm[js];
      gr[j]=tr;
      gi[j]=ti;
    }
  }
  T=dbMarkovPowers(Cnorm,Cdim,50);
  free(Cnorm);

  Fr=malloc(Cdim*2*sizeof(double));
  for(i=0; i<Cdim; i++){
    Fr[i]=0;
    Fr[i+Cdim]=0;
  }
  curr=T;
  while(curr->next != NULL){
    curr=curr->next;
  }
  while(curr->prev != NULL){
    for(i=0; i<Cdim; i++){
      for(j=0; j<dimg; j++){
	Fr[i]+=gr[j]*curr->data[i*Cdim*Cdim+j%(Cdim*Cdim)]; //hamdim=4
	Fr[i+Cdim]+=gi[j]*curr->data[i*Cdim*Cdim+j%(Cdim*Cdim)];
      }
    }
    curr=curr->prev;
  }
  //freeing lkmat
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
  free(gr);
  free(gi);
  return Fr;
}

void cholInv(double* Gred, int Gdim){
  double* L=NULL;
  double* Linv=NULL;
  int i,j,k,l;
  L=malloc(Gdim*Gdim*sizeof(double));
  //Cholesky decomposition
  for(i=0; i<Gdim; i++){
    if (Gred[i*(Gdim+1)]<=0){
      printf("Warning: Gred not positive at %d\n",i);
      printf("Value: %lf\n",Gred[i*(Gdim+1)]);
    }
    L[i*(Gdim+1)]=sqrt(Gred[i*(Gdim+1)]);
    for(j=1; j<Gdim-i; j++){
      L[i*(Gdim+1)+j]=Gred[i*(Gdim+1)+j]/sqrt(Gred[i*(Gdim+1)]);
    }
    for(k=1;k<Gdim-i; k++){
      for(l=1; l<Gdim-i; l++){
	Gred[(i+k)*Gdim+i+l]-=L[i*(Gdim+1)+k]*L[i*(Gdim+1)+l];
      }
    }
  }

  //inverse of the L matrix
  Linv=malloc(Gdim*Gdim*sizeof(double));
  for(i=0; i<Gdim; i++){
    Linv[i*(Gdim+1)]=1/L[i*(Gdim+1)];
    for(j=1; j<Gdim-i; j++){
      Linv[i*(Gdim+1)+j]=-L[i*(Gdim+1)+j]/L[i*(Gdim+1)]/L[j*(Gdim+1)];
    }
  }
  free(L);
  //inverse of the metric. From now on Gred contains its inverse
  for(i=0; i<Gdim; i++){
    for(j=0; j<Gdim; j++){
      Gred[i*(Gdim)+j]=0;
      for(k=max(i,j); k<Gdim; k++){
	Gred[i*(Gdim)+j]+=Linv[i*(Gdim)+k]*Linv[j*(Gdim)+k];
      }
    }
  }
  free(Linv);
}

void oneStep(double* Cr, double* Ci, int Cdim, double* ham, int hamdim, double step){
  //We must take large steps bc of computation time
  //We will use RK4 algorithm
  //double* Gred,Cnorm,k1,k2,k3,k4,kred,F,Ctr,Cti,U,V,UV,Fred;
  double* Gred=NULL;
  double* Cnorm=NULL;
  double* k1=NULL;
  double* k2=NULL;
  double* k3=NULL;
  double* k4=NULL;
  double* kred=NULL;
  double* F=NULL;
  double* Ctr=NULL;
  double* Cti=NULL;
  double* U=NULL;
  double* V=NULL;
  double* UV=NULL;
  double* Fred=NULL;
  int i,j,k;
  int d=log(0.9+Cdim)/log(2);
  double tval,tvalr,tvali;

  U=fullBasis(d);
  V=rotation(d);
  UV=malloc(Cdim*Cdim/2*sizeof(double));
  for(i=0; i<Cdim/2; i++){
    for(j=0; j<Cdim; j++){
      UV[i*Cdim+j]=0;
      for(k=0; k<Cdim; k++){
	UV[i*Cdim+j]+=U[(i+Cdim/2)*Cdim+k]*V[j*Cdim+k];
      }
    }
  }
  free(U);
  free(V);
  //I - Computation of k1
  //I.1 Compute the reduced curvature and its projected inverse
  Cnorm=malloc(Cdim*sizeof(double));
  for(i=0; i<Cdim; i++){
    Cnorm[i]=Cr[i]*Cr[i]+Ci[i]*Ci[i];
  }
  Gred=reducedCurvature(Cnorm,Cdim,0.00000000000001);
  
  
  free(Cnorm);
  cholInv(Gred,Cdim/2);
  for(i=0; i<Cdim*Cdim/4; i++){
    if(Gred[i]!=Gred[i]){
      printf("Nan error in curv for k1.\n");
      i=Cdim*Cdim/4;
    }
  }
  //I.2 Compute the force in the given basis
  F=force(Cr,Ci,Cdim,ham,hamdim);
  for(i=0; i<2*Cdim; i++){
    if(F[i]!=F[i]){
      printf("Nan error in force for k1.\n");
      i=2*Cdim;
    }
  }
  Fred=malloc(Cdim*sizeof(double));
  for(i=0; i<Cdim/2; i++){
    Fred[i]=0;
    for(j=0; j<Cdim; j++){
      Fred[i]+=UV[i*Cdim+j]*F[j];
    }

  }
  for(i=0; i<Cdim/2; i++){
    Fred[i+Cdim/2]=0;
    for(j=0; j<Cdim; j++){
      Fred[i+Cdim/2]+=UV[i*Cdim+j]*F[j+Cdim];
    }
  }
  for(i=0; i<Cdim; i++){
    if(Fred[i]!=Fred[i]){
      printf("Nan error in Fred for k1.\n");
      i=Cdim;
    }
  }
  
  free(F);
  //I.3 Compute the element of displacement
  kred=malloc(Cdim*sizeof(double));
  for(i=0; i<Cdim/2; i++){
    kred[i]=0;
    for(j=0; j<Cdim/2; j++){
      kred[i]+=Gred[i*Cdim/2+j]*Fred[j];
    }
  }
  for(i=0; i<Cdim/2; i++){
    kred[i+Cdim/2]=0;
    for(j=0; j<Cdim/2; j++){
      kred[i+Cdim/2]+=Gred[i*Cdim/2+j]*Fred[j+Cdim/2];
    }
  }
  if(kred[0]!=kred[0]){
    printf("Nan error in kred for k1.\n");
  }
  
  free(Fred);
  free(Gred);
  k1=malloc(2*Cdim*sizeof(double));
  for(i=0; i<Cdim; i++){
    k1[i]=0;
    for(j=0; j<Cdim/2; j++){
      k1[i]+=UV[j*Cdim+i]*kred[j];
    }
  }
  for(i=0; i<Cdim; i++){
    k1[i+Cdim]=0;
    for(j=0; j<Cdim/2; j++){
      k1[i+Cdim]+=UV[j*Cdim+i]*kred[j+Cdim/2];
    }
  }
  if(k1[0]!=k1[0]){
    printf("Nan error in k1.\n");
  }
  
  free(kred);
  //II - Computation of k2
  //II.1 : computation of the new position
  Ctr=malloc(Cdim*sizeof(double));
  Cti=malloc(Cdim*sizeof(double));
  for(i=0; i<Cdim; i++){
    Ctr[i]=Cr[i]*(1+step/2*k1[i+Cdim])+Ci[i]*step/2*k1[i];
    Cti[i]=Ci[i]*(1+step/2*k1[i+Cdim])-Cr[i]*step/2*k1[i];
  }
  if(Ctr[0]!=Ctr[0]){
    printf("Nan error in position for k2.\n");
  }
  
  //II.2 Compute the reduced curvature and its projected inverse
  Cnorm=malloc(Cdim*sizeof(double));
  for(i=0; i<Cdim; i++){
    Cnorm[i]=Ctr[i]*Ctr[i]+Cti[i]*Cti[i];
  }
  Gred=reducedCurvature(Cnorm,Cdim,0.00000000000001);
  if(Gred[0]!=Gred[0]){
    printf("Nan error in curv for k2.\n");
  }
  free(Cnorm);
  cholInv(Gred,Cdim/2);
  //II.3 Compute the force in the given basis
  F=force(Ctr,Cti,Cdim,ham,hamdim);
  if(F[0]!=F[0]){
    printf("Nan error in force for k2.\n");
  }
  free(Ctr);
  free(Cti);
  Fred=malloc(Cdim*sizeof(double));
  for(i=0; i<Cdim/2; i++){
    Fred[i]=0;
    for(j=0; j<Cdim; j++){
      Fred[i]+=UV[i*Cdim+j]*F[j];
    }
  }

  for(i=0; i<Cdim/2; i++){
    Fred[i+Cdim/2]=0;
    for(j=0; j<Cdim; j++){
      Fred[i+Cdim/2]+=UV[i*Cdim+j]*F[j+Cdim];
    }
  }
  free(F);
  //II.4 Compute the element of displacement
  kred=malloc(Cdim*sizeof(double));
  for(i=0; i<Cdim/2; i++){
    kred[i]=0;
    for(j=0; j<Cdim/2; j++){
      kred[i]+=Gred[i*Cdim/2+j]*Fred[j];
    }
  }
  for(i=0; i<Cdim/2; i++){
    kred[i+Cdim/2]=0;
    for(j=0; j<Cdim/2; j++){
      kred[i+Cdim/2]+=Gred[i*Cdim/2+j]*Fred[j+Cdim/2];
    }
  }  
  free(Fred);
  free(Gred);
  k2=malloc(2*Cdim*sizeof(double));
  for(i=0; i<Cdim; i++){
    k2[i]=0;
    for(j=0; j<Cdim/2; j++){
      k2[i]+=UV[j*Cdim+i]*kred[j];
    }
  }
  for(i=0; i<Cdim; i++){
    k2[i+Cdim]=0;
    for(j=0; j<Cdim/2; j++){
      k2[i+Cdim]+=UV[j*Cdim+i]*kred[j+Cdim/2];
    }
  }
  free(kred);

  //III - Computation of k3
  //III.1 : computation of the new position
  Ctr=malloc(Cdim*sizeof(double));
  Cti=malloc(Cdim*sizeof(double));
  for(i=0; i<Cdim; i++){
    Ctr[i]=Cr[i]*(1+step/2*k2[i+Cdim])+Ci[i]*step/2*k2[i];
    Cti[i]=Ci[i]*(1+step/2*k2[i+Cdim])-Cr[i]*step/2*k2[i];
  }
  //III.2 Compute the reduced curvature and its projected inverse
  Cnorm=malloc(Cdim*sizeof(double));
  for(i=0; i<Cdim; i++){
    Cnorm[i]=Ctr[i]*Ctr[i]+Cti[i]*Cti[i];
  }
  Gred=reducedCurvature(Cnorm,Cdim,0.00000000000001);
  if(Gred[0]!=Gred[0]){
    printf("Nan error in curv for k3.\n");
  }
  free(Cnorm);
  cholInv(Gred,Cdim/2);
  //III.3 Compute the force in the given basis
  F=force(Ctr,Cti,Cdim,ham,hamdim);
  if(F[0]!=F[0]){
    printf("Nan error in force for k3.\n");
  }
  free(Ctr);
  free(Cti);
  Fred=malloc(Cdim*sizeof(double));
  for(i=0; i<Cdim/2; i++){
    Fred[i]=0;
    for(j=0; j<Cdim; j++){
      Fred[i]+=UV[i*Cdim+j]*F[j];
    }
  }

  for(i=0; i<Cdim/2; i++){
    Fred[i+Cdim/2]=0;
    for(j=0; j<Cdim; j++){
      Fred[i+Cdim/2]+=UV[i*Cdim+j]*F[j+Cdim];
    }
  }
  free(F);
  //III.4 Compute the element of displacement
  kred=malloc(Cdim*sizeof(double));
  for(i=0; i<Cdim/2; i++){
    kred[i]=0;
    for(j=0; j<Cdim/2; j++){
      kred[i]+=Gred[i*Cdim/2+j]*Fred[j];
    }
  }
  for(i=0; i<Cdim/2; i++){
    kred[i+Cdim/2]=0;
    for(j=0; j<Cdim/2; j++){
      kred[i+Cdim/2]+=Gred[i*Cdim/2+j]*Fred[j+Cdim/2];
    }
  }  
  free(Fred);
  free(Gred);
  k3=malloc(2*Cdim*sizeof(double));
  for(i=0; i<Cdim; i++){
    k3[i]=0;
    for(j=0; j<Cdim/2; j++){
      k3[i]+=UV[j*Cdim+i]*kred[j];
    }
  }
  for(i=0; i<Cdim; i++){
    k3[i+Cdim]=0;
    for(j=0; j<Cdim/2; j++){
      k3[i+Cdim]+=UV[j*Cdim+i]*kred[j+Cdim/2];
    }
  }
  free(kred);

  //IIII - Computation of k4
  //IIII.1 : computation of the new position
  Ctr=malloc(Cdim*sizeof(double));
  Cti=malloc(Cdim*sizeof(double));
  for(i=0; i<Cdim; i++){
    Ctr[i]=Cr[i]*(1+step*k3[i+Cdim])+Ci[i]*step*k3[i];
    Cti[i]=Ci[i]*(1+step*k3[i+Cdim])-Cr[i]*step*k3[i];
  }
  //IIII.2 Compute the reduced curvature and its projected inverse
  Cnorm=malloc(Cdim*sizeof(double));
  for(i=0; i<Cdim; i++){
    Cnorm[i]=Ctr[i]*Ctr[i]+Cti[i]*Cti[i];
  }
  Gred=reducedCurvature(Cnorm,Cdim,0.00000000000001);
  if(Gred[0]!=Gred[0]){
    printf("Nan error in curv for k4.\n");
  }
  free(Cnorm);
  cholInv(Gred,Cdim/2);
  //IIII.3 Compute the force in the given basis
  F=force(Ctr,Cti,Cdim,ham,hamdim);
  if(F[0]!=F[0]){
    printf("Nan error in force for k4.\n");
  }
  free(Ctr);
  free(Cti);
  Fred=malloc(Cdim*sizeof(double));
  for(i=0; i<Cdim/2; i++){
    Fred[i]=0;
    for(j=0; j<Cdim; j++){
      Fred[i]+=UV[i*Cdim+j]*F[j];
    }
  }

  for(i=0; i<Cdim/2; i++){
    Fred[i+Cdim/2]=0;
    for(j=0; j<Cdim; j++){
      Fred[i+Cdim/2]+=UV[i*Cdim+j]*F[j+Cdim];
    }
  }
  free(F);
  //IIII.4 Compute the element of displacement
  kred=malloc(Cdim*sizeof(double));
  for(i=0; i<Cdim/2; i++){
    kred[i]=0;
    for(j=0; j<Cdim/2; j++){
      kred[i]+=Gred[i*Cdim/2+j]*Fred[j];
    }
  }
  for(i=0; i<Cdim/2; i++){
    kred[i+Cdim/2]=0;
    for(j=0; j<Cdim/2; j++){
      kred[i+Cdim/2]+=Gred[i*Cdim/2+j]*Fred[j+Cdim/2];
    }
  }  
  free(Fred);
  free(Gred);
  k4=malloc(2*Cdim*sizeof(double));
  for(i=0; i<Cdim; i++){
    k4[i]=0;
    for(j=0; j<Cdim/2; j++){
      k4[i]+=UV[j*Cdim+i]*kred[j];
    }
  }
  for(i=0; i<Cdim; i++){
    k4[i+Cdim]=0;
    for(j=0; j<Cdim/2; j++){
      k4[i+Cdim]+=UV[j*Cdim+i]*kred[j+Cdim/2];
    }
  }
  free(kred);

  //V - Yield the final value
  for(i=0; i<Cdim; i++){
    tvalr=Cr[i]*(1+step/6*(k1[i+Cdim]+2*k2[i+Cdim]+2*k3[i+Cdim]+k4[i+Cdim]))+Ci[i]*step/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
    tvali=Ci[i]*(1+step/6*(k1[i+Cdim]+2*k2[i+Cdim]+2*k3[i+Cdim]+k4[i+Cdim]))-Cr[i]*step/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
    
    Cr[i]=tvalr;
    Ci[i]=tvali;
  }
  //VI - Renormalize
  for(i=1; i<Cdim; i++){
    tvalr=(Cr[i]*Cr[0]+Ci[i]*Ci[0])/(Cr[0]*Cr[0]+Ci[0]*Ci[0]);
    tvali=(Ci[i]*Cr[0]-Cr[i]*Ci[0])/(Cr[0]*Cr[0]+Ci[0]*Ci[0]);
    Cr[i]=tvalr;
    Ci[i]=tvali;
  }
  Cr[0]=1;
  Ci[0]=0;
  free(UV);
  free(k1);
  free(k2);
  free(k3);
  free(k4);
}  
