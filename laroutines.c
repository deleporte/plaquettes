#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "laroutines.h"

int
dgghrd(char* COMPQ, char* COMPZ, int N, int ILO, int IHI,
       double* A, int LDA, double* B, int LDB, double* Q, int LDQ,
       double* Z, int LDZ)
{
  extern void dgghrd_(char* N, char* COMPZ, int *Np, int *ILOp,
		      int *IHIp,
		      double* A, int* LDAp, double* B, int* LDBp,
		      double* Q,  int* LDQp, double* Z, int *LDZp,
		      int* INFO);
  int INFO;
  dgghrd_(COMPQ,COMPZ,&N,&ILO,&IHI,A,&LDA,B,&LDB,Q,&LDQ,Z,&LDZ,&INFO);
  return INFO;
}

int
dhseqr(char* JOB, char* COMPZ, int N, int ILO, int IHI,
       double* H, int LDH, double* WR, double* WI, double* Z, int LDZ,
       double* WORK, int LWORK)
{
  extern void dhseqr_(char* JOB, char* COMPZ, int *Np, int *ILOp,
		      int *IHIp,
		      double* H, int *LDHp, double* WR, double* WI,
		      double* Z, int *LDZ,
		      double* WORK, int* LWORK, int* INFOp);
  int INFO;
  dhseqr_(JOB, COMPZ, &N, &ILO, &IHI, H, &LDH, WR, WI, Z, &LDZ, WORK, &LWORK,
	  &INFO);
  return INFO;
}

int
dhsein(char* SIDE, char* EIGSRC, char* INITV, int* SELECT,
       int N, double* H, int LDH, double* WR, double* WI,
       double* VL, int LDVL, double* VR, int LDVR, int MM, double* WORK,
       int* IFFAILL, int* IFFAILR)
{
  extern void dhsein_(char* SIDE, char* EIGSRC, char* INITV,
		      int* SELECT,
        int *Np, double* H, int *LDHp, double* WR, double* WI,
       double* VL, int *LDVLp, double* VR,
		       int *LDVRp, int *MMp, int *Mp,
		      double* WORK,
		      int* IFFAILL, int* IFFAILR, int* INFO);
  int INFO;
  int* M;
  M=malloc(sizeof(int));
  dhsein_(SIDE, EIGSRC, INITV, SELECT, &N, H, &LDH, WR, WI, VL, &LDVL,
	  VR, &LDVR, &MM, M, WORK, IFFAILL, IFFAILR,&INFO);
  free(M);
  return INFO;
}

int
dgeqrf(int M, int N, double* A, int LDA, double* TAU,
       double* WORK, int LWORK)
{
  extern void dgeqrf_(int* Mp, int* Np, double* A, int* LDAp,
		      double* TAU, double* WORK, int* LWORKp, int* INFOp);
  int INFO;
  dgeqrf_(&M,&N,A,&LDA,TAU,WORK,&LWORK,&INFO);
  return INFO;
}




int
dormqr(char* SIDE, char* TRANS, int M, int N, int K, double* A, int LDA,
       double* TAU, double* C, int LDC, double* WORK, int LWORK)
{
  extern void dormqr_(char* SIDE, char* TRANS, int* Mp, int* Np, int* Kp,
		      double* A, int* LDAp, double* TAU, double* C, int* LDCp,
		      double* WORK, int* LWORKp, int* INFOp);
  int INFO;
  dormqr_(SIDE,TRANS,&M,&N,&K,A,&LDA,TAU,C,&LDC,WORK,&LWORK,&INFO);
  return INFO;
}

int
deigvals(double* H, int dim, double* wr, double* wi, double prec){
  int info, converged, real;
  double* workspace;
  double* tau;
  double* Hmem;
  double* R;
  int workdim;
  double error,Delta;
  int i,j,k;

  tau=malloc(dim*dim*sizeof(double));
  R=malloc(dim*dim*sizeof(double));
  Hmem=malloc(dim*dim*sizeof(double));
  workspace=malloc(sizeof(double));
    //worspace dimension query
  dgeqrf(dim,dim,H,dim,tau,workspace,-1);
  workdim=workspace[0];
  free(workspace);
  workspace=malloc(workdim*workdim*sizeof(double));
  
  
  //We iteratively compute the QR factorisation until the subdiagonal entries
  //are small enough
  converged=0;
  k=0;
  while(converged==0){
    //decompose H_{k}=Q_{k}R_{k}
    for(i=0; i<dim*dim; i++){
      Hmem[i]=H[i];
    }
    dgeqrf(dim,dim,H,dim,tau,workspace,workdim);
    for(i=0; i<dim; i++){
      for(j=0; j<dim; j++){
	R[i*dim+j]=H[i*dim+j]; //row-major ?
      }
    }
    
    //Compute H_{k+1}=R_{k+1}Q_{k+1}
    dormqr("R","N",dim,dim,dim,R,dim,tau,H,dim,workspace,workdim);
    //check for convergence
    error=0;
    for(i=1; i<dim*dim; i++){
      error+=(H[i]-Hmem[i])*(H[i]-Hmem[i]);
    }
    if (error<prec*prec){
      converged=1;
    }
    k++;
  }
  printf("%d steps.\n",k);
  free(workspace);
  free(Hmem);
  free(tau);
  //now R holds the real Schur form
  //it must be close to diagonal except for complex pairs of eigenvalues
  k=0;
  for(i=0; i<dim; i++){
    real=1;
    if (isnan(R[i*dim+i])==0){
      for(j=i+1; j<dim; j++){
	if(abs(R[j*dim+i])>10*prec){  //col-major order !!
	  if(real==0){
	    wr[k]=NAN;
	    wi[k]=NAN;
	    wr[k+1]=NAN;
	    wi[k+1]=NAN;
	  }
	  else{
	    real=0;
	    Delta=pow(R[i*dim+i]+R[j*dim+j],2)
	      -4*(R[i*dim+i]*R[j*dim+j]-R[i*dim+j]*R[j*dim+i]);
	    if(Delta>=0){
	      wr[k]=NAN;
	      wi[k]=NAN;
	      wr[k+1]=NAN;
	      wi[k+1]=NAN;
	    }
	    else{
	      wr[k]=(R[i*dim+i]+R[j*dim+j])/2;
	      wr[k+1]=(R[i*dim+i]+R[j*dim+j])/2;
	      wi[k]=sqrt(-Delta)/2;
	      wi[k+1]=sqrt(-Delta)/2;
	    }
	  }
	  R[j*dim+j]=NAN;
	}
      }
      if(real=1){
	wr[k]=R[i*dim+i];
	wi[k]=0;
	k++;
      }
      if(real=0){
	k+=2;
      }
    }
  }
  if(k==dim){
    info=0;
  }
  else{
    info=1;
  }
  return info;
}
