#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "curvature.h"
#include "dynamics.h"

#define max(a,b) \
  ({ __typeof__ (a) _a = (a); \
  __typeof__ (b) _b = (b); \
  _a > _b ? _a : _b; })

#define DESCENT 1


void force(cx_vec& C, cx_vec& w, cx_mat& vel, cx_mat& ver, cx_mat& ham, cx_vec& F){
//double* Cr, double* Ci, int dim, double* eigvr, double* eigvi, double* eigveleft, double* eigveright, double* ham, int hamdim){
  //double* F;
  uvec indices=sort_index(abs(w),"descend");
  int i,j,k,l;
  vec deq_m;
  vec eq_m=abs(ver.col(indices(0)));
  //double tvalr,tvali,ttvalr,ttvali;
  complex tval;
  vec Cnorm=abs(C);
  double norm;
  int sqrthamdim=ham.n_cols;
  int dim=C.size();
  F.zeros(dim);
  deqm.zeros(dim*dim*sqrthamdim/4);
  //for(i=0; i<dim*2; i++){
  //F[i]=0;
  //}
  //Cnorm=malloc(dim*sizeof(double));.
  //eq_m=malloc(dim*sizeof(double));
  //deq_m=malloc(dim*dim*sqrthamdim/4*sizeof(double));
  //for(i=0; i<dim; i++){
  //  Cnorm[i]=Cr[i]*Cr[i]+Ci[i]*Ci[i];
  //}

      //normalize
  norm=0;
  for(j=0; j<dim; j++){
    norm+=eq_m[j];
  }
  for(j=0; j<dim; j++){
    eq_m[j]/=norm;
  }
  //compute the double equilibrium measure
  for(i=0; i<dim*dim*sqrthamdim/4; i++){
    deq_m(i)=eq_m(i/(dim*sqrthamdim/4));
    for(k=1; k<dim*sqrthamdim/4; k*=2){
      deq_m(i)*=Cnorm((i/k)%dim)/(Cnorm((((i/k)%dim)/2)*2)+Cnorm((((i/k)%dim)/2)*2+1));
    }
  }

  //compute the force
  for(j=0; j<sqrthamdim; j++){
    for(i=0; i<dim*dim*sqrthamdim/4; i++){
      tval=deq_m[i];
      tval*=ham[j*sqrthamdim+(i/(dim/2))%sqrthamdim];
      //tvali=0;
      l=i%(dim/2)+j*dim/2+(i/(dim*sqrthamdim/2))*(dim*sqrthamdim/2);
      for(k=1; k<=dim*sqrthamdim/4; k*=2){
	tval=tval*C((l/k)%dim)/C((i/k)%dim);
	//ttvalr=tvalr*Cr[(l/k)%dim]-tvali*Ci[(l/k)%dim];
	//ttvali=tvalr*Ci[(l/k)%dim]+tvali*Cr[(l/k)%dim];
	//tvalr=ttvalr;
	//tvali=ttvali;
	//ttvalr=tvalr*Cr[(i/k)%dim]+tvali*Ci[(i/k)%dim];
	//ttvali=-tvalr*Ci[(i/k)%dim]+tvali*Cr[(i/k)%dim];
	//tvalr=ttvalr/(pow(Cr[(i/k)%dim],2)+pow(Ci[(i/k)%dim],2));
	//tvali=ttvali/(pow(Cr[(i/k)%dim],2)+pow(Ci[(i/k)%dim],2));
      }
      //multiply by the probability of this event knowing the other (summed over the distances)
      for(k=0; k<dim; k++){
	//if(eigvr[k]*eigvr[k]+eigvi[k]*eigvi[k]>0.000001 && eigvr[k]<0.99999){
	if(k!=indices(0)){
	  if(1){
	    for(l=0; l<dim; l++){
	      F(l)+=tval*w(k)/(1.-w(k))*vel(l,k)*ver(k,(i/(dim*sqrthamdim/4)));
	      F(l)+=tval*w(k)/(1-w(k))*vel((i%dim),k)*ver(k,l)*eq_m(i%dim)/eq_m(l);
	    }
	    //for(l=0; l<dim; l++){
	    //F[l+dim]+=tvali*eigvr[k]/(1-eigvr[k])*eigveleft[l*dim+k]
	    //	*eigveright[k*dim+(i/(dim*sqrthamdim/4))];
	    //F[l+dim]+=tvali*eigvr[k]/(1-eigvr[k])*eigveleft[(i%dim)*dim+k]
	    //	*eigveright[k*dim+l]
	    //	*eq_m[i%dim]/eq_m[l];
	    //}
	  }
	  // else{
	//     for(l=0; l<dim; l++){
	//       F[l]+=tvalr*2*((eigvr[k]*(1-eigvr[k])-eigvi[k]*eigvi[k])
	// 		     *(eigveleft[l*dim+k]*eigveright[k*dim+(i/(dim*sqrthamdim/4))]
	// 		       -eigveleft[l*dim+k+1]*eigveright[(k+1)*dim+(i/(dim*sqrthamdim/4))])
	// 		     -eigvi[k]*
	// 		     (eigveleft[l*dim+k]*eigveright[(k+1)*dim+(i/(dim*sqrthamdim/4))]
	// 		      +eigveleft[l*dim+k+1]*eigveright[k*dim+(i/(dim*sqrthamdim/4))]))
	// 	/(pow(1-eigvr[k],2)+pow(eigvi[k],2));
	//       F[l]+=tvalr*2*((eigvr[k]*(1-eigvr[k])-eigvi[k]*eigvi[k])
	// 		     *(eigveleft[(i%dim)*dim+k]*eigveright[k*dim+l]
	// 		       -eigveleft[(i%dim)*dim+k+1]*eigveright[(k+1)*dim+l])
	// 		     -eigvi[k]*
	// 		     (eigveleft[(i%dim)*dim+k]*eigveright[(k+1)*dim+l]
	// 		      +eigveleft[(i%dim)*dim+k+1]*eigveright[k*dim+l]))
	// 	/(pow(1-eigvr[k],2)+pow(eigvi[k],2))
	// 	*eq_m[i%dim]/eq_m[l];
	//     }
	//     for(l=0; l<dim; l++){
	//       F[l+dim]+=tvali*2*((eigvr[k]*(1-eigvr[k])-eigvi[k]*eigvi[k])
	// 		     *(eigveleft[l*dim+k]*eigveright[k*dim+(i/(dim*sqrthamdim/4))]
	// 		       -eigveleft[l*dim+k+1]*eigveright[(k+1)*dim+(i/(dim*sqrthamdim/4))])
	// 		     -eigvi[k]*
	// 		     (eigveleft[l*dim+k]*eigveright[(k+1)*dim+(i/(dim*sqrthamdim/4))]
	// 		      +eigveleft[l*dim+k+1]*eigveright[k*dim+(i/(dim*sqrthamdim/4))]))
	// 	/(pow(1-eigvr[k],2)+pow(eigvi[k],2));
	//       F[l+dim]+=tvali*2*((eigvr[k]*(1-eigvr[k])-eigvi[k]*eigvi[k])
	// 		     *(eigveleft[(i%dim)*dim+k]*eigveright[k*dim+l]
	// 		       -eigveleft[(i%dim)*dim+k+1]*eigveright[(k+1)*dim+l])
	// 		     -eigvi[k]*
	// 		     (eigveleft[(i%dim)*dim+k]*eigveright[(k+1)*dim+l]
	// 		      +eigveleft[(i%dim)*dim+k+1]*eigveright[k*dim+l]))
	// 	/(pow(1-eigvr[k],2)+pow(eigvi[k],2))
	// 	*eq_m[i%dim]/eq_m[l];
	//     }
	//     i++;
	//   }
	// }
	}
      }
    }
  }
    //return F;  
}
  
  

double meanval(cx_vec& C, cx_vec& w, cx_mat& vel, cx_mat& ver, cx_mat& ham){
//(double* Cr, double* Ci, int dim, double* eigvr, double* eigvi, double* eigveleft, double* eigveright, double* ham, int hamdim){
  uvec indices=sort_index(abs(w),"descend");
  double val=0;
  double norm=0;
  int i,j,k,l;
  double* deq_m;
  vec eq_m=abs(vel.col(indices(0)));
  vec Cnorm=abs(C);
  double tvalr,tvali,ttvalr,ttvali;
  int sqrthamdim=(int)sqrt(hamdim);
  Cnorm=malloc(dim*sizeof(double));
  eq_m=malloc(dim*sizeof(double));
  deq_m=malloc(dim*dim*sqrthamdim/4*sizeof(double));
  for(i=0; i<dim; i++){
    Cnorm[i]=Cr[i]*Cr[i]+Ci[i]*Ci[i];
  }
  for(i=0; i<dim; i++){
    if(fabs(eigvi[i])<0.00001 && eigvr[i]>0.99999){
      for(j=0; j<dim; j++){
	eq_m[j]=eigveright[i*dim+j];
      }
      //positive
      if(eq_m[0]<0){
	for(j=0; j<dim; j++){
	  eq_m[j]=-eq_m[j];
	}
      }
      //normalize
      norm=0;
      for(j=0; j<dim; j++){
	norm+=eq_m[j];
      }
      for(j=0; j<dim; j++){
	eq_m[j]/=norm;
      }
    }
  }
  for(i=0; i<dim*dim*sqrthamdim/4; i++){
    deq_m[i]=eq_m[i/(dim*sqrthamdim/4)];
    for(k=1; k<dim*sqrthamdim/4; k*=2){
      deq_m[i]*=Cnorm[(i/k)%dim]/(Cnorm[(((i/k)%dim)/2)*2]+Cnorm[(((i/k)%dim)/2)*2+1]);
    }
  }
  free(eq_m);

  for(j=0; j<sqrthamdim; j++){
    for(i=0; i<dim*dim*sqrthamdim/4; i++){
      tvalr=deq_m[i];
      tvalr*=ham[j*sqrthamdim+(i/(dim/2))%sqrthamdim];
      tvali=0;
      l=i%(dim/2)+j*dim/2+(i/(dim*sqrthamdim/2))*(dim*sqrthamdim/2);
      for(k=1; k<=dim*sqrthamdim/4; k*=2){
	ttvalr=tvalr*Cr[(l/k)%dim]-tvali*Ci[(l/k)%dim];
	ttvali=tvalr*Ci[(l/k)%dim]+tvali*Cr[(l/k)%dim];
	tvalr=ttvalr;
	tvali=ttvali;
	ttvalr=tvalr*Cr[(i/k)%dim]+tvali*Ci[(i/k)%dim];
	ttvali=-tvalr*Ci[(i/k)%dim]+tvali*Cr[(i/k)%dim];
	tvalr=ttvalr/(pow(Cr[(i/k)%dim],2)+pow(Ci[(i/k)%dim],2));
	tvali=ttvali/(pow(Cr[(i/k)%dim],2)+pow(Ci[(i/k)%dim],2));
      }
      val+=tvalr;
    }
  }
  free(deq_m);
  free(Cnorm);
  return val;
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
  double* eigvr=NULL;
  double* eigvi=NULL;
  double* eigveleft=NULL;
  double* eigveright=NULL;
  int i,j,k;
  int d=log(0.9+Cdim)/log(2);
  double tvalr,tvali;

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

  eigvr=malloc(Cdim*sizeof(double));
  eigvi=malloc(Cdim*sizeof(double));
  eigveleft=malloc(Cdim*Cdim*sizeof(double));
  eigveright=malloc(Cdim*Cdim*sizeof(double));
  //I - Computation of k1
  printf("k1\n");
  //I.1 Compute the reduced curvature and its projected inverse
  Cnorm=malloc(Cdim*sizeof(double));
  for(i=0; i<Cdim; i++){
    Cnorm[i]=Cr[i]*Cr[i]+Ci[i]*Ci[i];
  }
  diag(Cnorm,Cdim,eigvr,eigvi,eigveleft,eigveright);
  Gred=reducedCurvature(eigvr,eigvi,eigveleft,eigveright,Cdim,0.00000000000001);
  free(Cnorm);
  for(i=0; i<Cdim*Cdim/4; i++){
    if(Gred[i]!=Gred[i]){
      printf("Nan error in curv for k1.\n");
      i=Cdim*Cdim/4;
    }
    printf("%lf\n",Gred[i]);
  }
  cholInv(Gred,Cdim/2);
  for(i=0; i<Cdim*Cdim/4; i++){
    if(Gred[i]!=Gred[i]){
      printf("Nan error in curv for k1.\n");
      i=Cdim*Cdim/4;
    }
    printf("%lf\n",Gred[i]);
  }
  //I.2 Compute the force in the given basis
  F=force(Cr,Ci,Cdim,eigvr,eigvi,eigveleft,eigveright,ham,hamdim);
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
  printf("k2\n");
  //II.1 : computation of the new position
  Ctr=malloc(Cdim*sizeof(double));
  Cti=malloc(Cdim*sizeof(double));
  for(i=0; i<Cdim; i++){
    if(DESCENT){
      Ctr[i]=Cr[i]*(1-step/2*k1[i])+Ci[i]*step/2*k1[i+Cdim];
      Cti[i]=Ci[i]*(1-step/2*k1[i])-Cr[i]*step/2*k1[i+Cdim];
    }
    else{
      Ctr[i]=Cr[i]*(1+step/2*k1[i+Cdim])+Ci[i]*step/2*k1[i];
      Cti[i]=Ci[i]*(1+step/2*k1[i+Cdim])-Cr[i]*step/2*k1[i];
    }
  }
  if(Ctr[0]!=Ctr[0]){
    printf("Nan error in position for k2.\n");
  }
  
  //II.2 Compute the reduced curvature and its projected inverse
  Cnorm=malloc(Cdim*sizeof(double));
  for(i=0; i<Cdim; i++){
    Cnorm[i]=Ctr[i]*Ctr[i]+Cti[i]*Cti[i];
  }
  diag(Cnorm,Cdim,eigvr,eigvi,eigveleft,eigveright);
  Gred=reducedCurvature(eigvr,eigvi,eigveleft,eigveright,Cdim,0.00000000000001);
  free(Cnorm);
  for(i=0; i<Cdim*Cdim/4; i++){
    if(Gred[i]!=Gred[i]){
      printf("Nan error in curv for k2.\n");
      i=Cdim*Cdim/4;
    }
    printf("%lf\n",Gred[i]);
  }
  cholInv(Gred,Cdim/2);
  for(i=0; i<Cdim*Cdim/4; i++){
    if(Gred[i]!=Gred[i]){
      printf("Nan error in curv for k2.\n");
      i=Cdim*Cdim/4;
    }
    printf("%lf\n",Gred[i]);
  }
  //II.3 Compute the force in the given basis
  F=force(Ctr,Cti,Cdim,eigvr,eigvi,eigveleft,eigveright,ham,hamdim);
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
  printf("k3\n");
  //III.1 : computation of the new position
  Ctr=malloc(Cdim*sizeof(double));
  Cti=malloc(Cdim*sizeof(double));
  for(i=0; i<Cdim; i++){
    if(DESCENT){
      Ctr[i]=Cr[i]*(1-step/2*k2[i])+Ci[i]*step/2*k2[i+Cdim];
      Cti[i]=Ci[i]*(1-step/2*k2[i])-Cr[i]*step/2*k2[i+Cdim];
    }
    else{
      Ctr[i]=Cr[i]*(1+step/2*k2[i+Cdim])+Ci[i]*step/2*k2[i];
      Cti[i]=Ci[i]*(1+step/2*k2[i+Cdim])-Cr[i]*step/2*k2[i];
    }
  }
  //III.2 Compute the reduced curvature and its projected inverse
  Cnorm=malloc(Cdim*sizeof(double));
  for(i=0; i<Cdim; i++){
    Cnorm[i]=Ctr[i]*Ctr[i]+Cti[i]*Cti[i];
  }
  diag(Cnorm,Cdim,eigvr,eigvi,eigveleft,eigveright);
  Gred=reducedCurvature(eigvr,eigvi,eigveleft,eigveright,Cdim,0.00000000000001);
  free(Cnorm);
  for(i=0; i<Cdim*Cdim/4; i++){
    if(Gred[i]!=Gred[i]){
      printf("Nan error in curv for k3.\n");
      i=Cdim*Cdim/4;
    }
    printf("%lf\n",Gred[i]);
  }
  cholInv(Gred,Cdim/2);
  for(i=0; i<Cdim*Cdim/4; i++){
    if(Gred[i]!=Gred[i]){
      printf("Nan error in curv for k3.\n");
      i=Cdim*Cdim/4;
    }
    printf("%lf\n",Gred[i]);
  }
  //III.3 Compute the force in the given basis
  F=force(Ctr,Cti,Cdim,eigvr,eigvi,eigveleft,eigveright,ham,hamdim);
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
  printf("k4\n");
  //IIII.1 : computation of the new position
  Ctr=malloc(Cdim*sizeof(double));
  Cti=malloc(Cdim*sizeof(double));
  for(i=0; i<Cdim; i++){
    if(DESCENT){
      Ctr[i]=Cr[i]*(1-step*k3[i])+Ci[i]*step*k3[i+Cdim];
      Cti[i]=Ci[i]*(1-step*k3[i])-Cr[i]*step*k3[i+Cdim];
    }
    else{
      Ctr[i]=Cr[i]*(1+step*k3[i+Cdim])+Ci[i]*step*k3[i];
      Cti[i]=Ci[i]*(1+step*k3[i+Cdim])-Cr[i]*step*k3[i];
    }
  }
  if(Ctr[0]!=Ctr[0]){
    printf("Nan error in position for k4.\n");
  }
  //IIII.2 Compute the reduced curvature and its projected inverse
  Cnorm=malloc(Cdim*sizeof(double));
  for(i=0; i<Cdim; i++){
    Cnorm[i]=Ctr[i]*Ctr[i]+Cti[i]*Cti[i];
  }
  diag(Cnorm,Cdim,eigvr,eigvi,eigveleft,eigveright);
  Gred=reducedCurvature(eigvr,eigvi,eigveleft,eigveright,Cdim,0.00000000000001);

  free(Cnorm);
  for(i=0; i<Cdim*Cdim/4; i++){
    if(Gred[i]!=Gred[i]){
      printf("Nan error in curv for k4.\n");
      i=Cdim*Cdim/4;
    }
    printf("%lf\n",Gred[i]);
  }
  cholInv(Gred,Cdim/2);
  for(i=0; i<Cdim*Cdim/4; i++){
    if(Gred[i]!=Gred[i]){
      printf("Nan error in curv for k4.\n");
      i=Cdim*Cdim/4;
    }
    printf("%lf\n",Gred[i]);
  }
  //IIII.3 Compute the force in the given basis
  F=force(Ctr,Cti,Cdim,eigvr,eigvi,eigveleft,eigveright,ham,hamdim);
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
    if(DESCENT){
      tvalr=Cr[i]*(1-step/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]))+Ci[i]*step/6*(k1[i+Cdim]+2*k2[i+Cdim]+2*k3[i+Cdim]+k4[i+Cdim]);
      tvali=Ci[i]*(1-step/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]))-Cr[i]*step/6*(k1[i+Cdim]+2*k2[i+Cdim]+2*k3[i+Cdim]+k4[i+Cdim]);
    }
    else{
      tvalr=Cr[i]*(1+step/6*(k1[i+Cdim]+2*k2[i+Cdim]+2*k3[i+Cdim]+k4[i+Cdim]))+Ci[i]*step/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
      tvali=Ci[i]*(1+step/6*(k1[i+Cdim]+2*k2[i+Cdim]+2*k3[i+Cdim]+k4[i+Cdim]))-Cr[i]*step/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
    }
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
  free(eigvr);
  free(eigvi);
  free(eigveleft);
  free(eigveright);
}  
