#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
//#include <atlas_enum.h>
#include "curvature.h"
//#include "clapack.h"
//lapack functions definitions
#include "python.h"

double* generate_C(int d, double m, double M){
  double* C=NULL;
  int i;

  C=malloc((int)pow(2,d)*sizeof(double));
  for(i=0; i<pow(2,d); i++){
    C[i]=((double)rand())/RAND_MAX*(M-m)+m;
  }
  return C;
}

double* Markov(double* C, int dim){
  double* T;
  int i;
  T=malloc(dim*dim*sizeof(double));
  for(i=0; i<dim; i++){
    T[((2*i)%dim)*dim+i]=C[(2*i)%dim]/(C[(2*i)%dim]+C[(2*i)%dim+1]);
    T[((2*i)%dim+1)*dim+i]=C[(2*i)%dim+1]/(C[(2*i)%dim]+C[(2*i)%dim+1]);
  }
  return T;
}

double* diag(double* C, int dim, double* wr, double* wi, double* vel, double* ver){
  double* T=NULL;
  T=Markov(C,dim);
  eigvals(T,dim,wr,wi,vel,ver);
  return T;
}

double* curvature(double* eigvr, double* eigvi, double* eigveleft, double* eigveright, int dim, double error){
  double* eq_m=NULL;
  double* G=NULL;
  int i,j,k,l;
  double norm;
  eq_m=malloc(dim*sizeof(double));
  if(1){
    if(1){
      //find the equilibrium state
      for(i=0; i<dim; i++){
	if(fabs(eigvi[i])<0.00001 && eigvr[i]>0.9999){
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
    }
  }
  G=malloc(dim*dim*sizeof(double));
  for(i=0; i<dim*dim; i++){
    if(i%(dim+1)){
      G[i]=-eq_m[i%dim]*eq_m[i/dim];
    }
    else{
      G[i]=eq_m[i%dim]*(1-eq_m[i/dim]);
    }
  }
  if(1){
    for(i=0; i<dim; i++){
      if(eigvr[i]*eigvr[i]+eigvi[i]*eigvi[i]>0.000001 && eigvr[i]<0.9999){
	if(fabs(eigvi[i])<0.00001){
	  for(j=0; j<dim; j++){
	    for(l=0; l<dim; l++){
	      G[j*dim+l]+=eigvr[i]/(1-eigvr[i])*eq_m[l]
		*(eigveleft[l*dim+i]*eigveright[i*dim+j]);
	      G[l*dim+j]+=eigvr[i]/(1-eigvr[i])*eq_m[l]
		*(eigveleft[l*dim+i]*eigveright[i*dim+j]);
	    }
	  }
	}
	else{
	  for(j=0; j<dim; j++){
	    for(l=0; l<dim; l++){
	      G[j*dim+l]+=2*((eigvr[i]*(1-eigvr[i])-eigvi[i]*eigvi[i])
		*(eigveleft[j*dim+i]*eigveright[i*dim+l]
		  -eigveleft[j*dim+i+1]*eigveright[(i+1)*dim+l])
		-(eigvi[i])
	        *(eigveleft[j*dim+i]*eigveright[(i+1)*dim+l]
		  +eigveleft[j*dim+i+1]*eigveright[i*dim+l]))
		/(pow((1-eigvr[i]),2)+pow(eigvi[i],2));
	      G[l*dim+j]+=2*((eigvr[i]*(1-eigvr[i])-eigvi[i]*eigvi[i])
		*(eigveleft[j*dim+i]*eigveright[i*dim+l]
		  -eigveleft[j*dim+i+1]*eigveright[(i+1)*dim+l])
				      -eigvi[i]
	        *(eigveleft[j*dim+i]*eigveright[(i+1)*dim+l]
		  +eigveleft[j*dim+i+1]*eigveright[i*dim+l]))
		/(pow((1-eigvr[i]),2)+pow(eigvi[i],2));
	    }
	  }
	  i++;
	}
      }
    }
  }
  free(eq_m);
  return G;
}

double* fullBasis(int d){
  double* U=NULL;
  int dim=pow(2,d);
  int line=0;
  int j,k,l,c,vect,len;
  double norm,val;

  U=malloc(dim*dim*sizeof(double));
  //The matrix should be full of zeros
  for(j=0; j<dim*dim; j++){
    U[j]=0;
  }
  //If d=2 the computation is different
  if(d==2){
    U[1]=sqrt(0.5);
    U[2]=-sqrt(0.5);
    U[7]=1;
    U[8]=1;
    U[13]=sqrt(0.5);
    U[14]=sqrt(0.5);
    return U;
  }
  else{    
    for(j=0; j<d-2; j++){
      for(c=0; c<pow(2,j); c++){
	for(k=0; k<d-2-j; k++){
	  for(l=0; l<d-2-j; l++){
	    vect=pow(2,d-1)-pow(2,l+j+2)+c*pow(2,l+1)+pow(2,l)-1;
	    val=sin(M_PI*(k+1)*(l+1)/(d-1-j));
	    U[(line+k)*dim+(int)pow(2,d-1)+vect]+=val;
	    U[(line+k)*dim+2*vect+1]-=val;
	  }
	  norm=0;
	  for(l=0; l<dim; l++){
	    norm+=U[(line+k)*dim+l]*U[(line+k)*dim+l];
	  }
	  for(l=0; l<dim; l++){
	    U[(line+k)*dim+l]/=sqrt(norm);
	  }
	}
	line+=d-2-j;
      }
    }
    for(k=0; k<d-1; k++){
      for(l=0; l<d-1; l++){
	vect=pow(2,d-1)-1-pow(2,l);
	val=sin(M_PI*(k+1)*(l+1)/d);
	U[(line+k)*dim+(int)pow(2,d-1)+vect]+=val;
	U[(line+k)*dim+2*vect+1]-=val;
      }
      norm=0;
      for(l=0; l<dim; l++){
	norm+=U[(line+k)*dim+l]*U[(line+k)*dim+l];
      }
      for(l=0; l<dim; l++){
	U[(line+k)*dim+l]/=sqrt(norm);
      }
    }
    line+=d-1;
    U[line*dim+(int)pow(2,d)-1]=1;
    line+=1;
    for(c=0; c<pow(2,d-2); c++){
      U[(line+c)*dim+2*c]=1;
    }
    line+=pow(2,d-2);
    for(c=0; c<pow(2,d-3); c++){
      U[(line+c)*dim+2*c+(int)pow(2,d-1)]=sqrt(0.5);
      U[(line+c)*dim+4*c+1]=sqrt(0.5);
    }
    line+=pow(2,d-3);
    for(c=0; c<pow(2,d-3); c++){
      vect=2*c+pow(2,(d-1))+pow(2,(d-2));
      len=1;
      while(vect>=pow(2,d-1)){
	vect=(2*vect)%(int)pow(2,d)+vect/(int)pow(2,d-1);
	len++;
      }
      vect=2*c+pow(2,(d-1))+pow(2,(d-2));
      for(k=0; k<len; k++){
	U[(line+c)*dim+vect]=1./sqrt(len);
	vect=(2*vect)%(int)pow(2,d)+vect/(int)pow(2,d-1);
      }
    }
    return U;
  }
}

double* rotation(int d){
  double* V=NULL;
  int c,k,l,sign;
  int dim=pow(2,d);
  V=malloc(dim*dim*sizeof(double));
  for(c=0; c<dim; c++){
    for(k=0; k<dim; k++){
      sign=0;
      for(l=0; l<d; l++){
	sign+=(1-(c%(int)pow(2,l+1))/(int)pow(2,l))*(1-(k%(int)pow(2,l+1))/(int)pow(2,l));
      }
      V[c*dim+k]=pow(-1,sign)*pow(2,-0.5*d);
    }
  }
  return V;
}

double* reducedCurvature(double* wr, double* wi, double* vel, double* ver, int dim, double error){
  double* U=NULL;
  double* V=NULL;
  double* G=NULL;
  double* Gred=NULL;
  double* UV=NULL;
  double* UVG=NULL;
  int d,i,j,k;

  d=(log(0.5+dim)/log(2));
  G=curvature(wr,wi,vel,ver,dim,error);
  V=rotation(d);
  U=fullBasis(d);

  UV=malloc(dim*dim/2*sizeof(double));
  for(i=0; i<dim/2; i++){
    for(j=0; j<dim; j++){
      UV[i*dim+j]=0;
      for(k=0; k<dim; k++){
	UV[i*dim+j]+=U[(i+dim/2)*dim+k]*V[j*dim+k]; //V is symmetric
      }
    }
  }
  free(V);
  free(U);
  UVG=malloc(dim*dim/2*sizeof(double));
  for(i=0; i<dim/2; i++){
    for(j=0; j<dim; j++){
      UVG[i*dim+j]=0;
      for(k=0; k<dim; k++){
	UVG[i*dim+j]+=UV[i*dim+k]*G[j*dim+k]; //G is symmetric
      }
    }
  }
  free(G);
  Gred=malloc(dim*dim/4*sizeof(double));
  for(i=0; i<dim/2; i++){
    for(j=0; j<dim/2; j++){
      Gred[i*dim/2+j]=0;
      for(k=0; k<dim; k++){
	Gred[i*dim/2+j]+=UVG[i*dim+k]*UV[j*dim+k];
	//UV[i*dim/2+l]*G[k*dim+l]*UV[j*dim/2+k]
      }
    }
  }
  free(UV);
  free(UVG);
  return Gred;
}
