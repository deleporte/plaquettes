#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "lapack/lapacke.h"
#include "curvature.h"

double* generate_C(int d, double m, double M){
  double* C=NULL;
  int i;

  C=malloc((int)pow(2,d)*sizeof(double));
  for(i=0; i<pow(2,d); i++){
    C[i]=((double)rand())/RAND_MAX*(M-m)+m;
  }
  return C;
}

double* sparse_mult(double* mat, double* C, int dim){
  double* prod=NULL;
  int i=0;
  int j=0;

  prod=malloc(dim*dim*sizeof(double));
  for(i=0; i<dim; i++){
    for(j=0; j<dim; j++){
      prod[i*dim+j]=mat[i*dim+(2*j)%dim]*C[(2*j)%dim];
      prod[i*dim+j]+=mat[i*dim+(2*j+1)%dim]*C[(2*j+1)%dim];
      prod[i*dim+j]/=C[(2*j)%dim]+C[(2*j+1)%dim];
      //printf("Value on line %d and column %d: %f\n", i,j,prod[i*dim+j]);
    }
  }
  return prod;
}

double** Markov_powers(double* C, int dim, double error){
  double** powerseq;
  char converged=0;
  int i;
  int j;
  double dist;

  powerseq=malloc(100*sizeof(double*));
  //The first matrix is the identity
  powerseq[0]=malloc(dim*dim*sizeof(double));
  for(i=0; i<dim; i++){
    for(j=0; j<dim; j++){
      if(i==j){
	powerseq[0][i*dim+j]=1;
      }
      else{
	powerseq[0][i*dim+j]=0;
      }
    }
  }
  //Test: writing the identity matrix
  
  i=0;
  while (converged==0 && i<99){
    i++;
    powerseq[i]=sparse_mult(powerseq[i-1],C,dim);
    dist=0;
    for(j=0; j<dim; j++){
      dist+= fabs(powerseq[i][j*dim]-powerseq[i][j*dim+1]);
    }
    if (dist<error){
      converged=1;
      //printf("Markov chain converges after %d steps.\n", i);
    }
    if (i==99){
      //printf("Warning: slow convergence.\n");
    }
  }
  for (j=i; j<100; j++){
    powerseq[j]=powerseq[i];
  }
  return powerseq;
}

int isDangerous(double* C){
  return 0;
}

double* curvature(double* C, int dim, double error){
  double** powerseq=NULL;
  double* eq_m=NULL;
  double* G=NULL;
  double* eigveright=NULL;
  double* eigveleft=NULL;
  double* eigvr=NULL;
  double* eigvi=NULL;
  double* Q=NULL;
  double* Z=NULL;
  double* Id=NULL;
  double* H=NULL;
  double* Hmem=NULL;
  lapack_logical* select=NULL;
  double* workspace=NULL;
  int* failedr=NULL;
  int* failedl=NULL;
  int i,j,k,l;
  int layout=LAPACK_ROW_MAJOR;
  double norm,r,th;
  powerseq=Markov_powers(C,dim,error);
  eq_m=malloc(dim*sizeof(double));
  if(isDangerous(C)){
    //spectral decomposition of C
    Id=malloc(dim*dim*sizeof(double));
    for(i==0; i<dim*dim; i++){
      if(i%(dim+1)==0){
	Id[i]=1;
      }
      else{
	Id[i]=0;
      }
    }
    Q=malloc(dim*dim*sizeof(double));
    Z=malloc(dim*dim*sizeof(double));
    Hmem=malloc(dim*dim*sizeof(double));
    for(i=0; i<dim*dim; i++){
      Hmem[i]=powerseq[1][i];
    }
    LAPACKE_dgghrd(layout,'I','I',dim,1,dim,Hmem,dim,Id,dim,Q,dim,Z,dim);
    free(Id);
    H=malloc(dim*dim*sizeof(double));
    for(i=0; i<dim*dim; i++){
      H[i]=Hmem[i];
    }
    LAPACKE_dhseqr(layout,'E','N',dim,1,dim,H,dim,eigvr,eigvi,NULL,1);
    for(i=0; i<dim*dim; i++){
      H[i]=Hmem[i];
    }
    select=malloc(dim*sizeof(lapack_logical));
    
    //find eigenvectors with large eigenvalues
    k=0;
    for(i=0; i<dim; i++){
      if(eigvr[i]*eigvr[i]+eigvi[i]*eigvi[i]>0.95){
        select[i]=1;
	k++;
      }
      else{
	select[i]=0;
      }
    }
    
    eigveright=malloc(k*dim*sizeof(double));
    eigveleft=malloc(k*dim*sizeof(double));
    failedr=malloc(k*sizeof(int));
    failedl=malloc(k*sizeof(int));
    LAPACKE_dhsein(layout, //int matrix_layout
		   'B',    //char job
		   'N',    //char eigsrc
		   'N',    //char initv
		   select, //lapack_logical* select
		   dim,    //lapack_int n
		   H,      //double* h
		   dim,    //int ldh
		   eigvr,  //double* wr
		   eigvi,  //double* wi
		   eigveleft, //double* vl
		   dim,    //int ldvl
		   eigveright, //double* vr
		   dim,     //int ldvr
		   k,       //int mm
		   NULL,       //int* m
		   failedl,   //int* iffaill
		   failedr);   //int* iffailr
    free(workspace);
    free(select);
    free(H);
    workspace=malloc(dim*sizeof(double));
    k=0;
    //multiply by Q the left eigenvectors and by Z' the right eigenvectors
    for(i=0; i<dim; i++){
      if(eigvr[i]*eigvr[i]+eigvi[i]*eigvi[i]>0.95){
	for(j=0; j<dim; j++){
	  workspace[j]=0;
	  for(l=0; l<dim; l++){
	    workspace[j]+=Q[j*dim+l]*eigveleft[k*dim+l];
	  }
	}
	for(j=0; j<dim; j++){
	  eigveleft[k*dim+j]=workspace[j];
	}
	for(j=0; j<dim; j++){
	  workspace[j]=0;
	  for(l=0; l<dim; l++){
	    workspace[j]+=Z[j*dim+l]*eigveright[k*dim+l];
	  }
	}
	for(j=0; j<dim; j++){
	  eigveright[k*dim+j]=workspace[j];
	}
	k++;
      }
    }
    free(Q);
    free(Z);
    free(workspace);
    free(Hmem);
    free(failedr);
    free(failedl);
    //normalize in linfty the right eigenvectors and l1 the left eigenvectors
    k=0;
    for(i=0; i<dim; i++){
      if(eigvr[i]*eigvr[i]+eigvi[i]*eigvi[i]>0.95){
	//if the eigenvalue is real, there is only one column
	if(abs(eigvi[i])<0.00001){
	  norm=0;
	  for(j=0; j<dim; j++){
	    if (abs(eigveright[k*dim+j])>norm){
	      norm=abs(eigveright[k*dim+j]);
	    }
	  }
	  for(j=0; j<dim; j++){
	    eigveright[k*dim+j]/=norm;
	  }
	  norm=0;
	  for(j=0; j<dim; j++){
	    norm+= abs(eigveleft[k*dim+j]);
	  }
	  for(j=0; j<dim; j++){
	    eigveleft[k*dim+j]/=norm;
	  }
	  k++;
	}
	else{
	  norm=0;
	  for(j=0; j<dim; j++){
	    if (sqrt(eigveright[k*dim+j]*eigveright[k*dim+j]
		     +eigveright[(k+1)*dim+j]*eigveright[(k+1)*dim+j])>norm){
	      norm=sqrt(eigveright[k*dim+j]*eigveright[k*dim+j]
			+eigveright[(k+1)*dim+j]*eigveright[(k+1)*dim+j]);
	    }
	  }
	  for(j=0; j<2*dim; j++){
	    eigveright[k*dim+j]/=norm;
	  }
	  norm=0;
	  for(j=0; j<dim; j++){
	    norm+= sqrt(eigveleft[k*dim+j]*eigveleft[k*dim+j]
			+eigveleft[(k+1)*dim+j]*eigveleft[(k+1)*dim+j]);
	  }
	  for(j=0; j<2*dim; j++){
	    eigveleft[k*dim+j]/=norm;
	  }
	  k+=2;
	}
      }
    }
    //find the equilibrium state
    k=0;
    for(i=0; i<dim; i++){
      if(eigvr[i]*eigvr[i]+eigvi[i]*eigvi[i]>0.95){
	k++;
	if(abs(eigvi[i])<0.00001 && eigvr[i]>0.99999){
	  for(j=0; j<dim; j++){
	    eq_m[j]=eigveleft[k*dim+j];
	  }
	}
      }
    }
  }
  else{
    for(i=0; i<dim; i++){
      eq_m[i]=powerseq[99][i*dim];
    }
  }
  G=malloc(dim*dim*sizeof(double));
  //in the general case, G depends only on the equilibrium measure
  //and the data of the first powers
  for(i=0; i<dim; i++){
    for(j=0; j<dim; j++){
      G[i*dim+j]=0;
      for(k=0; k<100; k++){
	G[i*dim+j]+=eq_m[j]*(powerseq[k][i*dim+j]-eq_m[i]);
	if (k!=0){
	  G[i*dim+j]+=eq_m[i]*(powerseq[k][j*dim+i]-eq_m[j]);
	}
      }
    }
  }
  if(isDangerous(C)){
    k=0;
    for(i=0; i<dim; i++){
      if(eigvr[i]*eigvr[i]+eigvi[i]*eigvi[i]>0.95){
	if(abs(eigvi[i])<0.00001){
	  for(j=0; j<dim; j++){
	    for(l=0; l<dim; l++){
	      G[j*dim+l]+=pow(eigvr[k],100)/(1-eigvr[k])*eq_m[l]
		*(eigveleft[k*dim+j]*eigveright[k*dim+l]);
	      G[l*dim+j]+=pow(eigvr[k],100)/(1-eigvr[k])*eq_m[l]
		*(eigveleft[k*dim+j]*eigveright[k*dim+l]);
	    }
	  }
	  k++;
	}
	else{
	  for(j=0; j<dim; j++){
	    for(l=0; l<dim; l++){
	      r=sqrt(eigvr[k]*eigvr[k]+eigvi[k]*eigvi[k]);
	      th=acos(eigvr[k]/r);
	      G[j*dim+l]+=2*pow(r,100)*((cos(100*th)*(1-eigvr[k])
				      -sin(100*th)*eigvi[k])
		*(eigveleft[k*dim+j]*eigveright[k*dim+l]
		  -eigveleft[(k+1)*dim+j]*eigveright[(k+1)*dim+l])
				      +(sin(100*th)*eigvr[k]
					+cos(100*th)*(1-eigvr[k]))
	        *(eigveleft[k*dim+j]*eigveright[(k+1)*dim+l]
		  -eigveleft[(k+1)*dim+j]*eigveright[k*dim+l]))
		/(pow((1-eigvr[k]),2)+pow(eigvi[k],2));
	      G[l*dim+j]+=2*pow(r,100)*((cos(100*th)*(1-eigvr[k])
				      -sin(100*th)*eigvi[k])
		*(eigveleft[k*dim+j]*eigveright[k*dim+l]
		  -eigveleft[(k+1)*dim+j]*eigveright[(k+1)*dim+l])
				      +(sin(100*th)*eigvr[k]
					+cos(100*th)*(1-eigvr[k]))
	        *(eigveleft[k*dim+j]*eigveright[(k+1)*dim+l]
		  -eigveleft[(k+1)*dim+j]*eigveright[k*dim+l]))
		/(pow((1-eigvr[k]),2)+pow(eigvi[k],2));
	    }
	  }
	  k+=2;
	}
      }
    }
  }
  
  for(k=0; k<100; k++){
    if(k==99){
      free(powerseq[k]);
    }
    else{
      if (powerseq[k]==powerseq[k+1]){
	k=99;
      }
      free(powerseq[k]);
    }
  }

  

  free(powerseq);
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

double* reducedCurvature(double* C, int dim, double error){
  double* U=NULL;
  double* V=NULL;
  double* G=NULL;
  double* Gred=NULL;
  double* UV=NULL;
  double* UVG=NULL;
  int d,i,j,k,l,m,n;
  double val;

  d=(log(0.5+dim)/log(2));
  G=curvature(C,dim,error);
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
