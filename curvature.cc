
#include <iostream>
#include <iomanip>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "curvature.h"

//lapack functions definitions

using namespace std;
using namespace arma;

void generate_C(int d, double m, double M, vec& vC){
  int i;

  vC.resize((int)pow(2,d));
  //C=malloc((int)pow(2,d)*sizeof(double));
  for(i=0; i<pow(2,d); i++){
    vC(i)=((double)rand())/RAND_MAX*(M-m)+m;
  }
  //return tableau(vC);
}

void Markov(vec& C, mat& T){
  int i;
  int dim=C.size();
  //cout<<dim<<endl;
  //T=malloc(dim*dim*sizeof(double));
  T.zeros(dim,dim);
  for(i=0; i<dim; i++){
    T((2*i)%dim,i)=C((2*i)%dim)/(C((2*i)%dim)+C((2*i)%dim+1));
    T((2*i)%dim+1,i)=C((2*i)%dim+1)/(C((2*i)%dim)+C((2*i)%dim+1));
  }
  //return T;
}

void diag(vec& C, cx_vec& w, cx_mat& vel, cx_mat& ver, mat& T)
{
  Markov(C,T); //computing Markov chain
  //cout<<"Markov matrix computed"<<endl;
  eig_gen(w,vel,T); //eigenvalues and eigenvectors
  //uvec indices=sort_index(abs(w),"descend"); //sort eigenvalues according to their abs value
  //w=w(indices);
  //ver=ver.cols(indices); //sort right eigenvectors accordingly
  //cout<<det(vel)<<endl;
  if(!inv(ver,vel)){
    cout<<"Error while inverting eigenmatrix of the Markov chain"<<endl;
    cout<<"Determinant: "<<real(det(vel))<<endl;
  }
  //cout<<"Diagcheck:"<<vel*diagmat(w)*ver-T<<endl;
  //eigvals(T,dim,wr,wi,vel,ver); //à changer pour utiliser plutot armadillo
  //return T;
}

void curvature(cx_vec& w, cx_mat& vel, cx_mat& ver, mat& G){
  uvec indices=sort_index(abs(w),"descend"); //sort eigenvalues according to their abs value
  vec eq_m=abs(vel.col(indices(0)));
  int i,j,k,l;
  double norm;

  //double* G=NULL;
  int dim=w.size();
  norm=0;
  for(i=0; i<dim; i++){
    norm+=eq_m(i);
  }
  for(i=0; i<dim; i++){
    eq_m(i)/=norm;
  }
  //cout<<"equilibrium measure: "<<endl<<eq_m<<endl;
  //eq_m=malloc(dim*sizeof(double));
  // if(1){
  //   if(1){
  //     //find the equilibrium state
  //     for(i=0; i<dim; i++){
  // 	if(fabs(eigvi(i))<0.00001 && eigvr(i)>0.9999){
  // 	  for(j=0; j<dim; j++){
  // 	    eq_m(j)=eigveright(i*dim+j);
  // 	  }
  // 	  //positive
  // 	  if(eq_m(0)<0){
  // 	    for(j=0; j<dim; j++){
  // 	      eq_m(j)=-eq_m(j);
  // 	    }
  // 	  }
  // 	  //normalize
  // 	  norm=0;
  // 	  for(j=0; j<dim; j++){
  // 	    norm+=eq_m(j);
  // 	  }
  // 	  for(j=0; j<dim; j++){
  // 	    eq_m(j)/=norm;
  // 	  }
  // 	}
  //     }
  //   }
  // }
    
  //G=malloc(dim*dim*sizeof(double));
  G.zeros(dim,dim);
  
  // for(i=0; i<dim*dim; i++){
  //   if(i%(dim+1)){
  //     G(i)=-eq_m(i%dim)*eq_m(i/dim);
  //   }
  //   else{
  //     G(i)=eq_m(i%dim)*(1-eq_m(i/dim));
  //   }
  // }
  for(i=0; i<dim; i++){
    G(i,i)=eq_m(i);
  }
  for(i=0; i<dim; i++){
    for(j=0; j<dim; j++){
      G(i,j)-=eq_m(i)*eq_m(j);
    }
  }
  
  for(i=0; i<dim; i++){
    if(i!=indices(0)){
      for(j=0; j<dim; j++){
	for(l=0; l<dim; l++){
	  G(j,l)+=real(w(i)/(1.-w(i))*eq_m(j)*vel(l,i)*ver(i,j));
	  G(l,j)+=real(w(i)/(1.-w(i))*eq_m(j)*vel(l,i)*ver(i,j)); //on peut etre + efficaces ?
	}
      }
    }
  }
  //cout<<"Curvature computed !"<<endl;
}

void fullBasis(int d, mat& U){
  int dim=pow(2,d);
  int line=0;
  int j,k,l,c,vect,len;
  double norm,val;

  U.zeros(dim,dim);

  //U=malloc(dim*dim*sizeof(double));
  //The matrix should be full of zeros
  //for(j=0; j<dim*dim; j++){
  //U(j)=0;
  //}
  //If d=2 the computation is different
  if(d==2){
    U(0,1)=sqrt(0.5);
    U(0,2)=-sqrt(0.5);
    U(1,3)=1;
    U(2,0)=1;
    U(3,1)=sqrt(0.5);
    U(3,2)=sqrt(0.5);
    //return U;
  }
  else{    
    for(j=0; j<d-2; j++){
      for(c=0; c<pow(2,j); c++){
	for(k=0; k<d-2-j; k++){
	  for(l=0; l<d-2-j; l++){
	    vect=pow(2,d-1)-pow(2,l+j+2)+c*pow(2,l+1)+pow(2,l)-1;
	    val=sin(M_PI*(k+1)*(l+1)/(d-1-j));
	    U((line+k),(int)pow(2,d-1)+vect)+=val;
	    U((line+k),2*vect+1)-=val;
	  }
	  norm=0;
	  for(l=0; l<dim; l++){
	    norm+=U((line+k),l)*U((line+k)*dim+l);
	  }
	  for(l=0; l<dim; l++){
	    U((line+k),l)/=sqrt(norm);
	  }
	}
	line+=d-2-j;
      }
    }
    for(k=0; k<d-1; k++){
      for(l=0; l<d-1; l++){
	vect=pow(2,d-1)-1-pow(2,l);
	val=sin(M_PI*(k+1)*(l+1)/d);
	U((line+k),(int)pow(2,d-1)+vect)+=val;
	U((line+k),2*vect+1)-=val;
      }
      norm=0;
      for(l=0; l<dim; l++){
	norm+=U((line+k),l)*U((line+k),l);
      }
      for(l=0; l<dim; l++){
	U((line+k),l)/=sqrt(norm);
      }
    }
    line+=d-1;
    U(line*dim+(int)pow(2,d)-1)=1;
    line+=1;
    for(c=0; c<pow(2,d-2); c++){
      U((line+c),2*c)=1;
    }
    line+=pow(2,d-2);
    for(c=0; c<pow(2,d-3); c++){
      U((line+c),2*c+(int)pow(2,d-1))=sqrt(0.5);
      U((line+c),4*c+1)=sqrt(0.5);
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
	U((line+c),vect)=1./sqrt(len);
	vect=(2*vect)%(int)pow(2,d)+vect/(int)pow(2,d-1);
      }
    }
    //return U;
  }
  U.shed_rows(0,(int)pow(2,d-1)-1); //needs to be made more efficiently
}

void rotation(int d,mat& V){
  int c,k,l,sign;
  int dim=pow(2,d);
  V.zeros(dim,dim);
  for(c=0; c<dim; c++){
    for(k=0; k<dim; k++){
      sign=0;
      for(l=0; l<d; l++){
	sign+=(1-(c%(int)pow(2,l+1))/(int)pow(2,l))*(1-(k%(int)pow(2,l+1))/(int)pow(2,l));
      }
      V(c,k)=pow(-1,sign)*pow(2,-0.5*d);
    }
  }

  //return V;
}

void reducedCurvature(cx_vec& w, cx_mat& vel, cx_mat& ver, mat& Gred){
  mat G,V,U;
  int d,i,j,k;
  int dim=w.size();

  d=(log(0.5+dim)/log(2));
  curvature(w,vel,ver,G);
  rotation(d,V);
  fullBasis(d,U);

  // UV=malloc(dim*dim/2*sizeof(double));
  // for(i=0; i<dim/2; i++){
  //   for(j=0; j<dim; j++){
  //     UV(i*dim+j)=0;
  //     for(k=0; k<dim; k++){
  // 	UV(i*dim+j)+=U((i+dim/2)*dim+k)*V(j*dim+k); //V is symmetric
  //     }
  //   }
  // }
  // free(V);
  // free(U);
  // UVG=malloc(dim*dim/2*sizeof(double));
  // for(i=0; i<dim/2; i++){
  //   for(j=0; j<dim; j++){
  //     UVG(i*dim+j)=0;
  //     for(k=0; k<dim; k++){
  // 	UVG(i*dim+j)+=UV(i*dim+k)*G(j*dim+k); //G is symmetric
  //     }
  //   }
  // }
  // free(G);
  // Gred=malloc(dim*dim/4*sizeof(double));
  // for(i=0; i<dim/2; i++){
  //   for(j=0; j<dim/2; j++){
  //     Gred(i*dim/2+j)=0;
  //     for(k=0; k<dim; k++){
  // 	Gred(i*dim/2+j)+=UVG(i*dim+k)*UV(j*dim+k);
  // 	//UV(i*dim/2+l)*G(k*dim+l)*UV(j*dim/2+k)
  //     }
  //   }
  // }
  // free(UV);
  // free(UVG);
  Gred=U*(V.t())*G*V*(U.t());
  //Gred=G;
}
