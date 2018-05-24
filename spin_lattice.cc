#include "curvature.h"
#include "dynamics.h"

#include "spin_lattice.h"


#include <stdlib.h> // pour system




//======
/*
eps=0,1
alpha=0,1,2,3
sigma[alpha][eps] = 0,1 : image du spin
TC[alpha][eps] coef complexe
*/
vector<vector<int>>  sigma = { {0,1} , {1,0} , {1,0}, {0,1} }; // modif du spin
vector<vector<Complex>> TC = { {1,1} , {1,1} , {I_, -I_}, {1,-1}};




//==========================
// Affiche vector
ostream & operator <<(ostream & sortie,  const vector<int> & V)
{

	sortie<<"["<<flush;

	for(int i=0;i<V.size();i++)
    {
		if(i>0)
			sortie<<","<<flush;
		sortie<<V[i]<<flush;
    }
	cout<<"]"<<endl;

    return sortie;
}

//=====
/*
Convertit entier >0 en ecriture binaire
bit faible est en position i=0.
 */
vector<int> int_to_binary(l_int e)
{
	vector<int> e_b;
	if(e<0)
		return e_b;

	while(e>0)
	{
		int b= e%2;
		e=e/2;
//		cout<<"b="<<b<<" e="<<e<<endl;
		e_b.push_back(b);
	}
	return e_b;
}


//=====
/*!
Convertit  ecriture binaire en entier
bit faible est en position i=0.
 */
l_int binary_to_int(vector<int> be)
{
	l_int e=0;
	for(int i=be.size()-1; i>=0; i--)
	{
		e=2*e;
		e= e+ be[i];
	}
	return e;
}



//============================
void Reseau::documentation_en_lyx()
{
	int res =system("lyx rapport/spins-reseau.lyx");
}

//============================
void Reseau::documentation_en_pdf()
{
int res =	system("evince rapport/spins-reseau.pdf");
}


//================================
/*!
  Compute ep and C
  -------------------

 - entree: e, alpha, i 
 - sortie: ep,C

 */
void  Reseau::Calcule_ep_C(l_int e, vector<int> alpha,int i, l_int & ep, Complex & C)
{
	vector<int> be = int_to_binary(e);

	while(be.size()<N)
		be.push_back(0);

	vector<int> bep = be;

	int m=alpha.size();

	C=1; 

	for(int j=0; j<= m-1; j++)
	{
		int j2=(i+j)%N;

		bep[j2] = sigma[alpha[j]][be[j2]]; 

		C = C *TC[alpha[j]][be[j2]]; 
	}


	ep = binary_to_int(bep);




	if(verbose>=2)
	{
		cout<<"----------------"<<endl;
		cout<<" e="<<e<<" be="<<be<<endl;
		cout<<" alpha="<<alpha<<" i="<<i<<endl;
		cout<<" bep="<<bep<<" ep="<<ep<<" C="<<C<<endl;

	}

}


//================================
void Reseau::Remplit_matrice_H()
{
	if(verbose)
		cout<<"remplit matrice H.."<<flush;
	if(opt_sp == 0)
		H.zeros(pow(2,N), pow(2,N));
	else
		Hsp = sp_cx_mat(pow(2,N), pow(2,N) ); // taille matrice
	

	for(l_int e=0; e<pow(2,N); e++) // elt de matrice epsilon
	{

		for(int ia =0; ia < E.size(); ia++) // indice de  alpha
		{  
			auto alpha = E[ia];
			double Ai = A[ia];


			for(int i =0; i<N; i++) // decalage
			{

				l_int ep;
				Complex C;
				Calcule_ep_C(e,alpha,i, ep, C); // -> ep, C

				
				
				if(opt_sp == 0)
					H(ep,e) += Ai * C; 
				else
					Hsp(ep,e) += Ai * C; 


			}// for i
		} // for alpha

	}// for e


//-----
	if(verbose)
		cout<<" OK."<<endl;


}


//================================
void Reseau::Remplit_matrice_Hcste()
{
	Hcste.zeros(N+1,N+1);

	for(int n=0;n <N+1;n++)
	{
		//.. terme Ising
		Hcste(n,n) +=  A[0]* (N*N- (4*n+1)*N +4 * n*n) / (N-1);


		//.. terme ferrom
		if((n-1)>=0)
			Hcste(n-1,n) += A[1] * sqrt(n * (N-n+1)); 

		if((n+1)<= N)
			Hcste(n+1,n) += A[1] * sqrt((N-n) * (n+1)); 


	}

}


/*====================================
entree: J
sortie: spectre exact par la methode des fermions
 */
vec Reseau::Spectre_exact(double J)
{

	vec VE(pow(2,N));

//... calcul du spectre a une particule (1 fermion). Formule de Sachdev page 139.

	vec S0(N), S1(N);
	double g = (1-J)/J;
	
	for(int k=0;k<N;k++)
	{
		S0(k) = 2*J * sqrt(1+g*g-2*g*cos((2*M_PI*(k+0.5)  )/double(N)));
		S1(k) = 2*J * sqrt(1+g*g-2*g*cos((2*M_PI*k )/double(N)));

		/*
		if((2*k) == N)
			S0(k) =  2 *(1-2*J);

		
		if((2*k) == (N-1))
			S1(k) =  2 *(1-2*J);
		*/
	}
	
	//.... On deduit le spectre total par remplissage

	
	for(int n=0; n<pow(2,N); n++) // remplissage
	{
		double E=0;
		int n0=n; //copie
		//.......
		int parite = 0;
		for(int k=0;k<N;k++) //decalage
		{
			parite += n0%2;
			n0/=2;
		}
		
		parite = parite%2;
		
		//............
		n0=n; //copie
		if(parite==0)
		{
			for(int k=0;k<N;k++) //decalage
			{
			
				E+= S0(k) * ((n0%2) -0.5 );
				n0=n0/2;
			}
		}
		else
			for(int k=0;k<N;k++) //decalage
			{
			
				E+= S1(k) * ((n0%2) -0.5 );
				n0=n0/2;
			}

		
		VE(n)=E;
	}

	return VE;
}

//================================
void Reseau::Dessin_spectre()
{
	int NJ=500; // discretisation de J

	TCanvas *c = new TCanvas("Spectre","Spectre",0,0, 500, 500);
	TH1F *h= c->DrawFrame(J1,-N,J2,N);
	h->SetXTitle("J");
	h->SetYTitle("E");

	for(int iJ=NJ-1;iJ>=0; iJ--)
	{
	
		double	J= (double)iJ/NJ*(J2-J1)+J1;
		cout<<"J="<<J<<endl;
		A[0]=-J; 
		A[1]=-(1-J);

		//.... Hamiltonien Total 
		Remplit_matrice_H();


        //-------


		if(verbose)
			cout<<"Diagonalise.."<<flush;

		cx_vec val_p;
		cx_mat vect_p;


		if(opt_sp == 0)
		{
			val_p = eig_gen( H);
			uvec indices = sort_index(real(val_p)); // liste des indices
			val_p=val_p(indices); // vecteur ordonné
		}

		else
			eigs_gen(val_p, vect_p, Hsp, Nval, "sr"); // sr: smallest real part


		for(int i=0;i<val_p.size(); i++)
		{
			TMarker *p = new TMarker(J,real(val_p[i]),6);
			p->SetMarkerColor(kRed);
			p->Draw();
		}

		if(verbose)
			cout<<"OK."<<endl;

		//... Hamiltonien approche
		/*
		Remplit_matrice_Hcste();
		vec val_p2 = eig_sym(Hcste);

		for(int i=0;i<val_p2.size(); i++)
		{
			double y =   val_p2[i];
//			double y =  val_p2[i]* real(val_p[0]) / val_p2[0] ;
			TMarker *p = new TMarker(J,y,8);
			p->SetMarkerColor(kBlue);
			//	p->Draw();
		}
		*/


		//... formule approchée de la premiere bande coté J=1
		/*
		for(int k1=0; k1<N; k1++)
			for(int k2=0; k2<N; k2++)
			{
				double E= -J*(N-4) + 2*(1-J)* (cos(2*M_PI*k1/(double(N))) + cos(2*M_PI*k2/(double(N))) );

				TMarker *p = new TMarker(J,E ,6);
				p->SetMarkerColor(kBlue);
				p->Draw();
			}
			
		*/


		//... spectre exact par la methode des Fermions


		vec E = Spectre_exact(J);
		for(int i=0;i<E.size(); i++)
		{
			TMarker *p = new TMarker(J,E(i),6);
			p->SetMarkerColor(kBlue);
			p->Draw();
		}
		

		//......
		if(iJ%100==0)
			c->Update();
	}// for iJ


}


//================================
void Reseau::Dessin_ecart()
{
	TCanvas *c = new TCanvas("Spectre","Spectre",500,500);
	TH1F *h= c->DrawFrame(N1,-10,N2,10);
	h->SetXTitle("N");
	h->SetYTitle("E");

	for(N=N1; N <= N2; N++)
	{
		double J=0.3;

		A[0]=-J; 
		A[1]=-(1-J);

		//.... Hamiltonien Total 
		Remplit_matrice_H();



		if(verbose)
			cout<<"Diagonalise.."<<flush;
		
		cx_vec val_p;

		if(opt_sp == 0)
		{
			val_p = eig_gen( H);
			uvec indices = sort_index(real(val_p)); // liste des indices
			val_p=val_p(indices); // vecteur ordonné
		}
		else
			eigs_gen(val_p,  Hsp, Nval, "sr"); // sr: smallest real part
			
		if(verbose)
			cout<<"OK."<<endl;


		//... Hamiltonien approche
		Remplit_matrice_Hcste();
		vec val_p2 = eig_sym(Hcste);


		double dE = val_p2[1] - real(val_p[1]);

		TMarker *p = new TMarker(N, dE, 8);
		p->SetMarkerColor(kRed);
		p->Draw();
		
		//......
		c->Update();
	}// for N


}




//=============================
void Reseau::Curvature_testing()
{

	
	srand(time(NULL));
	generate_C(SIZE,0.1,1.1,C);
	for(int i=0; i<pow(2,SIZE); i++){
		printf("C[%d]=%f\n",i,C[i]);
	}
	cout<<"C initialized."<<endl;
	diag(C,w,vel,ver,T);
	cout<<"Markov matrix diagonalized."<<endl;
	reducedCurvature(w,vel,ver,Gred);
	curvature(w,vel,ver,G);
	cout<<"G="<<G<<endl;
	cout<<"spectrum="<<eig_sym(G)<<endl;
	cout<<"Gred="<<Gred<<endl;
	cout<<"spectrum="<<eig_sym(Gred)<<endl;
}

void Reseau::Dynamics_testing()
{
  srand(time(NULL));
  generate_C(SIZE,0.1,1.1,Cr);
  generate_C(SIZE,0.1,1.1,Ci);
  Cc=Cr+I_*Ci;
  H.zeros(4,4);
  H(0,0)=1;

  oneStep(Cc,H,0.001);
}
