#if !defined(_SPIN_LATTICE_H)
#define _SPIN_LATTICE_H


#include <iostream>
using namespace std;

#include <armadillo>
using namespace arma;


#include <vector>

typedef long long int l_int;
typedef complex<double> Complex;

#include <TCanvas.h>

#include <TMarker.h>
#include <TH1D.h>
#include <TApplication.h>
#include <boost/math/special_functions/binomial.hpp>
using namespace boost::math;



#include <bitset>
#define I_ Complex(0,1)

//================================
class Reseau
{
public:


	//... documentation
	void documentation_Fred_en_lyx(); // make_gui = M("Help","spins-reseau.lyx")
	void documentation_Fred_en_pdf(); // make_gui = M("Help","spins-reseau.pdf")
	void documentation_Alix_en_pdf(); // make_gui = M("Help","curvature.pdf")

	//.....................
	
	int N = 4;  //  make_gui =   N(ZC, ":N nombre de sites") help = "conseil N<=10"

	int verbose =0;  //  make_gui =   N(ZC, ":verbose") 

	//... definition de l'interaction, Hamiltonien sur reseau
	vector<vector<int>> E = {{3,3} , {1}}; // multiindices
	vector<double> A = {1,0} ; // amplitudes correspondantes


	//.. modele total....
	cx_mat H; // matrice du Hamiltonien dans l'espace total (taille 2^N * 2^N)
	sp_cx_mat Hsp; //idem sparse

	int opt_sp = 0; //  make_gui = nl  C(ZT("Spectre"), ": utilise matrice sparse")


	int Nval= 5*N+1; //  make_gui =   N(ZT("Spectre"), ":(si sparse) nbre de valeurs propres demandees") 

	
	//....
	void Remplit_matrice_H();
	void Calcule_ep_C(l_int e, vector<int> alpha,int i, l_int & ep, Complex & C);



// ---------... modele approche , etats constants
	mat Hcste; // Matrice


	void Remplit_matrice_Hcste();

	double J1=0;  // make_gui = nl  N(ZT("Spectre"), ":J1") 
	double J2=1;   // make_gui =   N(ZT("Spectre"), ":J2") 
	void Dessin_spectre(); // make_gui =   B(ZT("Spectre"), "Dessin du spectre") help = "Dessin du spectre. Donner J1,J2"

	int N1=2;  // make_gui = nl  N(ZT("Spectre"), ":N1") 
	int N2=10;  // make_gui =   N(ZT("Spectre"), ":N1") 
	void Dessin_ecart();  // make_gui =   B(ZT("Spectre"), "Dessin des ecarts") help = "Dessin des ecarts. Donner N1,N2"

	//---- Methode des fermions
	vec Spectre_exact(double J);



//-----  methode de Markov

	vec C;
	cx_vec w;
	mat Gred,T,G;
	cx_mat vel,ver;
	vec Cr,Ci;
	cx_vec Cc;
	//cx_mat H;

	int SIZE = 2;  //  make_gui =   N(ZT("Curvature"), ": SIZE") help = "conseil SIZE>=2"
	



	void Curvature_testing();  // make_gui =   B(ZT("Curvature"), "Curvature testing") help = "test. Donner ..."
	void Dynamics_testing(); // make_gui = B(ZT("Dynamics"), "Dynamics testing") help "Test."



};

#endif
