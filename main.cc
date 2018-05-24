// -*- mode:C++ ; compile-command: "./compile .  main.cc; ./main" -*-



// ce qu'on aimerait comme bouton
// choisir SIZE(>=2)
// retourner curvature(...) ou reducedCurvature(...) (attention elles ne sont pas de la meme taille)
// afficher l'état d'équilibre (eq_m calculé au début de la fonction curvature)

#define JJ -0.9


#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <complex>
//#include <dynamics.h>

using namespace std;


#include "spin_lattice.h"


#include "com.h"






int main(int argc, char* argv[])
{

    TApplication theApp("App", nullptr, nullptr);

	Reseau reseau;
		
	thread t(Lance_com, &reseau);  // lance fenetre de commande dans autre thread, -> pointeur p_com

	t.join(); // attend ici que t se termine



	
//----donne la main a l'utilisateur ---
	theApp.Run();



	
	return 0;
}
