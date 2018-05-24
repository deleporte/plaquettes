
#if !defined(__COM_MODEL_H)
#define __COM_MODEL_H
//CE CODE C++ com.h EST ECRIT AUTOMATIQUEMENT PAR: makef.cc
//============================================================

# include <TCanvas.h> 
# include <TGFrame.h> 
# include <TGButton.h> 
# include <TGNumberEntry.h> 
# include <TGMenu.h> 
# include <TGLabel.h> 
# include <TApplication.h> 
# include <TGTab.h> 
# include <TGSlider.h> 
# include <TGProgressBar.h> 
# include <TGComboBox.h>  
# include <thread> 
# include <chrono> 
using namespace std::chrono; 
using namespace std::this_thread; 
# include <iostream> 
using namespace std; 
# include <mutex> 
# include "spin_lattice.h"

//==================
void Lance_com(Reseau *pReseau);

//===================
class Com: public  TGMainFrame
{
public :

	Com(Reseau *pReseau);
	TGNumberEntry *Reseau_N;
	TGNumberEntry *Reseau_verbose;
	TGCheckButton *Reseau_opt_sp;
	TGNumberEntry *Reseau_Nval;
	TGNumberEntry *Reseau_J1;
	TGNumberEntry *Reseau_J2;
	TGTextButton *Reseau_Dessin_spectre;
	TGNumberEntry *Reseau_N1;
	TGNumberEntry *Reseau_N2;
	TGTextButton *Reseau_Dessin_ecart;
	TGNumberEntry *Reseau_SIZE;
	TGTextButton *Reseau_Curvature_testing;

	Bool_t ProcessMessage(Long_t msg, Long_t p1, Long_t p2);
	void  Met_a_jour();
};

//============

extern Com *p_com;  // pointeur sur l'objet (declare dans .cc)
extern mutex mtx;

#endif
