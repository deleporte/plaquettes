//CE CODE C++ com.cc EST ECRIT AUTOMATIQUEMENT PAR: makef.cc
//========================================================
 
# include "./com.h" 


//==================
Reseau *p_Reseau; // pointeur sur objet de la classe Reseau
mutex mtx; 
Com *p_com; 

//=================== 
// Lance la fenetre de commandes 
// 
void Lance_com(Reseau *pReseau)
{

 	sleep_for(milliseconds{100}); // attente pour etre sur que un autre root a le temps de se lancer si besoin. 
 	TApplication theApp("App2", nullptr, nullptr); 
	p_com=  new Com(pReseau);  // cree fenetre de control et l'associe a  l'objet present
 	theApp.Run(); //----donne la main a l'utilisateur --- 
}

//====== identifieurs (int) pour les messages de commandes
enum Command_id 
{
Id_Reseau_M_documentation_en_lyx, 
Id_Reseau_M_documentation_en_pdf, 
Id_Reseau_N_N, 
Id_Reseau_N_verbose, 
Id_Reseau_C_opt_sp, 
Id_Reseau_N_Nval, 
Id_Reseau_N_J1, 
Id_Reseau_N_J2, 
Id_Reseau_B_Dessin_spectre, 
Id_Reseau_N_N1, 
Id_Reseau_N_N2, 
Id_Reseau_B_Dessin_ecart, 
Id_Reseau_N_SIZE, 
Id_Reseau_B_Curvature_testing, 
Id_Reseau_B_Dynamics_testing
}; 

int X=300, Y=150; // taille en pixels

//===================
Com::Com(Reseau *pReseau)
 : TGMainFrame(gClient->GetRoot(), X,Y) 
{


	//-- bloque le mutex
	mtx.lock(); // si mtx est libre on le bloque, sinon on attend


	//---- pointeurs sur objets
     p_Reseau = pReseau;


	//-- regles de positionnement
	TGLayoutHints *fLH_TX =  new TGLayoutHints(kLHintsTop | kLHintsExpandX , 2, 2, 2, 5);
	TGLayoutHints *fLH_C = new TGLayoutHints(kLHintsCenterX | kLHintsExpandX | kLHintsExpandY,5,5,5,5);
	TGLayoutHints *fLH_L = new TGLayoutHints(kLHintsLeft,5,5,5,5);
	TGLayoutHints *fLH_R = new TGLayoutHints(kLHintsRight,5,5,5,5);



//--------- Menu ----------
	TGMenuBar *fMB = new TGMenuBar(this, X,Y);
	AddFrame(fMB, fLH_TX);

    TGPopupMenu *f_Help = fMB->AddPopup("Help");
    f_Help ->AddEntry("spins-reseau.lyx",Id_Reseau_M_documentation_en_lyx);
    f_Help ->AddEntry("spins-reseau.pdf",Id_Reseau_M_documentation_en_pdf);
    f_Help -> Associate(this);




	// zone commune superieure (ZC) ===================
		
		
	auto  * f_ZC = new TGVerticalFrame(this);
	AddFrame(f_ZC,  fLH_L);

    auto  * f_ZC_0 = new TGHorizontalFrame(f_ZC);
	f_ZC->AddFrame(f_ZC_0,  fLH_L);
	
	// zone tab ZT  ================


	auto *fZT = new TGTab(this);
	AddFrame(fZT);

	auto  *f_ZT_Spectre = fZT->AddTab("Spectre");
    auto  *f_ZT_Spectre_0 = new TGHorizontalFrame(f_ZT_Spectre, X, Y);
    f_ZT_Spectre->AddFrame(f_ZT_Spectre_0, fLH_L);
	auto  *f_ZT_Curvature = fZT->AddTab("Curvature");
    auto  *f_ZT_Curvature_0 = new TGHorizontalFrame(f_ZT_Curvature, X, Y);
    f_ZT_Curvature->AddFrame(f_ZT_Curvature_0, fLH_L);
	auto  *f_ZT_Dynamics = fZT->AddTab("Dynamics");
    auto  *f_ZT_Dynamics_0 = new TGHorizontalFrame(f_ZT_Dynamics, X, Y);
    f_ZT_Dynamics->AddFrame(f_ZT_Dynamics_0, fLH_L);

	//-- d'apres la ligne de la classe Reseau: 
    // int N = 4;  //  make_gui =   N(ZC, ":N nombre de sites") help = "conseil N<=10"
	Reseau_N = new TGNumberEntry(f_ZC_0, p_Reseau->N, 3, Id_Reseau_N_N, (TGNumberFormat::EStyle) 0); 
    Reseau_N->GetNumberEntry()->SetToolTipText("conseil N<=10");
    Reseau_N->Associate(this);
    f_ZC_0->AddFrame(Reseau_N,  fLH_L);

	auto *N_lab_Reseau_N = new  TGLabel(f_ZC_0 , ":N nombre de sites");
	f_ZC_0->AddFrame(N_lab_Reseau_N, fLH_L);

	//-- d'apres la ligne de la classe Reseau: 
    // int verbose =0;  //  make_gui =   N(ZC, ":verbose")
	Reseau_verbose = new TGNumberEntry(f_ZC_0, p_Reseau->verbose, 3, Id_Reseau_N_verbose, (TGNumberFormat::EStyle) 0); 
    Reseau_verbose->Associate(this);
    f_ZC_0->AddFrame(Reseau_verbose,  fLH_L);

	auto *N_lab_Reseau_verbose = new  TGLabel(f_ZC_0 , ":verbose");
	f_ZC_0->AddFrame(N_lab_Reseau_verbose, fLH_L);

	//-- d'apres la ligne de la classe Reseau: 
    // int opt_sp = 0; //  make_gui = nl  C(ZT("Spectre"), ": utilise matrice sparse")
    auto  *f_ZT_Spectre_1 = new TGHorizontalFrame(f_ZT_Spectre, X, Y);
    f_ZT_Spectre->AddFrame(f_ZT_Spectre_1, fLH_L);
	Reseau_opt_sp = new TGCheckButton(f_ZT_Spectre_1,  ": utilise matrice sparse",  Id_Reseau_C_opt_sp);
    Reseau_opt_sp->Resize(120, Reseau_opt_sp->GetDefaultHeight());
    if(p_Reseau->opt_sp)
       Reseau_opt_sp->SetState(kButtonDown);
    else
       Reseau_opt_sp->SetState(kButtonUp);
    Reseau_opt_sp->Associate(this);
    f_ZT_Spectre_1->AddFrame(Reseau_opt_sp);

	//-- d'apres la ligne de la classe Reseau: 
    // int Nval= 5*N+1; //  make_gui =   N(ZT("Spectre"), ":(si sparse) nbre de valeurs propres demandees")
	Reseau_Nval = new TGNumberEntry(f_ZT_Spectre_1, p_Reseau->Nval, 7, Id_Reseau_N_Nval, (TGNumberFormat::EStyle) 0); 
    Reseau_Nval->Associate(this);
    f_ZT_Spectre_1->AddFrame(Reseau_Nval,  fLH_L);

	auto *N_lab_Reseau_Nval = new  TGLabel(f_ZT_Spectre_1 , ":(si sparse) nbre de valeurs propres demandees");
	f_ZT_Spectre_1->AddFrame(N_lab_Reseau_Nval, fLH_L);

	//-- d'apres la ligne de la classe Reseau: 
    // double J1=0;  // make_gui = nl  N(ZT("Spectre"), ":J1")
    auto  *f_ZT_Spectre_2 = new TGHorizontalFrame(f_ZT_Spectre, X, Y);
    f_ZT_Spectre->AddFrame(f_ZT_Spectre_2, fLH_L);
	Reseau_J1 = new TGNumberEntry(f_ZT_Spectre_2, p_Reseau->J1, 3, Id_Reseau_N_J1, (TGNumberFormat::EStyle) 5); 
    Reseau_J1->Associate(this);
    f_ZT_Spectre_2->AddFrame(Reseau_J1,  fLH_L);

	auto *N_lab_Reseau_J1 = new  TGLabel(f_ZT_Spectre_2 , ":J1");
	f_ZT_Spectre_2->AddFrame(N_lab_Reseau_J1, fLH_L);

	//-- d'apres la ligne de la classe Reseau: 
    // double J2=1;   // make_gui =   N(ZT("Spectre"), ":J2")
	Reseau_J2 = new TGNumberEntry(f_ZT_Spectre_2, p_Reseau->J2, 3, Id_Reseau_N_J2, (TGNumberFormat::EStyle) 5); 
    Reseau_J2->Associate(this);
    f_ZT_Spectre_2->AddFrame(Reseau_J2,  fLH_L);

	auto *N_lab_Reseau_J2 = new  TGLabel(f_ZT_Spectre_2 , ":J2");
	f_ZT_Spectre_2->AddFrame(N_lab_Reseau_J2, fLH_L);

	//-- d'apres la ligne de la classe Reseau: 
    // void Dessin_spectre(); // make_gui =   B(ZT("Spectre"), "Dessin du spectre") help = "Dessin du spectre. Donner J1,J2"
	Reseau_Dessin_spectre = new TGTextButton(f_ZT_Spectre_2, "Dessin du spectre", Id_Reseau_B_Dessin_spectre);
    Reseau_Dessin_spectre->SetToolTipText("Dessin du spectre. Donner J1,J2");
    Reseau_Dessin_spectre->Associate(this);
    f_ZT_Spectre_2->AddFrame(Reseau_Dessin_spectre,  fLH_L);

	//-- d'apres la ligne de la classe Reseau: 
    // int N1=2;  // make_gui = nl  N(ZT("Spectre"), ":N1")
    auto  *f_ZT_Spectre_3 = new TGHorizontalFrame(f_ZT_Spectre, X, Y);
    f_ZT_Spectre->AddFrame(f_ZT_Spectre_3, fLH_L);
	Reseau_N1 = new TGNumberEntry(f_ZT_Spectre_3, p_Reseau->N1, 3, Id_Reseau_N_N1, (TGNumberFormat::EStyle) 0); 
    Reseau_N1->Associate(this);
    f_ZT_Spectre_3->AddFrame(Reseau_N1,  fLH_L);

	auto *N_lab_Reseau_N1 = new  TGLabel(f_ZT_Spectre_3 , ":N1");
	f_ZT_Spectre_3->AddFrame(N_lab_Reseau_N1, fLH_L);

	//-- d'apres la ligne de la classe Reseau: 
    // int N2=10;  // make_gui =   N(ZT("Spectre"), ":N1")
	Reseau_N2 = new TGNumberEntry(f_ZT_Spectre_3, p_Reseau->N2, 4, Id_Reseau_N_N2, (TGNumberFormat::EStyle) 0); 
    Reseau_N2->Associate(this);
    f_ZT_Spectre_3->AddFrame(Reseau_N2,  fLH_L);

	auto *N_lab_Reseau_N2 = new  TGLabel(f_ZT_Spectre_3 , ":N1");
	f_ZT_Spectre_3->AddFrame(N_lab_Reseau_N2, fLH_L);

	//-- d'apres la ligne de la classe Reseau: 
    // void Dessin_ecart();  // make_gui =   B(ZT("Spectre"), "Dessin des ecarts") help = "Dessin des ecarts. Donner N1,N2"
	Reseau_Dessin_ecart = new TGTextButton(f_ZT_Spectre_3, "Dessin des ecarts", Id_Reseau_B_Dessin_ecart);
    Reseau_Dessin_ecart->SetToolTipText("Dessin des ecarts. Donner N1,N2");
    Reseau_Dessin_ecart->Associate(this);
    f_ZT_Spectre_3->AddFrame(Reseau_Dessin_ecart,  fLH_L);

	//-- d'apres la ligne de la classe Reseau: 
    // int SIZE = 2;  //  make_gui =   N(ZT("Curvature"), ": SIZE") help = "conseil SIZE>=2"
	Reseau_SIZE = new TGNumberEntry(f_ZT_Curvature_0, p_Reseau->SIZE, 3, Id_Reseau_N_SIZE, (TGNumberFormat::EStyle) 0); 
    Reseau_SIZE->GetNumberEntry()->SetToolTipText("conseil SIZE>=2");
    Reseau_SIZE->Associate(this);
    f_ZT_Curvature_0->AddFrame(Reseau_SIZE,  fLH_L);

	auto *N_lab_Reseau_SIZE = new  TGLabel(f_ZT_Curvature_0 , ": SIZE");
	f_ZT_Curvature_0->AddFrame(N_lab_Reseau_SIZE, fLH_L);

	//-- d'apres la ligne de la classe Reseau: 
    // void Curvature_testing();  // make_gui =   B(ZT("Curvature"), "Curvature testing") help = "test. Donner ..."
	Reseau_Curvature_testing = new TGTextButton(f_ZT_Curvature_0, "Curvature testing", Id_Reseau_B_Curvature_testing);
    Reseau_Curvature_testing->SetToolTipText("test. Donner ...");
    Reseau_Curvature_testing->Associate(this);
    f_ZT_Curvature_0->AddFrame(Reseau_Curvature_testing,  fLH_L);

	//-- d'apres la ligne de la classe Reseau: 
    // void Dynamics_testing(); // make_gui = B(ZT("Dynamics"), "Dynamics testing") help "Test."
	Reseau_Dynamics_testing = new TGTextButton(f_ZT_Dynamics_0, "Dynamics testing", Id_Reseau_B_Dynamics_testing);
    Reseau_Dynamics_testing->SetToolTipText("Test.");
    Reseau_Dynamics_testing->Associate(this);
    f_ZT_Dynamics_0->AddFrame(Reseau_Dynamics_testing,  fLH_L);

	// fenetre  generale  ===============
	
	MapSubwindows();
	Resize();
	MapWindow();
	SetWindowName("Commandes de main.cc"); // nom de la fenetre

    //-- debloque le mutex
	mtx.unlock(); 


}

//=====================================
Bool_t Com::ProcessMessage(Long_t msg, Long_t p1, Long_t p2)
{

	int M = GET_MSG(msg), S=GET_SUBMSG(msg);
	//cout<<"Process_Message:  M = "<<M<<"  S = "<<S<<"  p1="<<p1<<"  p2="<<p2<<endl;
//-- M: Menu
    if(M==1 && S==1 && p1== Id_Reseau_M_documentation_en_lyx )
    {
         p_Reseau->documentation_en_lyx();
         p_com->Met_a_jour();
    }

//-- M: Menu
    if(M==1 && S==1 && p1== Id_Reseau_M_documentation_en_pdf )
    {
         p_Reseau->documentation_en_pdf();
         p_com->Met_a_jour();
    }

//-- N: Numerica entry
    if(M==4 && S==1 && p1== Id_Reseau_N_N )
         p_Reseau->N = Reseau_N->GetNumber();

//-- N: Numerica entry
    if(M==4 && S==1 && p1== Id_Reseau_N_verbose )
         p_Reseau->verbose = Reseau_verbose->GetNumber();

//-- C: Check Button
    if(M==1 && S==4 && p1== Id_Reseau_C_opt_sp )
         if(Reseau_opt_sp->IsDown())
            p_Reseau->opt_sp = 1;
         else
            p_Reseau->opt_sp = 0;

//-- N: Numerica entry
    if(M==4 && S==1 && p1== Id_Reseau_N_Nval )
         p_Reseau->Nval = Reseau_Nval->GetNumber();

//-- N: Numerica entry
    if(M==4 && S==1 && p1== Id_Reseau_N_J1 )
         p_Reseau->J1 = Reseau_J1->GetNumber();

//-- N: Numerica entry
    if(M==4 && S==1 && p1== Id_Reseau_N_J2 )
         p_Reseau->J2 = Reseau_J2->GetNumber();

//-- B: Bouton
    if(M==1 && S==3 && p1== Id_Reseau_B_Dessin_spectre )
    {
         p_Reseau->Dessin_spectre();
         p_com->Met_a_jour();
    }

//-- N: Numerica entry
    if(M==4 && S==1 && p1== Id_Reseau_N_N1 )
         p_Reseau->N1 = Reseau_N1->GetNumber();

//-- N: Numerica entry
    if(M==4 && S==1 && p1== Id_Reseau_N_N2 )
         p_Reseau->N2 = Reseau_N2->GetNumber();

//-- B: Bouton
    if(M==1 && S==3 && p1== Id_Reseau_B_Dessin_ecart )
    {
         p_Reseau->Dessin_ecart();
         p_com->Met_a_jour();
    }

//-- N: Numerica entry
    if(M==4 && S==1 && p1== Id_Reseau_N_SIZE )
         p_Reseau->SIZE = Reseau_SIZE->GetNumber();

//-- B: Bouton
    if(M==1 && S==3 && p1== Id_Reseau_B_Curvature_testing )
    {
         p_Reseau->Curvature_testing();
         p_com->Met_a_jour();
    }

//-- B: Bouton
    if(M==1 && S==3 && p1== Id_Reseau_B_Dynamics_testing )
    {
         p_Reseau->Dynamics_testing();
         p_com->Met_a_jour();
    }

return kTRUE;
}


// =========================================
//  Met a jour toutes les valeurs de la fenetre de Commandes,
//   a partir des donnees des classes
void  Com::Met_a_jour()
{
	mtx.lock();
	Reseau_N->SetNumber(p_Reseau->N);
	Reseau_verbose->SetNumber(p_Reseau->verbose);
    if(p_Reseau->opt_sp)
       Reseau_opt_sp->SetState(kButtonDown);
    else
       Reseau_opt_sp->SetState(kButtonUp);
	Reseau_Nval->SetNumber(p_Reseau->Nval);
	Reseau_J1->SetNumber(p_Reseau->J1);
	Reseau_J2->SetNumber(p_Reseau->J2);
	Reseau_N1->SetNumber(p_Reseau->N1);
	Reseau_N2->SetNumber(p_Reseau->N2);
	Reseau_SIZE->SetNumber(p_Reseau->SIZE);
	mtx.unlock();
}
