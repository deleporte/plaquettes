// -*- mode:C++ ; compile-command: "g++ -o makef makef.cc -lm -std=c++11" -*- 

/* 1997. version 2017.  Programme pour creer automatiquement le fichier makefile pour compiler du c++
   ensuite
   ce programme sert a creer un fichier makefile automatiquement.
   en partant d'un fichier "nomfichier.cc" qui contient le main().

   Pour l'utiliser, il faut
   ---------------
   1) ecrire le fichier makef.config

   2) lancer : makef repertoire nomfichier.cc  makef.config    <nom_executable>

   3) ensuite pour compiler (lancer le Makefile), il
   il faut faire :
   make clean   (pour detruire les fichiers .o precedants)
   make all


   voir par exemple le script "compp"



//-----------------------
Compiler avec: make all





//----------------------


Fonctionnement final des commandes:
-----------------
 - Une fenetre de commande gui est lancee dans un autre thread (process) "com" et utilise root.
 - Cette fenetre de commande contient tous les widgets des variables (demandees) des classes que utilise le projet main().
 
 - Communication: 
      - Si un widget est activé (changé), le process "com" change les variables  du main avec une protection prealable de type mutex.
	  
	  - depuis le process "main", on appelle  explicitement p_com->Met_a_jour() pour que l'affichage de "com" corresponde aux variables de  "main".

	  - optionnellement, le process "com" appelle Met_a_jour() de facon periodique (tous les 0.1 sec), pour eviter cet appel explicite.

	  - dans le programme "main", avant/apres d'acceder aux variables partagees, il faut mettre/ enlever le mutex mtx: 	  mtx.lock(); 	  mtx.unlock(); 

	  - si besoin, on a acces aux fenetres de commandes par le pointeur Com *p_com;

sortie
------
  fichiers c++: com.cc, com.h dans le meme repertoire que fichier.h  qui sont des classes c++ d'un  panneau de commandes,
    avec 
	- un menu, 
	- une zone commune: ZC
	- des zones de tabs: ZT1, ZT2, ...

	Les widgets sont placés à la suite horizontalement dans la zone demandee.
	On peut revenir à la ligne si on commence par l'instruction "nl":
	ex:  double x1=4; // make_gui = nl  N(ZT("tab1"), "x1=")


chaque zone peut contenir des widgets avec le format suivant:

     - N: Numerique entre/sortie  pour int, double
	    ex: 
		   int a; //  make_gui = N(ZC,"a=")
		   double x1=4; // make_gui = N(ZT("tab1"), "x1=")

	 - T : Texte entre/sortie      pour string
	    ex:
         string text = "hello" ; //make_gui = T(ZT("tab1"))

	 - HS: Horizontal Slider entre/sortie pour int/double
	    ex: 
		double x = 1 ; // make_gui = HS(ZT("tab1"), 0, 5)

	 - VS: Vertical Slider:   entre/sortie pour int/double

	 - PB: Progress Bar: sortie pour int,double
	   ex:
	    double y = 5 ; // make_gui = PB(ZT("tab2"), 0,5)


	 - C : Check entree/sortie pour int x=0,1
	   ex: int opt = 1 ;  //make_gui = C(ZT("tab2"), ":opt")

	 - L : Liste entre/sortie pour int x \in [0,N-1] 
	  ex:
	  int choix = 2 ; //make_gui = L(ZT("tab2"), {"a","b","c","d"}) 

	 - B : Bouton d'appel pour fonction void f()
	    void f() ; // make_gui = B(ZC,"C1_f")

	 - M : Menu d'appel pour fonction void f()
	    void g() ; // make_gui = M("menu2","C1_g")

projet: 
	 - H1 : Histogramme 1D pour pour vector<int>, vector<double>, vec



les valeurs de ces widgets sont en permanence couplées à la variable c++.
Pour cela il y a un pointeur p_C1 pour chaque classe ex:C1


dans la fonction main(), 
il faut  rajouter (typo):
 #includ "com.h"

Com *com=  new Com(&c1,&c2); 

//thread t(Lance_com, &c1, &c2);  // lance fenetre de commande dans autre thread. Donne pointeurs sur objets
 

rem: cela ne doit pas empecher le fonctionnement de fichier.h sans les commandes
(on doit pouvoir desactiver en option la compilation de com.cc, com.h)
et les fichiers existants ne doivent pas etre modifiés.


//---------------

*/

/*
Faire:


- pouvoir afficher du texte simple: //make_gui = TXT(ZT("tab1"),"texte")

 - pouvoir avoir plusieurs widget pour une variable, avec le symbole & 

- Pour les slider, afficher le nom de la variable (nom)



- dans lorenz, faire que l'on puisse arreter en cours: Checkbouton stop?


- comme c'est souvent utile, rajouter une instruction dans le menu pour lancer  "lyx ./doc.lyx" ou "evince ./doc.pdf", cad lancer une commande bash.



 */


int verbose = 1;

int opt_makeh = 0; //0 ,1  : si ecriture automatique des fichiers .h qui manquent
int opt_makegui = 0; //0 ,1  : si ecriture automatique d'une fenetre de commmandes gui en c++

//*****************************************************************
#include <iostream>
using namespace std;

#include <fstream>
#include <string>
#include <stdlib.h> // pour system
 

#include <vector>
#include <algorithm>
//#include "fichiers/fichiers.h"

#include <vector>
#include <sstream>
// declaration
string traite_cc(string nom_rep_cc, string  nom_fichier_cc,  int opt_rep);

// variables globales
string chaine_objets;

string src2; // ex: essai.cc
string rep_src, rep_src2;
string fliste; // liste des fichiers traites
vector<string> liste_files_h_cc; // liste des fichiers .h .cc du projets
vector<string> liste_files_main; // liste des fichiers du repertoire du main()
string makefile;

string comp_name; // compilateur "g++" ou "emcc"
string ext_o; // extension, sera ".o" ou ".bc" selon le compilateur

int compteur;


/*#################################################
  petite routine qui va pointer apres le prochain caractere '=' */
int prochain(ifstream & fich)
{
	char caract;

	do {
		fich.get(caract);
	}
	while ((caract!='=')&&(!fich.eof()));

	if (fich.eof())
	{
		//     cerr<<" erreur : fin de fichier trouvee..\n";
		return 0;
	}
	else
		return 1;
}

/*#################################################
  petite routine qui va pointer apres la prochaine chaine donnee
  renvoit 0 si pas trouve 1 si trouve  */
int  prochain(ifstream & fich, string chaine)
{// unsigned char caract;
	char caract;
	int n,i;

	n=chaine.length();

//cout<<"recherche dans fichier, la chaine:"<<chaine<<endl;
	do
	{
		// cherche le premier caractere de la chaine
		do {
			fich.get(caract);

			//	   cout<<"."<<caract<<flush;
		}
		while ((caract!=chaine[0])&&(!fich.eof()));


// verifie si il y a le reste de la chaine
		i=0;

		if(!fich.eof())
		{
			do
			{ i++;
				fich.get(caract);
				//cout<<"-"<<caract<<flush;
			}
			while((chaine[i]==caract) && (i<(n-1)));
		}
	}

	while(!(
			  fich.eof()  || (
				  (i==(n-1))  &&   (chaine[n-1]==caract)
				  )
			  ));

	if (fich.eof())
	{
		//     cerr<<" Resultat : fin de fichier trouvee..et chaine non trouvee ..\n";
		return 0;
	}
	else return 1;
}



//===============================
// fonctions pour make_gui

//===============
// enleve les espaces ou tabulation  du debut et fin de chaine
void Enleve_espaces(string &s)
{
	while((s[0] == ' ') || (s[0] == '\t'))
		s.erase(0,1);
	while((s[s.size()-1]==' ') || (s[s.size()-1]=='\t'))
		s.erase(s.size()-1,1);
}
//===================
// renvoit true si l est absent de la liste L
template<typename A>
bool is_not_in(A l, vector<A> L)
{
	return( find(L.begin(), L.end(), l) ==  L.end() );
}


//================
class Infos
{
public:
	string nom_fichier, nom_class, ligne;
	string type, nom, val;
	int opt_nl;
	string widget, zone, texte, nom_tab, xmin, xmax, help;
	string Id, nom_obj;
	vector<string> L_l; // pour Liste
};

//===================

vector<Infos> L_infos;
vector<string> L_noms_classes; // liste des noms de classes
vector<string> L_noms_fichiers; // liste des noms de fichiers contenant des infos
vector<string> L_noms_menus; // liste des noms de menus
vector<string> L_noms_tabs; // liste des noms des zones tab
int n_passage =1; // numero de passage.  il faudra 2 passages si opt_makegui=1.

//==============================
// fonction pour make_gui
/*
Entree:
 p_fich: pointeur dans fichier,

Sortie: 
   0: si rien fait
   1: si ok

Action
-----


Algo:
-----
dans le fichier on cherche le code c++ de la prochaine class

On extrait les infos pour chaque ligne contenant le texte  "make_gui"

 */
int  Extrait_infos_from_next_class(ifstream & p_fich,  string nom_fichier)
{

	
	int res = prochain(p_fich, "class "); // res =0 si pas trouve
	if(res==0)
		return 0;



	if(verbose)
	{
		cout<<"Debut de Fonction Extract_next_class()"<<endl;
		cout<<"----------------------------"<<endl;
	}


	//----- 1) extrait chaque ligne du code c++ de la classe -> L_class ---------------
	vector<string> L_class;

	string t_class; // temporaire pour une ligne
	
	if(verbose)
		cout<<"------"<<endl;

	int cpt=0, cpt_v=0; // compteur des accolades, pour trouver la fin de la classe
	while(p_fich.good())
	{
		char c;
		p_fich.get(c); // -> c
	   
//		cout<<c;
		if(c!= '\n')
			t_class.push_back(c);
		else // nouvelle ligne
		{
			L_class.push_back(t_class);
			t_class.clear();
		}

		//.. bilan des acccolades {, } pour detecter la fin de classe
		if(c=='{')
			cpt++;
		else if(c=='}')
			cpt--;
		if(cpt==0 && cpt_v==1)
			break;
		cpt_v=cpt;

	}

    //-----------------
	if(verbose >=2)
		for(auto s : L_class)
			cout<<s<<endl;

	//--- 2) on extrait les infos pertinentes de chaque ligne, a partir de L_class

	string nom_class = L_class[0];
	Enleve_espaces(nom_class);

		
	if(verbose >=2)
	{
		cout<<"........"<<endl;
		cout<<"Nom  = "<<nom_class<<endl;
		cout<<"........"<<endl;
	}
	int cpt_mk = 0; // nbre de "make_gui" trouves
		
	for(int k= 1; k < L_class.size(); k++) // boucle sur les lignes de code
	{
		string ligne = L_class[k];
		Enleve_espaces(ligne);

		if(verbose >=2)	cout<<"........"<<endl;

		//--- etude de la ligne 
		if(ligne.find("make_gui") != string::npos) // si ligne contient "make_gui"
		{
			cpt_mk++;
			
			string widget, zone, texte, nom_tab, xmin, xmax, help;
			vector<string> L_l; // pour Liste
			
			//ligne0: code c++ avant le signe //  ------------------
			string ligne0 =  ligne.substr(0, ligne.find("//"));  // extrait  avant //

			ligne0 =  ligne0.substr(0, ligne0.find(";"));  // extrait  avant ;

			Enleve_espaces(ligne0);// enleve espaces de debut et fin

			string type =  ligne0.substr(0, ligne0.find(" "));  // extrait  avant " "
			Enleve_espaces(type);

			
			ligne0 =  ligne0.substr(ligne0.find(" ")+1); // extrait apres " "
			
			string nom, val;
			if(type == "void") // c'est une fonction
				nom =  ligne0.substr(0, ligne0.find("("));  // extrait  avant ()
			else
			{
					if(ligne0.find("=") != string::npos) // si contient =
					{
						nom =  ligne0.substr(0, ligne0.find("="));  // extrait  avant =
						val =  ligne0.substr(ligne0.find("=")+1); // extrait apres =
					}
					else
					{
						nom =  ligne0;
					}
			}
			Enleve_espaces(nom);
			Enleve_espaces(val);
			

//			istringstream is(ligne0);
			// string type; is>>type;
			// string nom; is>>nom;

			// if(type=="void") // cas d'une fonction funct()
			// 	nom =  nom.substr(0, nom.find("()"));  // extrait  avant ()
				
			// string val;
			// if(ligne0.find("=") != string::npos) // si contient =
			// {
			// 	is>>val; is >>val;
			// }
			
			if(verbose >=2)
			{
				cout<<"ligne0 = "<<ligne0<<endl;
				cout<<"   type  = "<<type<<endl;
				cout<<"   nom  = "<<nom<<endl;
				cout<<"   val  = "<<val<<endl<<endl;
			}

			//ligne1 -----------------------
		
			
			string ligne1 =  ligne.substr(ligne.find("//"));  // extrait apres //
			ligne1 =  ligne1.substr(ligne1.find("=")+1); // extrait apres '='
			Enleve_espaces(ligne1);// enleve espaces de debut et fin
			
			if(verbose >=2)
				cout<<"ligne1 = "<<ligne1<<endl;

			//.....
			int opt_nl = 0;
			if(ligne1[0] == 'n' && ligne1[1] == 'l') // si "nl"
			{
				opt_nl = 1; // demande de newline
				ligne1.erase(0,2); // enleve "nl" de debut
			}
			
			//.....
			widget = ligne1.substr(0, ligne1.find("(")); // extrait avant (
			Enleve_espaces(widget);// enleve espaces de debut et fin
			if(verbose >=2)
				cout<<"widget="<<widget<<endl;
			

			if(ligne.find("help") != string::npos) // si ligne contient "help"
			{
				help = ligne1.substr(ligne1.find("help"));  // extrait apres "help"
				help = help.substr(help.find("\"")+1);  // extrait apres "
				help = help.substr(0,help.rfind("\""));  // et avant " de la fin
				
				ligne1 = ligne1.substr(0, ligne1.find("help")); // extrait avant "help"
			}

			if(verbose >=2)
				cout<<"ligne1 = "<<ligne1<<endl;
				
			ligne1 = ligne1.substr(ligne1.find("(")+1);  // extrait apres (
			ligne1 = ligne1.substr(0,ligne1.rfind(")"));  // et avant ) de la fin
			if(verbose >=2)
				cout<<"ligne1 = "<<ligne1<<endl;

			//.. lecture de zone 
			if(widget == "M") // si menu
			{
				nom_tab = ligne1.substr(0, ligne1.find(","));  // extrait avant ,
				nom_tab = nom_tab.substr(nom_tab.find("\"")+1);  // extrait apres "
				nom_tab = nom_tab.substr(0, nom_tab.find("\"")); // extrait avant "
				if(is_not_in(nom_tab,  L_noms_menus)) // si nom_tab  pas deja present dans liste L_noms_menus
					L_noms_menus.push_back(nom_tab);
				
			}
			else // sinon
			{
				zone = ligne1.substr(0, ligne1.find(","));  // extrait avant ,
				Enleve_espaces(zone);
				if(zone[1] == 'T') // si ZT
				{
					nom_tab = zone.substr(zone.find("\"")+1);  // extrait apres "
					nom_tab = nom_tab.substr(0, nom_tab.find("\"")); // extrait avant "
					if(is_not_in(nom_tab,  L_noms_tabs))// si pas deja present
						L_noms_tabs.push_back(nom_tab);

					zone = "ZT";
				}
			}
			
			//........ suite
			
			if(widget=="N" || widget=="C" || widget=="B"  || widget=="M"  ) // Numerical or Check or Button or Menu
			{//-- on cherche info: zone, titre
				ligne1 = ligne1.substr(ligne1.find(",")+1);  // extrait apres ,
				Enleve_espaces(ligne1);
				ligne1.erase(0,1); // enleve " de debut
				ligne1.erase(ligne1.size()-1,1); // enleve " de fin
				texte = ligne1;
				
			}
			
			else if(widget=="T") //Text
			{
				
			}

			else if(widget=="HS" || widget=="VS" || widget=="PB" ) // Slider or ProgressBar
			{
				ligne1 = ligne1.substr(ligne1.find(",")+1);  // extrait apres ,
				xmin = ligne1.substr(0, ligne1.find(","));  // extrait avant ,
				Enleve_espaces(xmin);
				xmax = ligne1.substr(ligne1.find(",")+1);  // extrait apres ,
				Enleve_espaces(xmax);
			}
			
			else if(widget=="L") // List
			{
				ligne1 = ligne1.substr(ligne1.find(",")+1);  // extrait apres ,
				while(ligne1.find("\"") != string::npos) // si contient "
				{
					ligne1 = ligne1.substr(ligne1.find("\"")+1);  // extrait apres "
					string l = ligne1.substr(0, ligne1.find("\""));  // extrait avant "
					L_l.push_back(l);
					ligne1 = ligne1.substr(ligne1.find("\"")+1);  // extrait apres "
				
				}
			}

			if(verbose >=2)
			{
				cout<<"   widget  = \""<<widget<<"\""<<endl;
				cout<<"   zone  = \""<<zone<<"\""<<endl;
				cout<<"   nom_tab  = \""<<nom_tab<<"\""<<endl;
				cout<<"   texte  = \""<<texte<<"\""<<endl;
				cout<<"   xmin  = \""<<xmin<<"\""<<endl;
				cout<<"   xmax  = \""<<xmax<<"\""<<endl;
				for(int i=0; i<L_l.size(); i++)
					cout<<"   L_l["<<i<<"]  = \""<<L_l[i]<<"\""<<endl;
				
				cout<<endl;
			}
			
			Infos infos;
			infos.nom_fichier = nom_fichier;
			infos.nom_class = nom_class;
			infos.ligne = ligne;
			infos.type = type;
			infos.nom =nom;
			infos.val = val;
			infos.opt_nl = opt_nl;
			infos.widget = widget;
			infos.zone=zone;
			infos.nom_tab=nom_tab;
			infos.texte=texte;
			infos.xmin=xmin;
			infos.xmax=xmax;
			infos.L_l=L_l;
			infos.help = help;
			
			infos.Id = "Id_" + nom_class + "_" + widget + "_" + nom; //  -> nom Identificateur
			infos.nom_obj = nom_class + "_" +  nom; // -> nom de l'objet c++
			
			L_infos.push_back(infos);
				
		}
	}

    //--------------------------------
	if(cpt_mk>=1)
	{
		L_noms_classes.push_back(nom_class); 
		if(is_not_in(nom_fichier, L_noms_fichiers))
			L_noms_fichiers.push_back(nom_fichier);
	}
	//---------------------
	if(verbose)
	{
		cout<<" Fin de Fonction Extract_next_class()"<<endl;
		cout<<"----------------------------"<<endl<<endl;;
	}
	return 1;
	
}
//=================================
void Bilan_collecte_infos()
{
	
	
	//-------------------------------
	
	cout<<"Bilan de la collecte des infos:"<<endl;
	cout<<endl;

	
	for(auto nc : L_noms_classes)
		cout<<"  Classe: "<<nc<<endl; 

	for(auto nc : L_noms_menus)
		cout<<"  Menu: "<<nc<<endl; 

	for(auto nc : L_noms_tabs)
		cout<<"  Tab: "<<nc<<endl; 
	
	
	cout<<endl;
	cout<<" L_infos.size() = "<<L_infos.size()<<endl;

	for(auto inf : L_infos)
	{
		cout<<"---------------"<<endl;
		cout<<"   nom_fichier  = \""<<inf.nom_fichier<<"\""<<endl;
		cout<<"   nom_class  = \""<<inf.nom_class<<"\""<<endl;
		cout<<"   ligne  = \""<<inf.ligne<<"\""<<endl;
		cout<<"   ..............."<<endl;
		cout<<"   type  = \""<<inf.type<<"\""<<endl;
		cout<<"   nom  = \""<<inf.nom<<"\""<<endl;
		cout<<"   val  = \""<<inf.val<<"\""<<endl;
		cout<<"   ..............."<<endl;
		cout<<"   widget  = \""<<inf.widget<<"\""<<endl;
		cout<<"   zone  = \""<<inf.zone<<"\""<<endl;
		cout<<"   texte  = \""<<inf.texte<<"\""<<endl;
		cout<<"   nom_tab  = \""<<inf.nom_tab<<"\""<<endl;
		cout<<"   xmin  = \""<<inf.xmin<<"\""<<endl;
		cout<<"   xmax  = \""<<inf.xmax<<"\""<<endl;
	
		for(int i=0; i<inf.L_l.size(); i++)
			cout<<"   L_l["<<i<<"]  = \""<<inf.L_l[i]<<"\""<<endl;
		cout<<"   help  = \""<<inf.help<<"\""<<endl;
		
		cout<<"   ..............."<<endl;
		cout<<endl;
	}



	
}

//============================
// fonction pour make_gui
/*
  a partir des informations contenues dans
  L_infos
*/
void 	Ecriture_fichiers_cpp()
{
	
	cout<<"##== Part2: Ecriture des fichiers com.cc com.h  =========="<<endl;
	cout<<"=========================="<<endl<<endl;
	cout<<"Debut de Fonction Ecriture_fichiers_cpp()"<<endl;
	cout<<"----------------------------"<<endl;

	string nom_com = rep_src2 + "/com.cc";
	
	//---fichier com.cc ---------------------------------
	cout<<"Ecriture du fichier :"<<nom_com<<endl;
	ofstream f(nom_com);

	f<<"//CE CODE C++ com.cc EST ECRIT AUTOMATIQUEMENT PAR: makef.cc"<<endl;
	f<<"//========================================================"<<endl;

	//---- ecriture des entetes
	f<<" \n\
# include \"./com.h\" \n\
";



	f<<endl<<endl<<"//=================="<<endl;

	for(auto nc : L_noms_classes)
		f<<nc<<" *p_"<<nc<<"; // pointeur sur objet de la classe "<<nc<<endl;

	f<<"mutex mtx; \n\
Com *p_com; \n\
\n\
//=================== \n\
// Lance la fenetre de commandes \n\
// \n\
void Lance_com(";

	for(int i=0; i<L_noms_classes.size(); i++)
	{
		auto nc = L_noms_classes[i];
		f<<nc<<" *p"<<nc;
		if(i<L_noms_classes.size()-1)
			f<<", ";
	}
	f<<")"<<endl;

		
	f<<"{"<<endl<<"\n\
 	sleep_for(milliseconds{100}); // attente pour etre sur que un autre root a le temps de se lancer si besoin. \n \
	TApplication theApp(\"App2\", nullptr, nullptr); \n\
	p_com=  new Com(";

	for(int i=0; i<L_noms_classes.size(); i++)
	{
		auto nc = L_noms_classes[i];
		f<<"p"<<nc;
		if(i<L_noms_classes.size()-1)
			f<<", ";
	}
	f<<");  // cree fenetre de control et l'associe a  l'objet present"<<endl;

	f<<" 	theApp.Run(); //----donne la main a l'utilisateur --- \n\
}"<<endl;
	
	f<<endl;
	f<<"//====== identifieurs (int) pour les messages de commandes\n\
enum Command_id \n\
{"<<endl;


	for(int j=0;j<L_infos.size(); j++)
	{
		f<<L_infos[j].Id;
		if(j<(L_infos.size()-1))
			f<<", "<<endl;
	}
	f<<endl<<"}; \n\
\n\
int X=300, Y=150; // taille en pixels\n\
\n\
//===================\n\
Com::Com(";

	for(int i=0; i<L_noms_classes.size(); i++)
	{
		auto nc = L_noms_classes[i];
		f<<nc<<" *p"<<nc;
		if(i<L_noms_classes.size()-1)
			f<<", ";
	}
	f<<")"<<endl;	

	f<<" : TGMainFrame(gClient->GetRoot(), X,Y) \n\
{\n\
\n\
\n\
	//-- bloque le mutex\n\
	mtx.lock(); // si mtx est libre on le bloque, sinon on attend\n\
\n\
\n\
	//---- pointeurs sur objets\n\
";
	for(auto nc : L_noms_classes)
		f<<"     p_"<<nc<<" = p"<<nc<<";"<<endl;


	f<<"\n\
\n\
	//-- regles de positionnement\n\
	TGLayoutHints *fLH_TX =  new TGLayoutHints(kLHintsTop | kLHintsExpandX , 2, 2, 2, 5);\n\
	TGLayoutHints *fLH_C = new TGLayoutHints(kLHintsCenterX | kLHintsExpandX | kLHintsExpandY,5,5,5,5);\n\
	TGLayoutHints *fLH_L = new TGLayoutHints(kLHintsLeft,5,5,5,5);\n\
	TGLayoutHints *fLH_R = new TGLayoutHints(kLHintsRight,5,5,5,5);\n\
\n\
\n\
\n";


	if(L_noms_menus.size()>0)
	{
		f<<"//--------- Menu ----------\n\
	TGMenuBar *fMB = new TGMenuBar(this, X,Y);\n\
	AddFrame(fMB, fLH_TX);\n\
\n\
";

		for(auto menu : L_noms_menus)
		{
			f<<"    TGPopupMenu *f_"<<menu<<" = fMB->AddPopup(\""<<menu<<"\");"<<endl;
			//... boucle sur les objets a mettre au menu
			for(auto inf : L_infos) // boucle sur les objets
				if(inf.widget == "M" && inf.nom_tab == menu)
					f<<"    f_"<<menu<<" ->AddEntry(\""<<inf.texte<<"\","<<inf.Id<<");"<<endl;
		
			f<<"    f_"<<menu<<" -> Associate(this);"<<endl;
			f<<endl;
		}
	}


	int cpt_n_ZC = 0; // numero de la frame en cours de zone commune.
	
	
	f<<"\n\
\n\
\n\
	// zone commune superieure (ZC) ===================\n\
		\n\
		\n\
	auto  * f_ZC = new TGVerticalFrame(this);\n\
	AddFrame(f_ZC,  fLH_L);\n\
\n\
    auto  * f_ZC_"<<cpt_n_ZC<<" = new TGHorizontalFrame(f_ZC);\n\
	f_ZC->AddFrame(f_ZC_"<<cpt_n_ZC<<",  fLH_L);\n\
	"<<endl;

	if(L_noms_tabs.size()>0)
	{
		f<<"	// zone tab ZT  ================\n\
\n\
\n\
	auto *fZT = new TGTab(this);\n\
	AddFrame(fZT);\n\
"<<endl;
	}

	vector<int> cpt_n_ZT;
	
	for(auto nom_tab : L_noms_tabs)
	{
		string var_tab0 = "f_ZT_" + nom_tab;
		string var_tab = var_tab0 + "_0";
		cpt_n_ZT.push_back(0);
		
		f<<"	auto  *"<<var_tab0<<" = fZT->AddTab(\""<<nom_tab<<"\");"<<endl;
		f<<"    auto  *"<<var_tab<<" = new TGHorizontalFrame("<<var_tab0<<", X, Y);"<<endl;
		f<<"    "<<var_tab0<<"->AddFrame("<<var_tab<<", fLH_L);"<<endl;
	}

	
	f<<endl;
		
	//----boucle sur les objets a mettre dans les zones  ZC ou ZT
	for(auto inf : L_infos) // boucle sur les objets
	{
		if(inf.widget == "M") // menu
			continue;
		
		string var_tab, var_tab0;
		if(inf.zone == "ZT")
		{
			int pos = find(L_noms_tabs.begin(), L_noms_tabs.end(), inf.nom_tab) - L_noms_tabs.begin();
			if(inf.opt_nl == 1) // newline
				cpt_n_ZT[pos] = cpt_n_ZT[pos]+1;
			
			var_tab = "f_" + inf.zone + "_" + inf.nom_tab + "_" + to_string(cpt_n_ZT[pos]);
			var_tab0 = "f_ZT_" +  inf.nom_tab; // parent
			
			
		}
		else if(inf.zone == "ZC")
		{
			var_tab0 = "f_ZC"; // parent
			if(inf.opt_nl == 1) // newline
				cpt_n_ZC++;
			
			var_tab = "f_ZC_" + to_string(cpt_n_ZC);
		}
	
			
	
		f<<"	//-- d'apres la ligne de la classe "<<inf.nom_class<<": "<<endl;
		f<<"    // "<<inf.ligne<<endl;

		if(inf.opt_nl == 1) // newline
		{
			f<<"    auto  *"<<var_tab<<" = new TGHorizontalFrame("<<var_tab0<<", X, Y);"<<endl;
			f<<"    "<<var_tab0<<"->AddFrame("<<var_tab<<", fLH_L);"<<endl;
		}
		
		//----------------------------
		if(inf.widget == "N")
		{

			
			string style = "(TGNumberFormat::EStyle) 0"; // pour type = int
			if(inf.type=="double")
				style = "(TGNumberFormat::EStyle) 5"; // cf ligne 33 de https://root.cern.ch/doc/master/TGNumberEntry_8h_source.html

			int size = 2 + (inf.val).size();
			
			f<<"	"<<inf.nom_obj<<" = new TGNumberEntry("<<var_tab<<", p_"<<inf.nom_class<<"->"<<inf.nom<<", "<<size<<", "<<inf.Id<<", "<<style<<"); "<<endl;
			if((inf.help).size()>0)
				f<<"    "<<inf.nom_obj<<"->GetNumberEntry()->SetToolTipText(\""<< inf.help <<"\");"<<endl;
			f<<"    "<<inf.nom_obj<<"->Associate(this);"<<endl;
			f<<"    "<<var_tab<<"->AddFrame("<<inf.nom_obj<<",  fLH_L);"<<endl;
			f<<endl;

			//.... texte 
			string var2 = "N_lab_" + inf.nom_class + "_" + inf.nom;
			f<<"	auto *" << var2 << " = new  TGLabel("<<var_tab<<" , \""<<inf.texte<<"\");"<<endl;
			f<<"	"<<var_tab<<"->AddFrame("<<var2<<", fLH_L);"<<endl;
			f<<endl;

				
		}
		//----------------------------
		else 	if(inf.widget == "B")
		{
		
			f<<"	" << inf.nom_obj << " = new TGTextButton("<<var_tab<<", \""<<inf.texte<<"\", "<<inf.Id<<");"<<endl;
			if((inf.help).size()>0)
				f<<"    "<<inf.nom_obj<<"->SetToolTipText(\""<< inf.help <<"\");"<<endl;

			f<<"    "<<inf.nom_obj<<"->Associate(this);"<<endl;
			f<<"    "<<var_tab<<"->AddFrame("<<inf.nom_obj<<",  fLH_L);"<<endl;
			f<<endl;
				
		}
		//----------------------------
		else 	if(inf.widget == "T")
		{
		
			f<<"	" << inf.nom_obj << " = new TGTextEntry("<<var_tab<<",  p_"<<inf.nom_class<<"->"<<inf.nom<<".c_str(), "<<inf.Id<<");"<<endl;
			f<<"    "<<inf.nom_obj<<"->Associate(this);"<<endl;
			f<<"    "<<var_tab<<"->AddFrame("<<inf.nom_obj<<");"<<endl;
			f<<endl;
				
		}
		//----------------------------
		else 	if(inf.widget == "HS")
		{
		
			f<<"	" << inf.nom_obj << " = new TGHSlider("<<var_tab<<", 150, kSlider1|kScaleDownRight, "<<inf.Id<<");"<<endl;
			f<<"    "<<inf.nom_obj<<"->SetRange("<<inf.xmin<<", "<<inf.xmax<<");"<<endl;
			f<<"    "<<inf.nom_obj<<"->SetPosition(p_"<<inf.nom_class<<"->"<<inf.nom<<");"<<endl;
			f<<"    "<<inf.nom_obj<<"->Associate(this);"<<endl;
			f<<"    "<<var_tab<<"->AddFrame("<<inf.nom_obj<<");"<<endl;
			f<<endl;
				
		}
		//----------------------------
		else 	if(inf.widget == "PB")
		{
	
			f<<"	" << inf.nom_obj << " = new TGHProgressBar("<<var_tab<<", TGProgressBar::kFancy, 200);"<<endl;
			f<<"    "<<inf.nom_obj<<"->ShowPosition(kTRUE,kFALSE,\"y: %.0f\");"<<endl;
			f<<"    "<<inf.nom_obj<<"->SetBarColor(\"green\");"<<endl;
			f<<"    "<<inf.nom_obj<<"->SetRange("<<inf.xmin<<", "<<inf.xmax<<");"<<endl;
			f<<"    "<<inf.nom_obj<<"->SetPosition(p_"<<inf.nom_class<<"->"<<inf.nom<<");"<<endl;
			//		f<<"    "<<inf.nom_obj<<"->Associate(this);"<<endl;
			f<<"    "<<var_tab<<"->AddFrame("<<inf.nom_obj<<");"<<endl;
			f<<endl;
				
		}
		
		//----------------------------
		else 	if(inf.widget == "C")
		{
		
			f<<"	" << inf.nom_obj << " = new TGCheckButton("<<var_tab<<",  \""<<inf.texte<<"\",  "<<inf.Id<<");"<<endl;
			if((inf.help).size()>0)
				f<<"    "<<inf.nom_obj<<"->SetToolTipText(\""<< inf.help <<"\");"<<endl;
			
			f<<"    "<<inf.nom_obj<<"->Resize(120, "<<inf.nom_obj<<"->GetDefaultHeight());"<<endl;
			f<<"    if(p_"<<inf.nom_class<<"->"<<inf.nom<<")"<<endl;
			f<<"       "<<inf.nom_obj<<"->SetState(kButtonDown);"<<endl;
			f<<"    else"<<endl;
			f<<"       "<<inf.nom_obj<<"->SetState(kButtonUp);"<<endl;
			f<<"    "<<inf.nom_obj<<"->Associate(this);"<<endl;
			f<<"    "<<var_tab<<"->AddFrame("<<inf.nom_obj<<");"<<endl;
			f<<endl;
				
		}
			
		//----------------------------
		else 	if(inf.widget == "L")
		{
	
			f<<"	" << inf.nom_obj << " = new TGComboBox("<<var_tab<<", "<<inf.Id<<");"<<endl;
			for(int j=0; j<inf.L_l.size(); j++)
				f<<"    "<<inf.nom_obj<<"->AddEntry(\""<<inf.L_l[j]<<"\","<<j<<");"<<endl;
			f<<"    "<<inf.nom_obj<<"->Resize(150, 20);"<<endl;
			f<<"    "<<inf.nom_obj<<"->Select( p_"<<inf.nom_class<<"->"<<inf.nom<<");"<<endl;
			f<<"    "<<inf.nom_obj<<"->Associate(this);"<<endl;
			f<<"    "<<var_tab<<"->AddFrame("<<inf.nom_obj<<");"<<endl;
			f<<endl;
		}
	}


	f<<"	// fenetre  generale  ===============\n\
	\n\
	MapSubwindows();\n\
	Resize();\n\
	MapWindow();\n\
	SetWindowName(\"Commandes de "<< src2 <<"\"); // nom de la fenetre\n\
\n\
    //-- debloque le mutex\n\
	mtx.unlock(); \n\
\n\
\n\
}"<<endl;

	f<<"\n\
//=====================================\n\
Bool_t Com::ProcessMessage(Long_t msg, Long_t p1, Long_t p2)\n\
{\n\
\n\
	int M = GET_MSG(msg), S=GET_SUBMSG(msg);\n\
	//cout<<\"Process_Message:  M = \"<<M<<\"  S = \"<<S<<\"  p1=\"<<p1<<\"  p2=\"<<p2<<endl;\n\
";

	for(auto inf : L_infos) // boucle sur les objets
	{
		if(inf.widget == "N")
		{
			f<<"//-- N: Numerica entry"<<endl;
			f<<"    if(M==4 && S==1 && p1== "<<inf.Id<<" )"<<endl;
			f<<"         p_"<<inf.nom_class<<"->"<<inf.nom<<" = "<<inf.nom_obj<<"->GetNumber();"<<endl;
			f<<endl;
		}
		else	if(inf.widget == "T")
		{
			f<<"//-- T:  Text entry"<<endl;
			f<<"    if(M==4 && S==1 && p1== "<<inf.Id<<" )"<<endl;
			f<<"         p_"<<inf.nom_class<<"->"<<inf.nom<<" = "<<inf.nom_obj<<"->GetText();"<<endl;
			f<<endl;
		}
		else if(inf.widget == "B")
		{
			f<<"//-- B: Bouton"<<endl;
			f<<"    if(M==1 && S==3 && p1== "<<inf.Id<<" )"<<endl;
			f<<"    {"<<endl;
			f<<"         p_"<<inf.nom_class<<"->"<<inf.nom<<"();"<<endl;
			f<<"         p_com->Met_a_jour();"<<endl;
			f<<"    }"<<endl;
			f<<endl;
		}
		else if(inf.widget == "M")
		{
			f<<"//-- M: Menu"<<endl;
			f<<"    if(M==1 && S==1 && p1== "<<inf.Id<<" )"<<endl;
			f<<"    {"<<endl;
			f<<"         p_"<<inf.nom_class<<"->"<<inf.nom<<"();"<<endl;
			f<<"         p_com->Met_a_jour();"<<endl;
			f<<"    }"<<endl;
			f<<endl;
		}
		else if(inf.widget == "HS")
		{
			f<<"//-- HS: Slider"<<endl;
			f<<"    if(M==6 && S==1 && p1== "<<inf.Id<<" )"<<endl;
			f<<"         p_"<<inf.nom_class<<"->"<<inf.nom<<" = "<<inf.nom_obj<<"->GetPosition();"<<endl;
			f<<endl;
		}
		else if(inf.widget == "C")
		{
			f<<"//-- C: Check Button"<<endl;
			f<<"    if(M==1 && S==4 && p1== "<<inf.Id<<" )"<<endl;
			f<<"         if("<<inf.nom_obj<<"->IsDown())"<<endl;
			f<<"            p_"<<inf.nom_class<<"->"<<inf.nom<<" = 1;"<<endl;
			f<<"         else"<<endl;
			f<<"            p_"<<inf.nom_class<<"->"<<inf.nom<<" = 0;"<<endl;
			f<<endl;
		}
		else if(inf.widget == "L")
			{
			f<<"//-- L: List Combo"<<endl;
			f<<"    if(M==1 && S==7 && p1== "<<inf.Id<<" )"<<endl;
			f<<"         p_"<<inf.nom_class<<"->"<<inf.nom<<" = "<<inf.nom_obj<<"->GetSelected();"<<endl;
			f<<endl;
		}
	}
	
	f<<"return kTRUE;"<<endl;
	f<<"}"<<endl<<endl;

	

	
	f<<"\n\
// =========================================\n\
//  Met a jour toutes les valeurs de la fenetre de Commandes,\n\
//   a partir des donnees des classes\n\
void  Com::Met_a_jour()\n\
{\n\
	mtx.lock();\n\
";
	for(auto inf : L_infos) // boucle sur les objets
	{
		if(inf.widget == "N")
			f<<"	"<<inf.nom_obj<<"->SetNumber(p_"<<inf.nom_class<<"->"<<inf.nom<<");"<<endl;
		else if(inf.widget == "T")
			f<<"	"<<inf.nom_obj<<"->SetText(p_"<<inf.nom_class<<"->"<<inf.nom<<".c_str());"<<endl;
		else if(inf.widget == "HS")
			f<<"    "<<inf.nom_obj<<"->SetPosition(p_"<<inf.nom_class<<"->"<<inf.nom<<");"<<endl;
		else if(inf.widget == "PB")
			f<<"    "<<inf.nom_obj<<"->SetPosition(p_"<<inf.nom_class<<"->"<<inf.nom<<");"<<endl;
		else if(inf.widget == "C")
		{
			f<<"    if(p_"<<inf.nom_class<<"->"<<inf.nom<<")"<<endl;
			f<<"       "<<inf.nom_obj<<"->SetState(kButtonDown);"<<endl;
			f<<"    else"<<endl;
			f<<"       "<<inf.nom_obj<<"->SetState(kButtonUp);"<<endl;
		}
		else if(inf.widget == "L")
			f<<"    "<<inf.nom_obj<<"->Select( p_"<<inf.nom_class<<"->"<<inf.nom<<");"<<endl;
	}
	f<<"	mtx.unlock();"<<endl;
	f<<"}"<<endl;
	f.close(); // ferme le fichier com.cc

	
	//=======================================
	nom_com = rep_src2 + "/com.h";
	
	//---fichier com.h ---------------------------------
	cout<<"Ecriture du fichier :"<<nom_com<<endl;
	f.open(nom_com);
	//-----------------------------------------------------------------
	f<<"\n\
#if !defined(__COM_MODEL_H)\n\
#define __COM_MODEL_H\n\
//CE CODE C++ com.h EST ECRIT AUTOMATIQUEMENT PAR: makef.cc\n\
//============================================================\n\
\n\
# include <TCanvas.h> \n\
# include <TGFrame.h> \n\
# include <TGButton.h> \n\
# include <TGNumberEntry.h> \n\
# include <TGMenu.h> \n\
# include <TGLabel.h> \n\
# include <TApplication.h> \n\
# include <TGTab.h> \n\
# include <TGSlider.h> \n\
# include <TGProgressBar.h> \n\
# include <TGComboBox.h>  \n\
# include <thread> \n\
# include <chrono> \n\
using namespace std::chrono; \n\
using namespace std::this_thread; \n\
# include <iostream> \n\
using namespace std; \n\
# include <mutex> \n";

	for(auto nf : L_noms_fichiers)
		f<<"# include \""<<nf<<"\"\n";		
		


	f<<"\n\
//==================\n\
void Lance_com(";

	for(int j=0; j<L_noms_classes.size(); j++)
	{
		auto nc=L_noms_classes[j];
		f<<nc<<" *p"<<nc;
		if(j<(L_noms_classes.size()-1))
			f<<", ";
		else
			f<<");"<<endl;
	}
	f<<"\n\
//===================\n\
class Com: public  TGMainFrame\n\
{\n\
public :\n\
\n\
	Com(";
	for(int j=0; j<L_noms_classes.size(); j++)
	{
		auto nc=L_noms_classes[j];
		f<<nc<<" *p"<<nc;
		if(j<(L_noms_classes.size()-1))
			f<<", ";
		else
			f<<");"<<endl;
	}

	for(auto inf : L_infos) // boucle sur les objets
	{
		if(inf.widget == "N")
			f<<"	TGNumberEntry *"<<inf.nom_obj<<";"<<endl;
		else if(inf.widget == "T")
			f<<"	TGTextEntry *"<<inf.nom_obj<<";"<<endl;
		else if(inf.widget == "HS")
			f<<"	TGHSlider *"<<inf.nom_obj<<";"<<endl;
		else if(inf.widget == "PB")
			f<<"	TGHProgressBar *"<<inf.nom_obj<<";"<<endl;
		else if(inf.widget == "C")
			f<<"	TGCheckButton *"<<inf.nom_obj<<";"<<endl;
		else if(inf.widget == "L")
			f<<"	TGComboBox *"<<inf.nom_obj<<";"<<endl;
		else if(inf.widget == "B")
			f<<"	TGTextButton *"<<inf.nom_obj<<";"<<endl;
	}
	
	f<<"\n\
	Bool_t ProcessMessage(Long_t msg, Long_t p1, Long_t p2);\n\
	void  Met_a_jour();\n\
};\n\
\n\
//============\n\
\n\
extern Com *p_com;  // pointeur sur l'objet (declare dans .cc)\n\
extern mutex mtx;\n\
\n\
#endif\n\
";









	
//------------------------------	
	cout<<" Fin de Fonction Ecriture_fichiers_cpp()"<<endl;
	cout<<"----------------------------"<<endl<<endl;
}


//===========================================
// fonction pour make_gui
void Extrait_infos_from_file(string nom_rep, string nom_fich)
{
	string nom_rep_fichier = nom_rep + "/" + nom_fich; 
	ifstream p_fich(nom_rep_fichier);
	if(!p_fich.good())
		cout<<"Erreur ouverture du fichier "<<nom_rep_fichier<<" dans Extrait_infos_from_file()"<<endl;

	else
		cout<<" One extrait les infos depuis le fichier: "<<nom_rep_fichier<<":"<<endl<<endl;



	
	int res=1;
	do
		res = Extrait_infos_from_next_class(p_fich, nom_fich);
	while (res!=0);
	
}


//=================================================
/* recherche dans le fichier de nom fliste si le fichier nom_fichier s y trouve  renvoit 1:oui  0:non */

int deja_traite(string & nom_fichier)
{
	if(verbose>=2)
		cout<<"Fonction: deja_traité() pour le fichier: "<<nom_fichier<<endl;

	ifstream pliste(fliste.c_str());

	int res=prochain(pliste,nom_fichier);
	//cout<<"resultat de deja_traite pour fichier:"<< nom_fichier<<" res="<<res<<endl;
	return res;
}

//=================================================
/* rajoute le fichier nom_fichier a la liste de nom fliste*/

void rajoute_liste(string & nom_fichier)
{
	if(verbose>=2)
		cout<<"Rajoute à la liste"<<endl;

	ofstream pliste(fliste.c_str(),ios::app);
	pliste<<compteur<<" : "<<nom_fichier<<"\n ";
}

//===========================================================
/* fonction qui recherche la chaine 
 1):  #include "chemin/toto.h" ou 
 2):  #include "./chemin/toto.h"
dans le fichier pointe par p_fichier_cc (recherche a partir du pointeur).
   si trouve: renvoit 1 (ou 2)   et pointeur sur la chaine "chemin/toto.h"
   sinon renvoit 0.
*/
int recherche_h(ifstream & p_fichier_cc,string & chaine)
{

	
	int res;
	do
    {
		res=prochain(p_fichier_cc,"#include"); //-> res=0 ou 1
		if(res)  // #include trouve
		{
			p_fichier_cc>>chaine;

			
			
			if (chaine[0]=='\"')  // candidat recherche
			{
				chaine.erase(0,1);  // enleve les caracteres " de debut et fin
				chaine.erase(chaine.size()-1,1);

				
				if(chaine[0] == '.') 
					if(chaine[1] == '/') // si "./chemin/toto.h"
					{
						res=2;
						chaine.erase(0,2);  // enleve les caracteres ./ du debut
					}
			}
			else
				res=3; //include non recherche: ex: #include <iostream>
		}
    }
	while (res == 3);
	return res;
}
//==========================================================
/* fonction qui traite le cas d'un fichier toto.h cad operation (B)*/
string traite_h(string nom_rep_h, string  nom_fichier_h, int opt_rep)
{
	
	string text;

	string nom_rep_fichier_h = nom_rep_h + "/" + nom_fichier_h;

	
	if (deja_traite(nom_rep_fichier_h))  // si le fichier a deja ete traite
		return text;
	else
	{
		rajoute_liste(nom_rep_fichier_h); //le rajoute a la liste

		if(opt_rep == 1)
			liste_files_h_cc.push_back(nom_fichier_h);
		else if(opt_rep == 2)
			liste_files_main.push_back(nom_fichier_h);
		
	}

	compteur++;

	if(verbose)
	{
		for (int i=1;i<=compteur;i++) cout<<" ";
		cout<<compteur<<" Traitement du fichier :"<<nom_rep_fichier_h<<"\n";
	}

	{   // debut partie 1
// ouvre le fichier...........................
		ifstream p_fichier_h(nom_rep_fichier_h.c_str());
		if(!p_fichier_h)
			cerr<<"erreur ouverture du fichier :"<<nom_rep_fichier_h<<" dans (***)\n";


// trouve la liste des .cc inclus ... pour (A)
		if(prochain(p_fichier_h,"liste des .cc ="))  // liste des .cc trouvee
		{
 			int nombre;
			p_fichier_h>>nombre;
			for (int i=1;i<=nombre;i++)
			{
				prochain(p_fichier_h);
				string nom_fichier_cc;
				p_fichier_h>>nom_fichier_cc;

				int res=1;
				if(nom_fichier_cc[0] == '.') 
					if(nom_fichier_cc[1] == '/') // si "./chemin/toto.cc"
					{
						res=2;
						nom_fichier_cc.erase(0,2);  // enleve les caracteres ./ du debut
					}

				if(res==1)
					text += traite_cc(rep_src, nom_fichier_cc, 1);       // appelle (A)
				else if (res==2)
				{
					if(n_passage >=2 || opt_makegui==0 || nom_fichier_cc !="com.cc")  // com.h ne doit pas etre traite si existe pas
						text += traite_cc(rep_src2, nom_fichier_cc, 2);       // appelle (A)
				}	

			}
		}
		else
		{

			string chaine_rep2=nom_fichier_h.substr(0,nom_fichier_h.find(".h"))+".cc"; // on remplace le suffixe .h par .cc
			if(n_passage >=2 || opt_makegui==0 || chaine_rep2 !="com.cc")  // com.h ne doit pas etre traite si existe pas
				text += traite_cc(nom_rep_h, chaine_rep2, opt_rep);
		}
	} // fin partie 1


	{ // debut partie 2

// ouvre le fichier...........................
		ifstream p_fichier_h(nom_rep_fichier_h.c_str());
		if(!p_fichier_h)
			cerr<<"erreur ouverture du fichier :"<<nom_rep_fichier_h<<" dans (*) \n";


// recherche les fichiers titi.h  pour (2)
		int test(0);

		string chaine_h;
		string nom_fichier2_h;
		
		while(int res=recherche_h(p_fichier_h,nom_fichier2_h)) // -> nom_fichier2_h et pointeur p_fichier_h
		{


			if(verbose>=2)
				cout<<"res="<<res<<" nom_fichier2_h="<<nom_fichier2_h<<endl;

			if(n_passage == 1  && opt_makegui==1 &&  nom_fichier2_h =="com.h")
				continue;
			
			test=1;     // a trouve
			if(res==1)
			{
				text += traite_h(rep_src, nom_fichier2_h, 1);       // appelle (B)
				chaine_h += "$(SRC)/"+nom_fichier2_h + " ";  // pour (2)
			}
			else if(res==2) // dans repertoire du main
			{
			
					text += traite_h(rep_src2, nom_fichier2_h, 2);       // appelle (B)
					chaine_h += "$(SRC2)/"+nom_fichier2_h + " ";  // pour (2)
			
			}
			
			
		}

// rajoute la dependance de fichier.h .... pour (2)
		if(test)
		{
			text += "$(SRC)/" + nom_fichier_h + " : " + chaine_h + "\n";
			text += "\ttouch $@ \n\n";  //remet la date a jour
		}
	} //fin partie 2

	compteur--;

}

//==========================================================
/* fonction qui traite le cas d'un fichier toto.cc cad operation (A)
rep= 1: si repertoire src
rep =2 si repertoire src2
*/
string traite_cc(string nom_rep_cc, string  nom_fichier_cc, int opt_rep)
{
	string text;
	string chaine_o;
	string chaine_rep2;

	string nom_rep_fichier_cc = nom_rep_cc + "/" + nom_fichier_cc;
	if(verbose)
		cout<<" dans traite_cc:  nom_rep_fichier_cc="<< nom_rep_fichier_cc<<endl;


	if (deja_traite(nom_rep_fichier_cc))  // si le fichier a deja ete traite
	{
		cout<<"deja traite"<<endl;
		return text;
	}
	else
	{
		rajoute_liste(nom_rep_fichier_cc); //le rajoute a la liste
		if(opt_rep == 1)
			liste_files_h_cc.push_back(nom_fichier_cc);
		else if(opt_rep == 2)
			liste_files_main.push_back(nom_fichier_cc);
	}

	compteur++;

	if(verbose)
	{
		for (int i=1;i<=compteur;i++)
			cout<<" ";
		cout<<compteur<<" Traitement du fichier :"<<nom_rep_fichier_cc<<"\n";
	}

// rajoute a liste OBJET , pour (1)......
	string rep =  "$(SRC)/";
	if(opt_rep == 2) // si fichier contenant le main()
		rep = "$(SRC2)/";

	chaine_objets += rep  +  nom_fichier_cc.substr(0,nom_fichier_cc.find(".cc")) + ext_o + " "; // on remplace le suffixe .cc par .o ou .bc
	

// ouvre le fichier...........................
	ifstream p_fichier_cc(nom_rep_fichier_cc.c_str());
	if(!p_fichier_cc)
		cerr<<"erreur ouverture du fichier "<<nom_rep_fichier_cc<<" dans (**)\n";

//...pour (3)....

	chaine_o= rep + nom_fichier_cc.substr(0,nom_fichier_cc.find(".cc"))+ ext_o + " : " + rep  + nom_fichier_cc+" ";
	
// recherche les fichiers .h
	string nom_fichier_h;
	while(int res=recherche_h(p_fichier_cc, nom_fichier_h))
	{	

		if(verbose>=2)
				cout<<"res="<<res<<" nom_fichier_h="<<nom_fichier_h<<endl;

		if(n_passage ==1 && opt_makegui==1 &&  nom_fichier_h =="com.h")
				continue;
			
		if(res==1)
		{
			text += traite_h(rep_src, nom_fichier_h, 1);       // appelle (B)
			chaine_o += "$(SRC)/"+nom_fichier_h + " ";  // pour (3)
		}
		else if(res==2)
		{
//			cout<<"rep_src2="<<rep_src2<<endl;
			text += traite_h(rep_src2, nom_fichier_h, 2);       // appelle (B)
			chaine_o += "$(SRC2)/"+nom_fichier_h + " ";  // pour (3)
		}

	}

// rajoute la dependance de fichier .o ou .bc.... pour (3)
	text += chaine_o + "\n";
	text += "\t$(CC) $(CFLAGS)  -c  $< -o $@ \n\n"; // rem: $*  reconnait .o mais pas .bc. Alors on a mit $<
	//  rem : il faudrait  -c $<  -o $@

	compteur--;
}



//===================================================================
//Commande:   makef repertoire nomfichier.cc  makef.config    <nom_executable>
main(int argc,char *argv[])
{

	
	cout<<"=========================="<<endl;
	cout<<"MAKEF"<<endl;
	cout<<"--------"<<endl;


	
	compteur=0;
	rep_src2= string(argv[1]); // ex: /home/rep_essai
	src2 = string(argv[2]); // ex: essai.cc 
	string chaine_exec = src2.substr(0,src2.find(".cc")); // ex: essai
	string chaine_source = rep_src2 + "/" + src2;// ex: /home/rep_essai/essai.cc


//	cout<<"rep_src2="<<rep_src2<<endl;

	
	string nom_lib;
	if (argc >= 5)
		cout<<"\nL'executable choisi s'appellera:"<<(nom_lib=argv[4])<<endl;
	else
		cout<<"\nL'executable (par defaut) s'appellera:"<<(nom_lib=chaine_source.substr(0,chaine_source.find(".cc")))<<endl;

	string nom_config = rep_src2 + "/makef.config";
	cout<<endl<<"Recherche du fichier de parametres:"<<endl;
	string com = "ls " + rep_src2 + "/makef.config";
	int res = system(com.c_str());
	if(res != 0) // si fichier non trouve
	{
		nom_config =  argv[3];
		cout<<"on va donc prendre le fichier propose  par defaut:"<<nom_config<<endl;
	}
	else
		cout<<"ok."<<endl<<endl;

//**********
// lecture du fichier de config
	string rep_makef,librairies,chaine_comp;
	
	ifstream f_config(nom_config.c_str()); // fichier de configurations

	prochain(f_config); f_config>>rep_src;  // ex: /home/faure/c++/Utils

	prochain(f_config); f_config>>rep_makef; // ex:  /home/faure/c++/makef
	prochain(f_config); f_config>>verbose;
	prochain(f_config); f_config>>opt_makeh;
	prochain(f_config); f_config>>opt_makegui;
	
	prochain(f_config,"librairies et include  ="); f_config>>librairies;
	prochain(f_config,"compilation ="); f_config>>chaine_comp;

	int p;
	while(((p=librairies.find("_"))>=0) && (p<=librairies.size()))
		librairies.replace(p,1," ");  // remplace tous les "_" par " "
	while(((p=chaine_comp.find("_"))>=0) && (p<=chaine_comp.size()))
		chaine_comp.replace(p,1," ");  // remplace tous les "_" par " "
	while(((p=librairies.find("&"))>=0) && (p<=librairies.size()))
		librairies.replace(p,1,"_");  // remplace tous les "&" par "_"
	while(((p=chaine_comp.find("&"))>=0) && (p<=chaine_comp.size()))
		chaine_comp.replace(p,1,"_");  // remplace tous les "&" par "_"


	cout<<" opt_makeh= "<<opt_makeh<<endl;
	cout<<" opt_makegui= "<<opt_makegui<<endl;



    //--- detection du compilateur demandé : g++ ou emcc
	string ext_exec="";
	comp_name = chaine_comp.substr(0, chaine_comp.find(" "));
	comp_name = comp_name.substr(comp_name.rfind("/")+1); // -> g++ ou emcc
//	cout<<"comp_name="<<comp_name<<endl; 
	if(comp_name=="g++")
		ext_o = ".o";
	else if(comp_name=="emcc")
	{
		ext_o = ".bc";
		ext_exec = ".js";
	}
	else
	{
		cout<<"compilateur non reconnu:"<<comp_name<<endl;
		exit(0);
	}


	
//*****************

	string chaine1 = chaine_source.substr(0,chaine_source.rfind("/")) + "/Makefile";
	string chaine2 = chaine_source.substr(0,chaine_source.rfind("/")) + "/Makefile_tmp"; // pour zip
	string chaine3 = chaine_source.substr(0,chaine_source.rfind("/")) + "/README_tmp"; // pour zip

	

//	cout<<"\n On va creer un fichier Makefile:"<<endl<<chaine1<<endl;
	
	ofstream fichier_make(chaine1.c_str());
	ofstream fichier_make2(chaine2.c_str());

	ofstream fichier_README(chaine3.c_str());
	fichier_README<<"Pour compiler, ecrire dans un terminal: \tmake all \n";
	if(comp_name=="g++")
		fichier_README<<"Pour executer, ecrire dans un terminal: \t./" + chaine_exec + "\n\n";
	else if (comp_name=="emcc")
		fichier_README<<"Pour executer, ecrire dans un terminal: \tfirefox ./" + chaine_exec + ".html\n\n";
	fichier_README<<"Des informations sur le projet se trouvent dans le fichier: "<<src2<<"\n";
	fichier_README<<"\n\nConsignes pour installer c++ et les biblioteques necessaires sur votre ordinateur:\n";
	fichier_README<<" Voir: https://www-fourier.ujf-grenoble.fr/~faure/enseignement/c++/cours_c++/cours_1/cours_1.xhtml \n";
	fichier_README.close();

	
	//--------------------



	//..............
	string text1;
	
	text1 += "# Ce fichier Makefile a ete cree automatiquement. \n# Pour compiler, ecrire dans un terminal: \tmake all \n";
	text1 += "# Pour executer, ecrire dans un terminal:\t./" + chaine_exec + "\n\n";

	fichier_make << text1;
	fichier_make2 << text1;

	fichier_make << "# Pour fabriquer une librairie: \tmake lib \n";
	fichier_make << "# Pour exporter le projet: \tmake zip \n \n";

    //...........	
	text1 = "# -----repertoires (librairies, sources, et main) ----------------------\n\n";
	fichier_make << text1;
	fichier_make2 << text1;


	//......
	
	fichier_make<<"SRC ="<<rep_src<<endl;  // ex: /home/faure/c++/Utils
	fichier_make<<"SRC2 ="<<rep_src2<<endl<<endl; // ex: /home/rep_essai

	fichier_make2<<"SRC = ./Utils"<<endl;
	fichier_make2<<"SRC2 = ."<<endl;

	
	//..................
	string text2;
	text2 += "LIBR=" + librairies + "\n \n";

	text2 += "#-- instructions de commande,  utilitaires ---------------\n\nCC =" + chaine_comp + "  # compilateur\nCFLAGS =  -I$(SRC)       # options de compilation\n\n";

	text2 += "#============ cibles =======================================\n\n";


	
	//=======================================
	// si opt_makegui==1, il faut deux passages, car:
	// 1er passage: on ne rajoute pass com.h et com.cc a Makefile (car n'existent pas), mais on lit les infos et on ecrit les fichiers com.h et com.cc
	// 2eme passage: on  rajoute les fichiers com.h et com.cc à Makefile
	string text2_bis;
	
	for(n_passage =1; n_passage <= 1+opt_makegui; n_passage++) // 1 ou 2 passages selon opt_makegui
	{
		chaine_objets="OBJETS= ";
			
		//---- vide listes

		L_infos.clear();
		L_noms_classes.clear();
		L_noms_fichiers.clear();
		L_noms_menus.clear();
		L_noms_tabs.clear();
		liste_files_h_cc.clear();
		liste_files_main.clear();
		
//----
	
		
		//--- vide le fichier liste.txt ----------------------

		fliste = rep_makef+"/liste.txt";

		ofstream pliste(fliste.c_str()); // crée fichier
		pliste<<endl;
		pliste.close();
		
		//----------------------
	
		if(verbose)
			cout<<"####=== NUMERO DE PASSAGE:  n_passage = "<<n_passage<<endl;
		
		text2_bis = traite_cc(rep_src2, src2, 2);  // fichier .cc principal, contenant main()

		//================= Option on fabrique l'interface de commandes com.h et com.cc
		if(opt_makegui)
		{
			//-- etape 1: on extrait les informations des fichiers
			cout<<"###=========================="<<endl;
			cout<<"MAKE_GUI"<<endl;
			cout<<"--------"<<endl;
			cout<<"##== Part1: extraction des informations =========="<<endl;

			for(auto nom_fich : liste_files_h_cc)
				if(nom_fich != "com.cc" && nom_fich != "com.h")
					Extrait_infos_from_file(rep_src, nom_fich); 
	
			for(auto nom_fich : liste_files_main)
				if(nom_fich != "com.cc" && nom_fich != "com.h")
					Extrait_infos_from_file(rep_src2,  nom_fich); 

			if(verbose>=1)
				Bilan_collecte_infos();
		
			//-- etape 2  Ecriture des fichiers com.cc com.h
			Ecriture_fichiers_cpp(); // -> com.cc et com.h
		
		}
		for(auto nc: L_noms_classes)
			cout<<"**** nc="<<nc<<endl;

	} // for n_passage

	text2 += text2_bis;
	
	//===============================================
	text2 += chaine_objets + "\n";  // finit (1)



	text2 += "\n#-- fabrique executable final ------------- \n\nall:  $(OBJETS) $(SRC2)/" + chaine_exec + ext_o + "  \n\t$(CC) $(CFLAGS)  -o $(SRC2)/" + chaine_exec + ext_exec + "  $(OBJETS)     $(LIBR) \n";


	fichier_make<<text2;
	fichier_make2<<text2<<endl;

	fichier_make<<"\tmake zip\n\n"; // fabrique aussi le projet zip



	fichier_make << "\n#--- fabrique une librairie ------------------ \n\nlib: $(OBJETS)\n\tar r $(SRC2)/" << chaine_exec << "  $(OBJETS) \n\n";



	
    //======= Fabrication de l'exportation project.zip============
	string nom_zip = "project_c++_" + chaine_exec + ".zip"; // ex: project_c++_essai.zip
	
	string text3;

	
	text3 += "\n#--- fabrique un projet exportable .zip \n\n";
	text3 += "zip: $(OBJETS)\n";
	//.. rem: && sert a enchainer les instructions comme ; mais garanti que le faire que si la precedente est OK
	//.. rem: \ permet de revenir à la ligne
	text3 += "\tcd $(SRC) && cd .. && \\\n";  // se place dans /home/faure/c++
	text3 += "\trm -rf tmp && \\\n"; // efface repertoire precedant par prevention
	text3 += "\tmkdir tmp && \\\n"; //cree nouveau repertoire

	//... rajoute les fichiers .h et .cc
	for(int i=0; i<liste_files_h_cc.size(); i++)
	{
		text3 += "\tzip tmp/" + nom_zip + " Utils/" + liste_files_h_cc[i];
		if(i<liste_files_h_cc.size()-1)
			text3 += " && \\\n";
		else
			text3 += "\n";
	}

	

	text3 += "\tcd $(SRC2) && \\\n";  // se place dans /home/rep_essai
	text3 += "\trm -rf " +  nom_zip + " && \\\n"; // detruit eventuel fichier par prevention

	if(liste_files_h_cc.size()>0)
		text3 += "\tmv $(SRC)/../tmp/"  + nom_zip + " " + rep_src2 + " && \\\n"; //deplace le projet 

	
	//... rajoute les fichiers du main()
	for(int i=0; i<liste_files_main.size(); i++)
	{
		text3 += "\tzip " + nom_zip + " " + liste_files_main[i];
		text3 += " && \\\n";
	}

	text3 += "\tzip " + nom_zip + " Makefile_tmp  && \\\n";

	text3 += "\tprintf \"@ Makefile_tmp\\n@=Makefile\\n\" | zipnote -w "  + nom_zip + " && \\\n";	 // change le nom Makefile_tmp -> Makefile dans l'archive projet.zip

	text3 += "\trm Makefile_tmp  && \\\n";

	text3 += "\tzip " + nom_zip + " README_tmp  && \\\n";

	text3 += "\tprintf \"@ README_tmp\\n@=README\\n\" | zipnote -w " + nom_zip + " && \\\n";	 // change le nom README_tmp -> README dans l'archive projet.zip

	text3 += "\trm README_tmp\n";
	
	fichier_make<<text3;

	
	//---------------------
	cout<<"Fabrication du Makefile terminee."<<endl<<endl;
	cout<<"=========================="<<endl;
	
	return 0;
}


