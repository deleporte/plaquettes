/* test_root.cpp


compilation : g++ test_root.cpp -o test_root -lm -I/usr/include/root `root-config --cflags` `root-config --libs` `root-config --glibs`

 */






#include <TApplication.h> // (A) 
#include <TCanvas.h> //(B) //Erreur Alix: non trouvée
#include <TEllipse.h> // (C) //Erreur Alix: non trouvée
 
//------ fonction main (A) --------------------
int main()
{
 TApplication theApp("App", nullptr, nullptr); // (A)
 
 TCanvas c("c","fenetre",400,400); // (B) objet fenetre graphique. Taille en pixels
 
 c.Range(0,0,5,5); // coordonnees de la fenetre c: x:0-5, y:0-5 (optionnel, par defaut ce sera 0-1)
 
//------ dessin d'une ellipse dans la fenetre c (C) -------
 TEllipse e(2,3,1); // on precise le centre (x=2,y=3) et le rayon=1
 e.Draw(); // dessine l'ellipse
 c.Update(); //(B) Montre le dessin
 
 theApp.Run(); // (A) garde la fenetre ouverte et permet a l'utilisateur d'interagir.
}
