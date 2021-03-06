#ifndef __INGANAVERTEX_HXX__
#define __INGANAVERTEX_HXX__

// ROOT includes
#include "TApplication.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TCanvas.h"
#include <TStyle.h>
#include "TString.h"
#include "TSystem.h"
#include "TSpectrum.h"
#include "TTree.h"
//C++ libraly includes
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip.h>
#include <algorithm>
#include <deque>
#include <vector>
#include <sys/stat.h>
#include <unistd.h> // using getopt      
using namespace std;
// INGRID includes


//__________________________________________________________

int pln_th  =2;
float ch_th =150;

void fIngVertex(int mod, int ZCriteria, double XYCriteria){//seules les traces ayant vtxtrk true pourront passer ensuite comme ingbasrec. On teste ici si 2 traces sont issues du meme vertex. A la base, toute trace 3D qui sort de ingana.hxx a vtxtrk=true (on a vérifié les conditions ou les plans de départ ne sont pas trop éloignés et que les plans d'arrivées de chaque trace 2D ne dépassent pas les plans de départ de l'autre trace 2D+cluster en temps). Mais on a pas testé de les assembler avec plusieurs autres traces!
  pln_th=ZCriteria;
  ch_th=XYCriteria;
  //cout<<"pln_th="<<pln_th<<", ch_th="<<ch_th<<endl;
  for(int i=0;i<ingtrack.size();i++){//déja ça montre bien ce qu'on a vu dans IngAna.hxx => chaque ingtrack[i]=1trace 3D => Ntrack=#traces 3D
    ingtrack[i].Ntrack=1;
  }

  for(int i=0;i<ingtrack.size();i++){//on parcourt toutes le traces du event/module/cycle
    for(int j=i+1;j<ingtrack.size();j++){//va rega&rder toutes les traces suivantes
      /*
      if(abs((ingtrack[i].startxpln)-(ingtrack[j].startxpln))>pln_th)continue;
      if(abs((ingtrack[i].startypln)-(ingtrack[j].startypln))>pln_th)continue;
      if(fabs((ingtrack[i].x)-(ingtrack[j].x))>ch_th)continue;
      if(fabs((ingtrack[i].y)-(ingtrack[j].y))>ch_th)continue;
      */
      
      if(abs((ingtrack[i].startxpln)-(ingtrack[j].startxpln))+abs((ingtrack[i].startypln)-(ingtrack[j].startypln))>pln_th)continue;
      if(fabs((ingtrack[i].y)-(ingtrack[j].y)+fabs((ingtrack[i].x)-(ingtrack[j].x)))>ch_th)continue;//la on teste la proximité entre les vertex => Si on a par exemple une intéraction très inelastique, on va avoir des hits partout, avec parfois des cases manquantes: on va parfois reconstruire 2 vertex différents dans ce cas alors que l'on avait la meme inté en fait => 2 IngBasRec différentes pour le meme vertex. Ici on regarde si 2 traces 3D sont suffisamment proches. Si ce n'est pas le cas, pourquoi ne met on pas vtxtrk=false? Car ce peut etre 2 traces 3D bien disctinctes?? donc on recontruit les 2 séparément!!!!!

      if(!(ingtrack[i].vtxtrk)||!(ingtrack[j].vtxtrk))continue;//cela marque si l'une est en gros véto et déja utilisées dans une autre trace

      if((ingtrack[i].vetowtracking||ingtrack[i].edgewtracking)&&//si les 2 traces 3D sont véto et si la trace j est la plus downstream => augmenter la trace la plus downstream d'une trace en plus dans ntrack, mais ne lui ajoute pas toutes les infos de la trace (on a ni ses hits, ni rien d'autre). Seule la trace la plus downstream sera reconstruite, avec un poids de 2 (ou pls) traces avec veretex commun, et ne contenant l'info que d'une seule trace
	 (ingtrack[j].vetowtracking||ingtrack[j].edgewtracking)&&
	 (ingtrack[i].endxpln+ingtrack[i].endypln)<(ingtrack[j].endxpln+ingtrack[j].endypln)){
	ingtrack[i].vtxtrk = false;
	for(int ihit=0;ihit<ingtrack[i].hitnum.size();ihit++){
	  ingtrack[j].hitnum.push_back(ingtrack[i].hitnum[0]);  
	}
	ingtrack[j].Ntrack += ingtrack[i].Ntrack; 
      }
      else if((ingtrack[i].vetowtracking||ingtrack[i].edgewtracking)&&//si les 2 traces sont véto mais que c'est la première la plus downstream(le elseif dit ça) alors c'est la plus downstream qui va etre représentante de la congrégation de traces issues du meme vertex
	      (ingtrack[j].vetowtracking||ingtrack[j].edgewtracking)){
	ingtrack[j].vtxtrk = false;
	ingtrack[i].Ntrack += ingtrack[j].Ntrack;
        for(int ihit=0;ihit<ingtrack[j].hitnum.size();ihit++){
          ingtrack[i].hitnum.push_back(ingtrack[j].hitnum[0]);
        }
      }     
      else if(ingtrack[j].vetowtracking||ingtrack[j].edgewtracking){//si j est upstreamvéto et/ou edgevéto => c'est elle qui va porter la voix de l'ensemble des traces. Pourquoi exactement cette mesure?????
	ingtrack[i].vtxtrk = false;
	ingtrack[j].Ntrack += ingtrack[i].Ntrack;
        for(int ihit=0;ihit<ingtrack[i].hitnum.size();ihit++){
          ingtrack[j].hitnum.push_back(ingtrack[i].hitnum[0]);
        }
      }
      else if(ingtrack[i].vetowtracking||ingtrack[i].edgewtracking){//de meme si c'est i
	ingtrack[j].vtxtrk = false;
	ingtrack[i].Ntrack += ingtrack[j].Ntrack;
	for(int ihit=0;ihit<ingtrack[j].hitnum.size();ihit++){
          ingtrack[i].hitnum.push_back(ingtrack[j].hitnum[0]);
        }
      }
      else if((ingtrack[i].endxpln+ingtrack[i].endypln)<(ingtrack[j].endxpln+ingtrack[j].endypln)){//dans le cas ou j est la plus downstream, c'est elle qui porte la culotte
	ingtrack[i].vtxtrk = false;
	ingtrack[j].Ntrack += ingtrack[i].Ntrack;
        for(int ihit=0;ihit<ingtrack[i].hitnum.size();ihit++){
          ingtrack[j].hitnum.push_back(ingtrack[i].hitnum[0]);
        }
      }
      else{//si c'est i la plus down, cette fois c'est elle qui porte la culotte
	ingtrack[j].vtxtrk = false;
	ingtrack[i].Ntrack += ingtrack[j].Ntrack;
	for(int ihit=0;ihit<ingtrack[j].hitnum.size();ihit++){
          ingtrack[i].hitnum.push_back(ingtrack[j].hitnum[0]);
        }
      }
      
    }
  }

  /*conclusions:*2 traces 3D dont les vertex ne matchent pas sont construites dans 2 ingbasrec différentes
                *2 traces 3D dont les vertex matchent: une seule est mise dans ingbasrec avec un NTrack=2 (la plus downstream ou la véto)
		*Cela se généralise a plus de 2 traces évidemment
donc dans chaque ingtrack[i] on a une trace 3D au début. Puis on va légerement changer la plus downstream des traces 3D en lui ajoutant les autres et en la mettant la seule a avoir vtxtrk=true
		*/


};



/*


void fIngVertex(int mod){
  vector<IngTrack>::iterator it1;
  vector<IngTrack>::iterator it2;

  for(it1=ingtrack.begin();it1<ingtrack.end();it1++){
    for(it2=it1+1;it2<ingtrack.end();it2++){
      if(abs((it1->startxpln)-(it2->startxpln))>pln_th)continue;
      if(abs((it1->startypln)-(it2->startypln))>pln_th)continue;
      if(fabs((it1->x)-(it2->x))>ch_th)continue;
      if(fabs((it1->y)-(it2->y))>ch_th)continue;

      if((it1->vetowtracking||it1->edgewtracking)&&
	 (it2->vetowtracking||it2->edgewtracking)&&
	 (it1->endxpln+it1->endypln)<(it2->endxpln+it2->endypln)){
	it1 = ingtrack.erase(it1);
	it1--;
      }
      else if((it1->vetowtracking||it1->edgewtracking)&&
	      (it2->vetowtracking||it2->edgewtracking)){
	it2 = ingtrack.erase(it2);
	it2--;
      }     
      else if(it2->vetowtracking||it2->edgewtracking){
	it1 = ingtrack.erase(it1);
	it1--;
      }
      else if(it1->vetowtracking||it1->edgewtracking){
	it2 = ingtrack.erase(it2);
	it2--;
      }
      else if((it1->endxpln+it1->endypln)<(it2->endxpln+it2->endypln)){
	it1 = ingtrack.erase(it1);
	it1--;
      }
      else{
	it2 = ingtrack.erase(it2);
	it2--;
      }
    }
  }
};
*/

#endif
