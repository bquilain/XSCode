//##### Standard C++ lib. ######
#include<iostream>
#include<sstream>
#include<fstream>
using namespace std;
#include <iomanip.h>
#include <sys/stat.h>
//##### Root Library ###########
#include <TROOT.h>
#include <TStyle.h>
#include <TH1.h>
#include <TH2.h>

#include <TApplication.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TObject.h>
#include <TEventList.h>
#include <TBranch.h>
#include <TSystem.h>
#include <TBrowser.h>
//##### INGRID Library #########
#include "INGRIDEVENTSUMMARY.h"
#include "IngridHitSummary.h"
#include "IngridSimHitSummary.h"
#include "IngridBasicReconSummary.h"
#include "IngridConstants.h"

TFile* fTFile;
TH1F*  fHVerX  [nMod];
TH1F*  fHVerY  [nMod];
TH1F*  fHVerXZ [nMod];
TH1F*  fHVerYZ [nMod];
TH1F*  fHVerXdiv [nMod][2];
TH1F*  fHVerYdiv [nMod][2];
TH1F*  fHVerXZdiv[nMod][2];
TH1F*  fHVerYZdiv[nMod][2];
void Book(string str){
  fTFile = new TFile(str.c_str(), "recreate");
  for(int imod=0; imod<nMod; imod++){
    fHVerX[imod]  = new TH1F( Form("fHVerX_Mod%02d", imod),
			      Form("vertex X of Mod%02d", imod),
			      24, 0, 24 );
    fHVerY[imod]  = new TH1F( Form("fHVerY_Mod%02d", imod),
			      Form("vertex Y of Mod%02d", imod),
			      24, 0, 24 );


    fHVerX [imod] -> SetMinimum(0);
    fHVerY [imod] -> SetMinimum(0);
    fHVerX [imod] -> SetXTitle ("vertex X[ch#]");
    fHVerY [imod] -> SetXTitle ("vertex Y[ch#]");
    fHVerX [imod] -> SetYTitle ("# of events");
    fHVerY [imod] -> SetYTitle ("# of events");

    for(int id=0; id<2; id++){
      char name[300];
      if(id==0)sprintf(name, "vertex z=0~3"); 
      if(id==1)sprintf(name, "vertex z=4~"); 
      fHVerXdiv[imod][id]  = new TH1F( Form("fHVerX_Mod%02d_%d", imod, id),
				       Form("vertex X of Mod%02d(%s)", imod, name),
				       24, 0, 24 );
      fHVerYdiv[imod][id]  = new TH1F( Form("fHVerY_Mod%02d_%d", imod, id),
				       Form("vertex Y of Mod%02d(%s)", imod, name),
				       24, 0, 24 );


      fHVerXdiv [imod][id] -> SetMinimum(0);
      fHVerYdiv [imod][id] -> SetMinimum(0);
      fHVerXdiv [imod][id] -> SetXTitle ("vertex X[ch#]");
      fHVerYdiv [imod][id] -> SetXTitle ("vertex Y[ch#]");
      fHVerXdiv [imod][id] -> SetYTitle ("# of events");
      fHVerYdiv [imod][id] -> SetYTitle ("# of events");

    }
  }
}

void Fill(IngridBasicReconSummary* basicrecon){
  int imod  = basicrecon -> hitmod;
  int verx  = basicrecon -> vertexx [0];
  int very  = basicrecon -> vertexy [0];
  int verxz = basicrecon -> vertexxz;
  int veryz = basicrecon -> vertexyz;
  fHVerX [ imod ] -> Fill( verx  );
  fHVerY [ imod ] -> Fill( very  );
  int id;
  if( basicrecon -> vertexxz < 4  )
    id = 0;
  else
    id = 1;
  fHVerXdiv [ imod ][ id ] -> Fill( verx  );
  fHVerYdiv [ imod ][ id ] -> Fill( very  );
  
  

}

void Write(){
  fTFile -> Write();
  fTFile -> Close();
}

bool FidX(IngridBasicReconSummary* basicrecon, int ich){

  if( 23 - ich >= basicrecon -> vertexx[0] &&
      0 + ich <= basicrecon -> vertexx[0] )
    return true;
  else
    return false;

}

bool FidY(IngridBasicReconSummary* basicrecon, int ich){

  if( 23 - ich >= basicrecon -> vertexy[0] &&
      0 + ich <= basicrecon -> vertexy[0] )
    return true;
  else
    return false;

}
