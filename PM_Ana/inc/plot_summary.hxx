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


#define  FILLCYC 100
const static long int   minpot  =  5e14;   //## Min. value of pot on History 2D plot
const static long int   maxpot  =  1.1e18; //## Max. value of pot on History 2D plot
const static long int  nbinpot  =   2e3;   //## # of Bins of pot axis of History 2D plot 
const static int        minEvt  =     0;   //## Min. value of events on History 2D plot
const static int        maxEvt  = 73000;   //## Max. value of events on History 2D plot
const static int       nbinEvt  =  2000;   //## # of Bins of event axis of History 2D plot 
const static long int nbintime  =  4000;   //## # of Bins of time axis of History 2D plot 
                                           //## Min. and Max. time value will be read in main()
const static int    minexptime  = -100;
const static int    maxexptime  =  100;
const static int   nbinexptime  =  101;


//####################
TFile*                        rfile;
TTree*                        rtree;
TBranch*                      EvtBr;
IngridBasicReconSummary* basicrecon;
