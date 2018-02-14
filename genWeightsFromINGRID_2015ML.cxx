//
//
// Example to show how to reweight events from INGRID file (SK tree)
//
// ./generateWeightsFromINGRID.exe -i ingrid_inputfile [optional nevents] 
// 
// M.Licciardi  - Sep.06/2017 . Write code based on genWeightsFromINGRID_2015 from B.Quilain with updated values from TN265

#include <stdlib.h>
#include <cstdlib>

#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TString.h"

#include "T2KReWeight.h"
#include "T2KSyst.h"

#include "T2KGenieReWeight.h" 
#include "T2KGenieUtils.h"

#include "T2KNeutReWeight.h"
#include "T2KNeutUtils.h"

#include "T2KJNuBeamReWeight.h"

#include "T2KNIWGReWeight.h"
#include "T2KNIWGUtils.h"

#include "SK__h1.h"
// For weight storer class
#include "T2KWeightsStorer.h"


using std::cout;
using std::cerr;
using std::endl;

using namespace t2krew;

int fNskEvts = -1;
TString fSKFileName;
TString fOutFileName;

void Usage();
int ParseArgs(int argc, char *argv[]);

int main(int argc, char *argv[])
{
 
  // process the arguments
  /* int args = ParseArgs(argc, argv);
     if(args != 0){
     std::cerr << "Usage " << std::endl;
     return 0;
     }*/
  int c=-1;
  bool tune=false;

  while ((c = getopt(argc, argv, "i:o:t")) != -1) {
    switch(c){
    case 'i':
      fSKFileName=optarg;
      break;
    case 't':
      tune=true;
      break;
    case 'o':
      fOutFileName=optarg;
      break;
    }
  }

  //T2KWeightsStorer *storer = new T2KWeightsStorer("testweights_skexample_nue.root"); //forweightstorer
  T2KWeightsStorer *storer = new T2KWeightsStorer(fOutFileName); //forweightstorer

  cout << "Starting to reweight NIWG events from SK file: " << fSKFileName << endl;


  // Load the SK "h1" tree
  TFile * SK_infile = new TFile(fSKFileName, "OPEN");
  if(!SK_infile){
    cerr << "Cannot open SK file!" << endl;
    exit(1);
  }
  TTree * SK_tree = (TTree*) SK_infile->Get("h1");
  if(!SK_tree){
    cerr << "Cannot find SK_tree!" << endl;
  }
 
  // Instantiate the reader of the SK tree class
  SK::SK__h1 *skVtxs = new SK::SK__h1(SK_tree,1);

  const int ndraws = 7*20 + 7; // -3 sigma ... +3 sigma (7 variations, 18 dials + 9 tunings - 7 individual, one with all CCQE tunings, and one where all tunings are applied)

  cout<<"before NIWG ReWeight"<<endl;

#ifdef __T2KRW_NIWG_ENABLED__
#ifdef __T2KRW_NEUT_ENABLED__
  cout<<"NIWG ReWeight"<<endl;
  // Make a t2kreweighting object and add a NIWG weighting engine. 
  t2krew::T2KReWeight rw; 

  rw.AdoptWghtEngine("niwg_rw", new t2krew::T2KNIWGReWeight());
  rw.AdoptWghtEngine("neut_rw", new t2krew::T2KNeutReWeight());

  // NIWG 2015
  // Tuning for CCQE
  rw.Systematics().Include(t2krew::kNIWG2014a_SF_RFG);
  rw.Systematics().Include(t2krew::kNIWG_rpaCCQE_norm);
  rw.Systematics().Include(t2krew::kNIWG_rpaCCQE_shape);

  // Uncertainties
  // CCQE:
  rw.Systematics().Include(t2krew::kNIWG2014a_pF_C12);
  rw.Systematics().Include(t2krew::kNIWG2014a_pF_O16);
  rw.Systematics().Include(t2krew::kNIWG2014a_Eb_C12);
  rw.Systematics().Include(t2krew::kNIWG2014a_Eb_O16);
  rw.Systematics().Include(t2krew::kNIWGMEC_Norm_C12);
  rw.Systematics().Include(t2krew::kNIWGMEC_Norm_O16);
  rw.Systematics().Include(t2krew::kNXSec_MaCCQE);
  rw.Systematics().Include(t2krew::kNXSec_VecFFCCQE);

  // CC and NC single pion resonance:
  rw.Systematics().Include(t2krew::kNXSec_CA5RES);
  rw.Systematics().Include(t2krew::kNXSec_MaNFFRES);
  rw.Systematics().Include(t2krew::kNXSec_BgSclRES);

  // nue/numu uncertainties -- not needed a priori
  //rw.Systematics().Include(t2krew::kNXSec_SCCVecQE);
  //rw.Systematics().Include(t2krew::kNXSec_SCCAxlQE);
  //rw.Systematics().Include(t2krew::kNXSec_PsFF);
  //rw.Systematics().Include(t2krew::kNIWG2012a_ccnueE0);

  // All other CC and NC
  rw.Systematics().Include(t2krew::kNIWG2012a_dismpishp);
  rw.Systematics().Include(t2krew::kNIWG2012a_cccohE0);
  rw.Systematics().Include(t2krew::kNIWG2012a_nccohE0);
  rw.Systematics().Include(t2krew::kNIWG2012a_ncotherE0);

  // Final State Interations
  rw.Systematics().Include(t2krew::kNCasc_FrAbs_pi);
  rw.Systematics().Include(t2krew::kNCasc_FrCExLow_pi);
  rw.Systematics().Include(t2krew::kNCasc_FrCExHigh_pi);
  rw.Systematics().Include(t2krew::kNCasc_FrInelLow_pi);
  rw.Systematics().Include(t2krew::kNCasc_FrInelHigh_pi);
  rw.Systematics().Include(t2krew::kNCasc_FrPiProd_pi);


  // Absolute tweak dials set the fractional uncertainty, instead of
  // in units of "sigma", defined in the code.
  // Useful so that you define the uncertainty within the code, as what is
  // hardcoded may not be the same as what is used for analysis.

  // CCQE:
  rw.Systematics().SetAbsTwk(t2krew::kNIWG2014a_pF_C12);
  rw.Systematics().SetAbsTwk(t2krew::kNIWG2014a_pF_O16);
  rw.Systematics().SetAbsTwk(t2krew::kNIWG2014a_Eb_C12);
  rw.Systematics().SetAbsTwk(t2krew::kNIWG2014a_Eb_O16);
  rw.Systematics().SetAbsTwk(t2krew::kNIWGMEC_Norm_C12);
  rw.Systematics().SetAbsTwk(t2krew::kNIWGMEC_Norm_O16);
  rw.Systematics().SetAbsTwk(t2krew::kNXSec_MaCCQE);

  // CC and NC single pion resonance:
  rw.Systematics().SetAbsTwk(t2krew::kNXSec_CA5RES);
  rw.Systematics().SetAbsTwk(t2krew::kNXSec_MaNFFRES);
  rw.Systematics().SetAbsTwk(t2krew::kNXSec_BgSclRES);

  // nue/numu uncertainties
  //rw.Systematics().SetAbsTwk(t2krew::kNXSec_SCCVecQE);
  //rw.Systematics().SetAbsTwk(t2krew::kNXSec_SCCAxlQE);
  //rw.Systematics().SetAbsTwk(t2krew::kNXSec_PsFF);
  //rw.Systematics().SetAbsTwk(t2krew::kNIWG2012a_ccnueE0);

  // All other CC and NC
  rw.Systematics().SetAbsTwk(t2krew::kNIWG2012a_dismpishp);
  rw.Systematics().SetAbsTwk(t2krew::kNIWG2012a_cccohE0);
  rw.Systematics().SetAbsTwk(t2krew::kNIWG2012a_nccohE0);
  rw.Systematics().SetAbsTwk(t2krew::kNIWG2012a_ncotherE0);

  // Final State Interations
  rw.Systematics().SetAbsTwk(t2krew::kNCasc_FrAbs_pi);
  rw.Systematics().SetAbsTwk(t2krew::kNCasc_FrCExLow_pi);
  rw.Systematics().SetAbsTwk(t2krew::kNCasc_FrCExHigh_pi);
  rw.Systematics().SetAbsTwk(t2krew::kNCasc_FrInelLow_pi);
  rw.Systematics().SetAbsTwk(t2krew::kNCasc_FrInelHigh_pi);
  rw.Systematics().SetAbsTwk(t2krew::kNCasc_FrPiProd_pi);


  // RPA tuning takes AbsTwk too
  rw.Systematics().SetAbsTwk(t2krew::kNIWG_rpaCCQE_norm);
  rw.Systematics().SetAbsTwk(t2krew::kNIWG_rpaCCQE_shape);

  for(int dial=0; dial<ndraws; dial++){

    cout<<"dial="<<dial<<endl;
    Double_t pfcdial = 0.;
    Double_t pfodial = 0.;
    Double_t ebcdial = 0.;
    Double_t ebodial = 0.;
    Double_t meccdial = 0.;
    Double_t mecodial = 0.;
    Double_t maqedial = 0.;

    Double_t cadial = 0.;
    Double_t manresdial = 0.;
    Double_t bkgdial = 0.;

    /*    Double_t sccvdial = 0.;
	  Double_t sccadial = 0.;
	  Double_t fpdial = 0.;
	  Double_t ccnuedial = 0.;*/

    Double_t dismpidial = 0.;
    Double_t cccohdial = 0.;
    Double_t nccohdial = 0.;
    Double_t ncothdial = 0.;

    Double_t absdial = 0.;
    Double_t cxlodial = 0.;
    Double_t cxhidial = 0.;
    Double_t inelodial = 0.;
    Double_t inehidial = 0.;
    Double_t piprodial = 0.;


    if(dial>=0 && dial<140){//redefined the nominal as the tune
      if(tune){
	pfcdial = 0.0276;
	maqedial = -0.0496;
	meccdial = -0.73;
      }
      rw.Systematics().SetTwkDial(t2krew::kNIWG2014a_SF_RFG,1); // SF->RFG tuning
      rw.Systematics().SetTwkDial(t2krew::kNIWG_rpaCCQE_norm,1); // add RPA
      rw.Systematics().SetTwkDial(t2krew::kNIWG_rpaCCQE_shape,0); // rel RPA
    }else{
      rw.Systematics().SetTwkDial(t2krew::kNIWG2014a_SF_RFG,0); 
      rw.Systematics().SetTwkDial(t2krew::kNIWG_rpaCCQE_norm,0);
      rw.Systematics().SetTwkDial(t2krew::kNIWG_rpaCCQE_shape,0);
    }


    if(dial>=0&&dial<=6){
      pfcdial += (dial-3-7*0.)*(tune?0.143:0.0737); //**** Nominal is 217 MeV/c for C12.  Fit result (TN192 Table13) was 223 +/- 31 MeV/c. Apply fractional change of 31/217 = 0.143
    }else if(dial>=7&&dial<=13){
      pfodial = (dial-3-7*1.)*(tune?0.143:0.0737); //**** Nominal is 225 MeV/c for O16 and same fractional error as C12. 
    }else if(dial>=14&&dial<=20){
      ebcdial = (dial-3-7*2.)*0.36; //Nominal is 25+/-9 MeV for C12. 9/25= 0.36 fractional error
    }else if(dial>=21&&dial<=27){
      ebodial = (dial-3-7*3.)*0.33333; //Nominal is 27+/-9 MeV for O16. 9/27= 0.3333 fractional error
    }else if(dial>=28&&dial<=34){
      meccdial += (dial-3-7*4.)*(tune?.73:1.); //***Nominal is 1. Use 73% error
    }else if(dial>=35&&dial<=41){
      mecodial = (dial-3-7*5.)*1; //****Nominal is 1. Use 100% error
    }else if(dial>=42&&dial<=48){
      maqedial += (dial-3-7*6.)*(tune?0.1565:0.1488); //****Nominal is 1.21. Fit result is 1.15+/-0.18. Apply fractional change of 0.18/1.15 = 0.1565 or 0.18/1.21 = 0.1488
    }else if(dial>=49&&dial<=55){
      cadial = (dial-3-7*7.)*0.1188; // Fit result is 1.01+/-0.12
    }else if(dial>=56&&dial<=62){
      manresdial = (dial-3-7*8.)*0.1579; // Fit result is 0.95+/-0.15. 
    }else if(dial>=63&&dial<=69){
      bkgdial = (dial-3-7*9.)*0.1538; // Fit result is 1.3+/-0.2. 
    }else if(dial>=70&&dial<=76){
      dismpidial = (dial-3-7*10.)*0.4; // Same as 2012/2013. 0+/-0.4
    }else if(dial>=77&&dial<=83){
      cccohdial = (dial-3-7*11.)*.3; //1+/-0.3
    }else if(dial>=84&&dial<=90){
      nccohdial = (dial-3-7*12.)*.3; // Same as 2012/2013. 1+/-0.30
    }else if(dial>=91&&dial<=97){
      ncothdial = (dial-3-7*13.)*.3; // Same as 2012/2013. 1+/-0.30
    }else if(dial>=98&&dial<=104){
      absdial = (dial-3-7*14)*0.5;
    }else if(dial>=105&&dial<=111){
      cxlodial = (dial-3-7*15)*0.5;
    }else if(dial>=112&&dial<=118){
      cxhidial = (dial-3-7*16)*.3;
    }else if(dial>=119&&dial<=125){
      inelodial = (dial-3-7*17)*.5;
    }else if(dial>=126&&dial<=132){
      inehidial = (dial-3-7*18)*.3;
    }else if(dial>=133&&dial<=139){
      piprodial = (dial-3-7*19)*.5;
    }

    // Special dials for applying NIWG 2015 tuning 
    if(dial>=144){
      maqedial = -0.0496; // MaQE tuning 1.21 -> 1.15
    }
    if(dial>=145){
      pfcdial = 0.0276; // pF_C tuning 217 -> 223
    }
    if(dial==146){
      meccdial = -.73; // MEC tuning 100 -> 27
    }
    
    /*
    // Special dials for applying NIWG 2015 tuning (Note: MEC_O and pF_O instead of MEC_C and pF_C, but using the same fractional change)
    else if(dial==129){
    maqedial = -0.15; // MaQE tuning: 1.2 -> 1.02 (1.02/1.2 - 1 = -0.15)
    }
    else if(dial==130){
    mecodial = -0.42; // MEC_norm_C tuning: 1 -> 0.58 (0.58/1 - 1 = -0.42)
    }
    else if(dial==131){
    pfodial = 0.1014; // pF_C tuning: 217 -> 239 (239/217 - 1 = 0.1014)
    }
    else if(dial==132){
    maqedial = -0.15; // All 3 of the above tunings applied together
    mecodial = -0.42;
    pfodial = 0.1014;
    }
    else if(dial==133){
    maqedial = -0.15; // All tunings applied (including the three below)
    mecodial = -0.42;
    pfodial = 0.1014;
    }
    */

    rw.Systematics().SetTwkDial(t2krew::kNIWG2014a_pF_C12, pfcdial);
    rw.Systematics().SetTwkDial(t2krew::kNIWG2014a_pF_O16, pfodial);
    rw.Systematics().SetTwkDial(t2krew::kNIWG2014a_Eb_C12, ebcdial);
    rw.Systematics().SetTwkDial(t2krew::kNIWG2014a_Eb_O16, ebodial);
    rw.Systematics().SetTwkDial(t2krew::kNIWGMEC_Norm_C12, meccdial);
    rw.Systematics().SetTwkDial(t2krew::kNIWGMEC_Norm_O16, mecodial);
    rw.Systematics().SetTwkDial(t2krew::kNXSec_MaCCQE, maqedial);
    rw.Systematics().SetTwkDial(t2krew::kNXSec_VecFFCCQE, 2); 

    rw.Systematics().SetTwkDial(t2krew::kNXSec_CA5RES, cadial);
    rw.Systematics().SetTwkDial(t2krew::kNXSec_MaNFFRES, manresdial);
    rw.Systematics().SetTwkDial(t2krew::kNXSec_BgSclRES, bkgdial);

    /*    rw.Systematics().SetTwkDial(t2krew::kNXSec_SCCVecQE, sccvdial);
	  rw.Systematics().SetTwkDial(t2krew::kNXSec_SCCAxlQE, sccadial);
	  rw.Systematics().SetTwkDial(t2krew::kNXSec_PsFF, fpdial);
	  rw.Systematics().SetTwkDial(t2krew::kNIWG2012a_ccnueE0, ccnuedial);*/
    
    rw.Systematics().SetTwkDial(t2krew::kNIWG2012a_dismpishp, dismpidial);
    rw.Systematics().SetTwkDial(t2krew::kNIWG2012a_cccohE0, cccohdial);
    rw.Systematics().SetTwkDial(t2krew::kNIWG2012a_nccohE0, nccohdial);
    rw.Systematics().SetTwkDial(t2krew::kNIWG2012a_ncotherE0, ncothdial);

    rw.Systematics().SetTwkDial(t2krew::kNCasc_FrAbs_pi, absdial);
    rw.Systematics().SetTwkDial(t2krew::kNCasc_FrCExLow_pi, cxlodial);
    rw.Systematics().SetTwkDial(t2krew::kNCasc_FrCExHigh_pi, cxhidial);
    rw.Systematics().SetTwkDial(t2krew::kNCasc_FrInelLow_pi, inelodial);
    rw.Systematics().SetTwkDial(t2krew::kNCasc_FrInelHigh_pi, inehidial);
    rw.Systematics().SetTwkDial(t2krew::kNCasc_FrPiProd_pi, piprodial);

    // now the individual tunings 140->146    
    if(dial>=141) rw.Systematics().SetTwkDial(t2krew::kNIWG2014a_SF_RFG,1); // SF->RFG tuning
    else if(dial==140) {
      rw.Systematics().SetTwkDial(t2krew::kNIWG2014a_SF_RFG,0);
      rw.Systematics().SetTwkDial(t2krew::kNXSec_VecFFCCQE,402); 
    }

    if(dial==141){ // no RPA
      rw.Systematics().SetTwkDial(t2krew::kNIWG_rpaCCQE_norm,0);
      rw.Systematics().SetTwkDial(t2krew::kNIWG_rpaCCQE_shape,0);
    }
    else if(dial==142){ //  non-relativistic RPA
      rw.Systematics().SetTwkDial(t2krew::kNIWG_rpaCCQE_norm,1);
      rw.Systematics().SetTwkDial(t2krew::kNIWG_rpaCCQE_shape,-1);
    }
    else if(dial>=143){ //  relativistic RPA
      rw.Systematics().SetTwkDial(t2krew::kNIWG_rpaCCQE_norm,1);
      rw.Systematics().SetTwkDial(t2krew::kNIWG_rpaCCQE_shape,0);
    }
    
    rw.Reconfigure();
    storer->NewSystSet(rw.Systematics()); // save the current twk dial values to the weight storer class


    // Loop over SK entries and calculate weights.
    if(fNskEvts < 0) fNskEvts = skVtxs->GetEntries();
    //fNskEvts=1;
    for(int i = 0; i < fNskEvts; i++){
      skVtxs->GetEntry(i);
      Double_t weight =1.0;
      cout<<"interaction mode="<<skVtxs->mode<<endl;
      if(skVtxs->mode>0){
	weight = rw.CalcWeight(skVtxs);//the weight is changed for each variation.
	int err=(int) dial/7;
	int sig=dial%7;
	/*     if(err>=2 && err<7 && sig>1 && sig<5){

	cout<<"*********************************************************************************************************************************************************************************************************************************************************************"<<endl;
	cout<<"abs="<<absdial<<", cxlodial="<<cxlodial<<", cxhidial="<<cxhidial<<", inelodial="<<inelodial<<", inehidial="<<inehidial<<endl;
	cout<<"Param="<<err<<", sig="<<sig<<", HERE IS THE WEIGHT="<<weight<<endl;
	}*/
      }
      storer->AddWeight(weight); // add weight for each
      cout<<"reweight ="<<weight<<endl;
    } // event loop

  } // index of tweak dial changes


  cout << "Outputfile"<< fOutFileName <<" has weight storer tree"  << endl;


  storer->SaveToFile(); // save the weights to a file
#endif // __T2KRW_NIWG_ENABLED__
#endif // __T2KRW_NEUT_ENABLED__
  delete storer; 

 
  return 0;
}

// Print the cmd line syntax
void Usage(){
  cout << "Cmd line syntax should be:" << endl;
  cout << "generateWeightsFromSK_NIWGexample.exe -s sk_inputfile [-e nevents]" << endl;
}
/*
  int ParseArgs(int argc, char **argv){

  while( (argc > 1) && (argv[1][0] == '-') ){

  switch(argv[1][1]){
  case 'i':
  fSKFileName    = argv[2];
  ++argv; --argc;
  if (argc <= 2) break;
  else if ( argv[2][0] != '-' ){
  fNskEvts = atoi(argv[2]);
  ++argv; --argc;
  }
  break;
  case 'o':
  fOutFileName    = argv[2];
  ++argv; --argc;
  break;
  }
  ++argv; --argc;
  }
  return 0;
  }
*/
