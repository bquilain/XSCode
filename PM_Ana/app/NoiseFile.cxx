#define DEBUG
//#define DEBUG_HIT

#include <iostream>

#include "TFile.h"
#include "TTree.h"

#include "INGRIDEVENTSUMMARY.h"
#include "IngridBasicReconSummary.h"
#include "IngridHitSummary.h"

#include "IngridConstants.h"
#include "RunInfo.h"
#include "RunInfo.cxx"
#include "noiseana.h"

int main(int argc, char *argv[]){
  int result;
  char input[200];
  char output[200];
  int run = 0;
  int i_bunch = bunch1st_cyc;
  int ncyc = nCyc;
  vector<int> reconmod;
  vector<int> reconcyc;

  while((result=getopt(argc,argv,"hr:i:o:"))!=-1){
    switch(result){
    case 'r':
      run = atoi(optarg);
      break;
    case 'i':
      sprintf(input,"%s",optarg);
      break;
    case 'o':
      sprintf(output,"%s",optarg);
      break;
    case 'h':
      std::cout << "r: MR run #" << std::endl;
      std::cout << "i: input root file" << std::endl;
      std::cout << "o: output root file" << std::endl;
      exit(1);
    }
  }
  RunInfo *ri = new RunInfo();
  int num_bunch = (int)ri->GetNBunches(run);
  int nmod      = (int)ri->GetNModules(run);
  int f_bunch   = i_bunch + num_bunch;

  // Input File ########################################################
  TFile *fIn = new TFile(input,"READ");
  if(fIn->IsZombie()){
    cout << "Error opening file" << endl;
    exit(1);
  }
  TTree *tIn = (TTree*)fIn->Get("tree");
  IngridEventSummary *evtIn = new IngridEventSummary();
  tIn->SetBranchAddress("fDefaultReco.", &evtIn);
  const int nevt = (int)tIn->GetEntries();
#ifdef DEBUG
  std::cout << input << std::endl;
  std::cout << "nevt = " << nevt << std::endl;
#endif
  //if(nevt==0){
  //std::cout << "NO EVENTS >> end" << std::endl;
  //return -1;
  //}

  // Output File #########################################################
  TFile *fOut = new TFile(output,"RECREATE");
  TTree *tOut = new TTree("tree","");
  noiseana *noise = new noiseana();
  tOut->Branch("evtid",&noise->evtid,"evtid/I");
  tOut->Branch("stime",&noise->stime,"stime/I");
  tOut->Branch("nhit",&noise->nhit,"nhit/I");
  tOut->Branch("ncyc",noise->ncyc,"ncyc[16]/I");
  tOut->Branch("mod",noise->mod,"mod[nhit]/I");
  tOut->Branch("cyc",noise->cyc,"cyc[nhit]/I");
  tOut->Branch("view",noise->view,"view[nhit]/I");
  tOut->Branch("pln",noise->pln,"pln[nhit]/I");
  tOut->Branch("ch",noise->ch,"ch[nhit]/I");
  tOut->Branch("pe",noise->pe,"pe[nhit]/F");
  tOut->Branch("lope",noise->lope,"lope[nhit]/F");
  tOut->Branch("time",noise->time,"time[nhit]/F");
  //tOut->Branch("nhitmod",noise->nhitmod,"nhitmod[16]/I");
  //tOut->Branch("nhitpln",noise->nhitpln,"nhitpln[16][2][15]/I");
  //tOut->Branch("nhitch",noise->nhitch,"nhitch[16][2][15][24]/I");
#ifdef DEBUG
      std::cout << output << " is recreated." << std::endl;
#endif

  // Main ################################################################
  for(int ievt=0; ievt<nevt; ievt++){
    // ### start get entry ###############################################
#ifdef DEBUG
    if(ievt%1000==0) std::cout << "++++++++++++++++ " << ievt << " ++++++++++++++++" << std::endl;;
#endif
    tIn->GetEntry(ievt);

    noise->evtid = (int)ievt;
    noise->stime = (int)evtIn->time;

    int nrecon = (int)evtIn->NIngridBasicRecons();

    // ### Fill module# and cyc# of events whitch have BasicRecon ##########################################
    reconmod.clear();
    reconcyc.clear();
    for(int irecon=0; irecon<nrecon; irecon++){
      IngridBasicReconSummary *recon = (IngridBasicReconSummary*)evtIn->GetBasicRecon(irecon);
      if((recon->hitcyc<i_bunch)||(f_bunch<=recon->hitcyc)){
	reconmod.push_back(recon->hitmod);
	reconcyc.push_back(recon->hitcyc);
#ifdef DEBUG_RECON
	printf("ievt: %3d| cyc: %3d| mod: %3d <-fill\n",ievt,recon->hitcyc,recon->hitmod);
#endif
      }///### if cycle ###
    }///### for irecon ###
#ifdef DEBUG_RECON
    if((reconcyc.size()!=0)&&(reconmod.size()!=0))
      std::cout << reconcyc.size() << ", " << reconmod.size() << std::endl;
    for(int j=0; j<reconcyc.size(); j++){
      printf("ievt: %3d| cyc: %3d| mod: %3d\n",ievt,reconcyc[j],reconmod[j]);
    }
#endif

    // fill noise info ###############################
    int allhit = 0;
    for(int imod=0; imod<nmod; imod++){
      // init ++++++++++++++++
      int nhitmod = 0;
      int nhitpln[2][15];
      int nhitch[2][15][24];
      for(int i=0; i<15; i++){
	for(int j=0; j<24; j++){
	  nhitch[0][i][j] = 0;
	  nhitch[1][i][j] = 0;
	}
	nhitpln[0][i] = 0;
	nhitpln[1][i] = 0;
      }
      
      //
      int noise_cyc = ncyc - num_bunch; //for debug
      for(int icyc=0; icyc<ncyc; icyc++){
	if((i_bunch<=icyc)&&(icyc<f_bunch)) continue;
	bool havehit = false; //for debug
	for(unsigned int j=0; j<reconcyc.size(); j++){
	  if((reconmod[j]==imod)&&(reconcyc[j]==icyc)){	
#ifdef DEBUG_RECON
	    printf("ievt: %3d| cyc: %3d| mod: %3d <-rejected\n",ievt,reconcyc[j],reconmod[j]);
#endif
	    havehit = true; //for debug
	  }
	}// ### cyc ###
	if(havehit){ //for debug
	  noise_cyc--;
#ifdef DEBUG_CYCLE
	  std::cout << "# of cycles whitch have noise = " << noise_cyc << std::endl; 
#endif
	  continue;
	}
	int nhit = (int)evtIn->NIngridModHits(imod, icyc);
	nhitmod += nhit;
#ifdef DEBUG_HIT
	std::cout << "---------------------------" << std::endl;
	std::cout << "imod = " << imod << std::endl;
	std::cout << "icyc = " << icyc << std::endl;
	std::cout << "nhit = " << nhit << std::endl;
#endif
	for(int ihit=0; ihit<nhit; ihit++){
	  IngridHitSummary *hit = (IngridHitSummary*)evtIn->GetIngridModHit(ihit, imod, icyc);
	  int view = (int)hit->view;
	  int  pln = (int)hit->pln;
	  int   ch = (int)hit->ch;
	  noise->mod[allhit] = (int)hit->mod;
	  noise->cyc[allhit] = (int)hit->cyc;
	  noise->view[allhit] = view;
	  noise->pln[allhit] = pln;
	  noise->ch[allhit] = ch;
	  noise->pe[allhit] = (float)hit->pe;
	  noise->lope[allhit] = (float)hit->lope;
	  noise->time[allhit] = (float)hit->time;
	  allhit++;
	  nhitpln[view][pln]++;
	  nhitch[view][pln][ch]++;
	}// ### ihit ###
      }// ### icyc ###
      // fill +++++++++++++++++++++++++++
      noise->ncyc[imod] = (int)noise_cyc;
      //noise->nhitmod[imod] = nhitmod;
      for(int i=0; i<15; i++){
	for(int j=0; j<24; j++){
	  //noise->nhitch[imod][0][i][j] = nhitch[0][i][j];
	  //noise->nhitch[imod][1][i][j] = nhitch[1][i][j];
	}
	//noise->nhitpln[imod][0][i] = nhitpln[0][i];
	//noise->nhitpln[imod][1][i] = nhitpln[1][i];
      }
    }// ### imod ###
    noise->nhit = allhit;
    tOut->Fill();
  }// ### ievt ###
  fOut->Write();
  fIn->Close();
  fOut->Close();
  delete ri;
  delete fIn;
  delete fOut;

  return 0;
}
