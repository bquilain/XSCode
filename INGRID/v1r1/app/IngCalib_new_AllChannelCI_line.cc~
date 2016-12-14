//##### Standard C++ lib. ######
#include<iostream>
#include<sstream>
#include<fstream>
#include <iomanip>
#include <sys/stat.h>
using namespace std;
//##### Root Library ###########
#include <TROOT.h>
#include <TStyle.h>
#include <TApplication.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TObject.h>
#include <TEventList.h>
#include <TBranch.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
//##### INGRID Library #########
#include "INGRIDEVENTSUMMARY.h"
#include "IngridHitSummary.h"
#include "IngridSimHitSummary.h"
//#define DEBUG

const int nMod  = 17;
const int nView =  2;
const int nPln  = 24;
const int nCh   = 32;
char hgraph[100], lgraph[100];
TFile * CI_File=NULL;

const int online_calib_run = 13047;
float     High_ADC_Cal[ nMod ][ nView ][ nPln ][ nCh ];
float     Low_ADC_Cal[ nMod ][ nView ][ nPln ][ nCh ];
float     gain[ nMod ][ nView ][ nPln ][ nCh ];
float     ped[ nMod ][ nView ][ nPln ][ nCh ];
float     logain[ nMod ][ nView ][ nPln ][ nCh ];
float     loped[ nMod ][ nView ][ nPln ][ nCh ];
float     gain_correction[ nMod ][ nView ][ nPln ][ nCh ];
float     Gain_Cal[ nMod ][ nView ][ nPln ][ nCh ];
float     logain_correction[ nMod ][ nView ][ nPln ][ nCh ];
float     Gain_Cal_Lo[ nMod ][ nView ][ nPln ][ nCh ];
float     gain_corrected[ nMod ][ nView ][ nPln ][ nCh ];
float     logain_corrected[ nMod ][ nView ][ nPln ][ nCh ];
float     LoGain[ nMod ][ nView ][ nPln ][ nCh ];
float value[ nMod ][ nView ][ nPln ][ nCh ];
float gain_line[ nMod ][ nView ][ nPln ][ nCh ][100];
int NPoints[ nMod ][ nView ][ nPln ][ nCh ][100];
float NHiADC[ nMod ][ nView ][ nPln ][ nCh ][100];
float NLoADC[ nMod ][ nView ][ nPln ][ nCh ][100];
float NPEHi[ nMod ][ nView ][ nPln ][ nCh ][100];
float logain_line[ nMod ][ nView ][ nPln ][ nCh ][100];
float logain_calib[ nMod ][ nView ][ nPln ][ nCh ];

TGraph *Hi_Graph;
TGraph *Lo_Graph;
TF1 * Hi_Fit[ nMod ][ nView ][ nPln ][ nCh ];
TF1 * Lo_Fit[ nMod ][ nView ][ nPln ][ nCh ];

void InitilizeCalibConst(){
  for(int mod=0; mod<nMod; mod++){
    for(int view=0; view<nView; view++){
      for(int pln=0; pln<nPln; pln++){
	for(int ch=0; ch<nCh; ch++){
	  gain    [mod][view][pln][ch] =  1e5;
	  logain  [mod][view][pln][ch] =  1e5;
	  ped     [mod][view][pln][ch] = -1e-5;
	  loped   [mod][view][pln][ch] = -1e-5;
	  /*###############added by B.Quilain#######################*/
	  Hi_Fit[mod][view][pln][ch]= new TF1("Hi_Fit","0",0,10000);
	  Lo_Fit[mod][view][pln][ch]= new TF1("Lo_Fit","0",0,10000);
	  LoGain   [mod][view][pln][ch] = -1e-5;
	  value [mod][view][pln][ch] = -1e-5;
	  logain_calib  [mod][view][pln][ch] =  1e5;
	  for(int i=0;i<100;i++){      
	    gain_line[mod][view][pln][ch][i] =0; 
	    NPoints[mod][view][pln][ch][i]=0;
	    NHiADC[mod][view][pln][ch][i]=0;
	    NLoADC[mod][view][pln][ch][i]=0;
	    NPEHi[mod][view][pln][ch][i]=0;
	    logain_line[mod][view][pln][ch][i]=0;
	  }

	}
      }
    }
  }
}


float HighADC[ nMod ][ nView ][ nPln ][ nCh ];
float HighGain[ nMod ][ nView ][ nPln ][ nCh ];
float HighPed [ nMod ][ nView ][ nPln ][ nCh ];
float LowADC[ nMod ][ nView ][ nPln ][ nCh ];
int   LowPed  [ nMod ][ nView ][ nPln ][ nCh ];
bool  GoodCh  [ nMod ][ nView ][ nPln ][ nCh ];
int   NDChan  [ nMod ][ nView ][ nPln ][ nCh ];

bool GetCalibConst(int irun, int isubrun){
  TFile* fTCalib   = new TFile( Form("/home/bquilain/Ingrid_Process/calib_root_file/ingrid_%08d_%04d_Calib00.root", irun, isubrun),"read" );
  TTree* calibtree = (TTree*)fTCalib->Get("calibtree");

  calibtree       -> SetBranchAddress("HighADC", HighADC);
  calibtree       -> SetBranchAddress("HighGain", HighGain);
  calibtree       -> SetBranchAddress("HighPed" , HighPed );
  calibtree       -> SetBranchAddress("LowADC"   , LowADC);
  calibtree       -> SetBranchAddress("LowPed"  , LowPed  );
  calibtree       -> SetBranchAddress("GoodCh"  , GoodCh  );
  calibtree       -> SetBranchAddress("NDChan"  , NDChan);
  calibtree       -> GetEntry(0);//vector are filled

  for(int imod=0; imod<nMod; imod++){
    if(imod==16) continue;//skip the PM
    for(int iview=0; iview<nView; iview++){
      for(int ipln=0; ipln<nPln; ipln++){
	for(int ich=0; ich<nCh ; ich++){
	  if( GoodCh[imod][iview][ipln][ich] ){
	    //so here, link channel,plan,view... with hte #channel in ND280 calib
	    //récupérer sa mean ADC en plus!
	    //?
	    High_ADC_Cal [imod][iview][ipln][ich] = HighADC[imod][iview][ipln][ich]; 
	    gain  [imod][iview][ipln][ich] = HighGain[imod][iview][ipln][ich];
	    logain[imod][iview][ipln][ich] = 0.1 * HighGain[imod][iview][ipln][ich];
	    Low_ADC_Cal [imod][iview][ipln][ich] = LowADC[imod][iview][ipln][ich];
	    ped   [imod][iview][ipln][ich] = HighPed[imod][iview][ipln][ich];
	    loped [imod][iview][ipln][ich] = LowPed[imod][iview][ipln][ich];
	    
	  }
	  else{
	    High_ADC_Cal [imod][iview][ipln][ich] = HighADC[imod][iview][ipln][ich]; 
	    Low_ADC_Cal [imod][iview][ipln][ich] = LowADC[imod][iview][ipln][ich];
	    gain  [imod][iview][ipln][ich] = 1000;
	    logain[imod][iview][ipln][ich] = 1000;
	    ped   [imod][iview][ipln][ich] = 1000;
	    loped [imod][iview][ipln][ich] = 1000;
	  }

	}
      }
    }
  }
    
  //### double use VETO 
  for(int i=1; i<=6; i++){
    for(int j=0; j<22; j++){
      double temp;
      gain  [i][1][11][j] = gain  [i-1][1][12][j];
      ped   [i][1][11][j] = ped   [i-1][1][12][j];
      logain[i][1][11][j] = logain[i-1][1][12][j];
      loped [i][1][11][j] = loped [i-1][1][12][j];
    }
  }
  for(int i=8; i<=13; i++){
    for(int j=0; j<22; j++){
      gain[i][0][13][j]   = gain[i-1][0][14][j];
      ped [i][0][13][j]   = ped [i-1][0][14][j];
      logain[i][0][13][j] = gain[i-1][0][14][j];
      loped [i][0][13][j] = ped [i-1][0][14][j];
    }
  }
  return true;
}


void GetLinearization(){
  /*TFile * CI_File=new TFile("/home/bquilain/Ingrid_Process/CISimpleOutFile.root");
  char Hi_Graph_Adress[100], Lo_Graph_Adress[100];
  int nchantot=2;
  for(int i=1;i<nchantot;i++){
    
    sprintf(Hi_Graph_Adress,"hiGraph_%d",i);
    TGraph *Hi_Graph=(TGraph*) CI_File->Get(Hi_Graph_Adress);
    sprintf(Lo_Graph_Adress,"loGraph_%d",i);
    TGraph *Lo_Graph=(TGraph*) CI_File->Get(Lo_Graph_Adress);
    }*/
  
 for(int imod=0; imod<nMod; imod++){
    for(int iview=0; iview<nView; iview++){
      for(int ipln=0; ipln<nPln; ipln++){
	for(int ich=0; ich<nCh ; ich++){
	  if( GoodCh[imod][iview][ipln][ich] ){
#ifdef DEBUG
	    //	    cout<<"new iteration = Module:"<<imod<<"  , Plane:"<<ipln<<"  , Channel:"<<ich<<endl;
#endif
	    //placer le highadccalib_value sur le graph puis corriger le gain en regardant ou est le highadcvalue (pour cela, on va surement
	    // obtenir la CI équivalente a High_ADC_Cal [imod][iview][ipln][ich]
	    //
	//Now we have to know what is the adc counts value when the charge was calibrated!? d'ou ça vient surtout?
	//Prendre la mean de la distri! Attention si!!! le pedestal a été enlevé dans tfbCalib déja, faut pas que ce soit le cas. Sinon prendre mean-pedestal
	//l'autre option, est de refitter pedestal et 1er pic (comparer les valeurs avec celles prcédentes dans chaque histo, ais comment savoir quel histo est associé à quelle channel?
	    //High_ADC_Cal[imod][iview][ipln][ich]
	    //en sortir la charge injectée équivalente: diviser l'un par l'autre
	    //faire de meme avec le #adc du run et sa charge injectée. Diviser l'un par l'autre=> gain a modifier par le rapport des deux
	    /*High_ADC_Graph = Hi_Graph->GetX();
	    High_CI_Graph = Hi_Graph->GetY();
	    Low_ADC_Graph = Lo_Graph->GetX();
	    Low_CI_Graph = Lo_Graph->GetY();*/
	    //en fait faut récupérer la fonction de fit, pas le graph
	    
	    
	    sprintf(hgraph,"hiGraph_%d",NDChan[imod][iview][ipln][ich]);
	    sprintf(lgraph,"loGraph_%d",NDChan[imod][iview][ipln][ich]);
#ifdef DEBUG
	    //   cout<<"High Graph = "<<hgraph<<"  ,Low Graph =  "<<lgraph<<endl;
#endif
	    Hi_Graph=(TGraph*) CI_File->Get(hgraph);
	    Lo_Graph=(TGraph*) CI_File->Get(lgraph);
	    
	    Hi_Fit[imod][iview][ipln][ich]=(TF1*)Hi_Graph->GetFunction("fFitHi");
	    Lo_Fit[imod][iview][ipln][ich]=(TF1*)Lo_Graph->GetFunction("fFitLo");
	    //cout<<"Lo_Fit is evaluated"<<endl;
	    //cout<<High_ADC_Cal[imod][iview][ipln][ich]<<"   "<<Low_ADC_Cal[imod][iview][ipln][ich]<<endl;
	    double Hi_CI_Cal=Hi_Fit[imod][iview][ipln][ich]->Eval(gain[imod][iview][ipln][ich]);/*High_ADC_Cal[imod][iview][ipln][ich]-ped[imod][iview][ipln][ich]*/
	    double Lo_CI_Cal=Lo_Fit[imod][iview][ipln][ich]->Eval(Low_ADC_Cal[imod][iview][ipln][ich]-loped[imod][iview][ipln][ich]);
#ifdef DEBUG
	    //cout<<"At INGRID CAlib Point, ADC in High Gain Ch is = "<< High_ADC_Cal[imod][iview][ipln][ich] - ped[imod][iview][ipln][ich]<<" , and in Low Gain Ch, ADC = "<<Low_ADC_Cal[imod][iview][ipln][ich] - loped[imod][iview][ipln][ich]<<endl;
	    //cout<<"At INGRID CAlib Point, equivalent CI in High Gain Channel is = "<<Hi_CI_Cal<<" , and Low Gain CI is = "<<Lo_CI_Cal<<endl;
#endif	
	    //double Hi_CI=Hi_Fit(High_ADC[imod][iview][ipln][ich]);
	    //double Lo_CI=Lo_Fit(Low_ADC[imod][iview][ipln][ich]);
	    
	    Gain_Cal[imod][iview][ipln][ich]=Hi_CI_Cal/gain[imod][iview][ipln][ich];/*/(High_ADC_Cal[imod][iview][ipln][ich]-ped[imod][iview][ipln][ich])*/
	    //Gain_Cal_Lo[imod][iview][ipln][ich]=Lo_CI_Cal/(Low_ADC_Cal[imod][iview][ipln][ich] - loped[imod][iview][ipln][ich]);
	    //LoGain[imod][iview][ipln][ich]=(High_ADC_Cal[imod][iview][ipln][ich]-ped[imod][iview][ipln][ich])*gain[imod][iview][ipln][ich]/(Low_ADC_Cal[imod][iview][ipln][ich] - loped[imod][iview][ipln][ich]);
#ifdef DEBUG
	    //cout<<"Low Gain Original Value = "<<logain[imod][iview][ipln][ich]<<" ,  Low Gain Value = "<<LoGain[imod][iview][ipln][ich]<<endl;
	    //cout<<Gain_Cal[imod][iview][ipln][ich]<<"   "<<Gain_Cal_Lo[imod][iview][ipln][ich]<<endl<<endl;
#endif
	    //double Gain=Hi_CI/High_ADC[imod][iview][ipln][ich];
	    //gain_correction[imod][iview][ipln][ich] =Gain/Gain_Cal;
	    //gain[imod][iview][ipln][ich]*= gain_correction[imod][iview][ipln][ich];
	  }
	}
      }
    }
 }
}

const static Int_t   cTimeCorrectionBase = 24;//Time correction for difference
                                          //of cable length

int main(int argc,char *argv[]){
  TROOT        root ("GUI","GUI");
  TApplication app  ("App",0,0);
  int run_number, sub_run_number;
  int c=-1;
  FileStat_t fs;
  bool cosmic      = false;
  bool PModule     = false;
  bool semioffline = false; //if true, use the reference table
  char Output[300], FileName[300];
  int chantot=0;
  TH1D * OldGain = new TH1D("OldGain","Gain Value",100,0,30);
  TH1D * NewGain = new TH1D("NewGain","New Gain Value",100,0,30);
  TH1D * GNL     = new TH1D("GNL","Gain non Linearity",200,0,1000);
  TH1D * Evt     = new TH1D("Evt","Evt by ADC",100,0,100);
  TH1D * Evt2     = new TH1D("Evt2","Evt by ADC",100,0,100);
  TH1D * Evt0     = new TH1D("Evt0","Evt by ADC",100,0,100);
  TH1D * Evt20     = new TH1D("Evt20","Evt by ADC",100,0,100);
  TH1D * LoHi0     = new TH1D("LoHi0","Lo/Hi p.e ",100,0,100);
  TH1D * LoHi     = new TH1D("LoHi","Lo/Hi p.e ",100,0,100);
  TH1D * Lo     = new TH1D("Lo","Lo p.e ",100,0,100);
  TH1D * Lo0     = new TH1D("Lo0","Lo p.e ",100,0,100);
  TH2D * HiLine     = new TH2D("HiLine","Lo/Hi p.e ",100,0,30,100,0,30);
  TH2D * LoHipe     = new TH2D("LoHipe","Lo/Hi p.e ",100,0,100,100,0,100);
  TH2D * LoHipeStd     = new TH2D("LoHipeStd","Lo/Hi p.e ",100,0,100,100,0,100);
  double CI_ADC=0;
  //float NPoints[100];
  //for(int i=0;i<100;i++){
  //NPoints[i]=0;
  //    }

  while ((c = getopt(argc, argv, "r:s:cqpf:o:")) != -1) {
    switch(c){
    case 'f':
      sprintf(FileName,"%s",optarg);
      break;
    case 'o':
      sprintf(Output,"%s",optarg);
      break;
    case 'r':
      run_number=atoi(optarg);
      break;
    case 's':
      sub_run_number=atoi(optarg);
      break;
    case 'c':
      cosmic = true;
      break;
    case 'q':
      semioffline = true;
      break;
    case 'p':
      PModule = true;
      break;
    }
  }
  cout<<"######################### Careful, the PM is not implemented yet ###########################"<<endl;
  InitilizeCalibConst(); //### Initialize gain, pedestal, and so on


  //#### Read file before calibration ######
  //########################################
  IngridEventSummary* summary = new IngridEventSummary();
  /*
  if(!cosmic&&!PModule)
    sprintf(buff,"/home/daq/data/dst/ingrid_%08d_%04d.root",
	    run_number, sub_run_number);
  if(cosmic&&!PModule)
    sprintf(buff,"/home/daq/data/cosmic/ingrid_%08d_%04d.root",
	    run_number, sub_run_number);
  if(!cosmic&&PModule)
    sprintf(buff,"/home/daq/data/PM/ingrid_%08d_%04d.root",
	    run_number, sub_run_number);
  if(cosmic&&PModule)
    sprintf(buff,"/home/daq/data/PM_cosmic/ingrid_%08d_%04d.root",
	    run_number, sub_run_number);
  */
  if(gSystem->GetPathInfo(FileName,fs)){
    cout<<"Cannot open file "<< FileName <<endl;
    exit(1);
  }

  /*
  std::cout << "Calib Read file number: " << run_number 
	    << "_" << sub_run_number << std::endl;
  */
  TFile*            rfile = new TFile(FileName,"read");
  TTree*             tree = (TTree*)rfile -> Get("tree");
  TBranch*          EvtBr = tree->GetBranch("fDefaultReco.");
  EvtBr                   ->SetAddress(&summary);
  tree                    ->SetBranchAddress("fDefaultReco.", &summary);
  int                nevt = (int)tree -> GetEntries();
  cout << "Total # of events = " << nevt <<endl;

  if(!semioffline)
    GetCalibConst(run_number,sub_run_number);
  else if(semioffline)
    GetCalibConst(online_calib_run,0);

  //#### Create file after calibration ######
  //#########################################
  /*
  if(!cosmic&&!PModule)
    sprintf(buff,"/home/daq/data/dst/ingrid_%08d_%04d_calib.root",
	    run_number, sub_run_number);
  if(cosmic&&!PModule)
    sprintf(buff,"/home/daq/data/cosmic/ingrid_%08d_%04d_calib.root",
	    run_number, sub_run_number);
  if(!cosmic&&PModule)
    sprintf(buff,"/home/daq/data/PM/ingrid_%08d_%04d_calib.root",
	    run_number, sub_run_number);
  if(cosmic&&PModule)
    sprintf(buff,"/home/daq/data/PM_cosmic/ingrid_%08d_%04d_calib.root",
	    run_number, sub_run_number);
  */
#ifdef DEBUG
  //cout<<" Calib non linearized is done"<<endl;
#endif
  CI_File=new TFile("/home/bquilain/Ingrid_Process/CISimpleOutFile.root");

  Hi_Graph=(TGraph*) CI_File->Get("hiGraph_0");
  Lo_Graph=(TGraph*) CI_File->Get("loGraph_0");
  
  GetLinearization();
#ifdef DEBUG
  //cout<<"Calib linearized is done"<<endl;
#endif
  TFile*            wfile = new TFile(Output,"recreate");
  TTree*            wtree = new TTree("tree","tree");
  IngridEventSummary* wsummary = new IngridEventSummary(); 
  wtree              -> Branch   ("fDefaultReco.","IngridEventSummary", 
				 &wsummary,  64000,   99);
  for(int ievt = 0; ievt < nevt; ++ievt){
    if(ievt%1000==0)cout<<"Run#"<<run_number<<"\tcalib event:"<<ievt<<endl;
    summary -> Clear();
    wsummary-> Clear();
    tree    -> GetEntry(ievt);
    for(int mod=0; mod<nMod; mod++){
      if(mod==16) continue; //for now, skip the PM                                                                                                     
      for(int cyc=0; cyc<23; cyc++){
	int ninghit = summary -> NIngridModHits(mod, cyc);
#ifdef DEBUG
	//cout<<endl;                                                                                                                                  
	//cout<<"Module:"<<mod<<"   , Cycle:"<<cyc<<endl;                                                                                              
#endif
	    for(int i=0; i<ninghit; i++){

	      IngridHitSummary *inghitsum;
	      inghitsum   = (IngridHitSummary*) (summary -> GetIngridModHit(i, mod, cyc) );

	      int view = inghitsum -> view;
	      int pln  = inghitsum -> pln;
	      int ch   = inghitsum -> ch;
#ifdef DEBUG
	      //  cout<<"Plan = "<<pln<<" , Channel = "<<ch<<" , View = "<<view<<endl;
#endif
	      CI_ADC=Hi_Fit[mod][view][pln][ch]->Eval(inghitsum->adc-ped[mod][view][pln][ch]);
	      gain_correction[mod][view][pln][ch]=Gain_Cal[mod][view][pln][ch]/(CI_ADC/(inghitsum->adc-ped[mod][view][pln][ch]));//R2/R1            
	      gain_corrected[mod][view][pln][ch]=gain[mod][view][pln][ch]*gain_correction[mod][view][pln][ch];
	      int line = (int) (inghitsum->adc-ped[mod][view][pln][ch])/10;
	      gain_line[mod][view][pln][ch][line] +=gain_corrected[mod][view][pln][ch];
	      NPoints[mod][view][pln][ch][line]++;
	      NHiADC[mod][view][pln][ch][line]+=inghitsum->adc-ped[mod][view][pln][ch];
	      NLoADC[mod][view][pln][ch][line]+=inghitsum->loadc-loped[mod][view][pln][ch];
	      NPEHi[mod][view][pln][ch][line]+=(inghitsum->adc-ped[mod][view][pln][ch])/gain_corrected[mod][view][pln][ch];

	      if((inghitsum->adc-ped[mod][view][pln][ch]<80) && (inghitsum->adc-ped[mod][view][pln][ch])>value[mod][view][pln][ch]){
		LoGain[mod][view][pln][ch]=gain_corrected[mod][view][pln][ch]*(inghitsum->loadc-loped[mod][view][pln][ch])/(inghitsum->adc-ped[mod][view][pln][ch]);
		value[mod][view][pln][ch]=inghitsum->adc-ped[mod][view][pln][ch];
		//toeactiveifpointGain_Cal_Lo[mod][view][pln][ch]=Lo_Fit[mod][view][pln][ch]->Eval(inghitsum->loadc-loped[mod][view][pln][ch])/(inghitsum->loadc-loped[mod][view][pln][ch]);
	      }
	    }
	    //	    for(int ilin=0;ilin<100;ilin++){
	    // gain_line[mod][view][pln][ch][ilin]/=NPoints[ilin];
	    // cout<<ilin<<"  "<<gain_line[mod][view][pln][ch][ilin]<<endl;
	    //}
      }
    }
  }
  
  int Cal=0;
  for(int mod=0;mod<nMod;mod++){
  for(int view=0;view<2;view++){
    for(int pln=0;pln<nPln;pln++){
      for(int ch=0;ch<nCh;ch++){
	for(int ilin=0;ilin<100;ilin++){
	  if(NPoints[mod][view][pln][ch][ilin]==0) continue;
	  else {
	    gain_line[mod][view][pln][ch][ilin] /= NPoints[mod][view][pln][ch][ilin];
	    NHiADC[mod][view][pln][ch][ilin] /= NPoints[mod][view][pln][ch][ilin];
	    NLoADC[mod][view][pln][ch][ilin] /= NPoints[mod][view][pln][ch][ilin];
	    NPEHi[mod][view][pln][ch][ilin] /= NPoints[mod][view][pln][ch][ilin];
	    logain_line[mod][view][pln][ch][ilin]=NLoADC[mod][view][pln][ch][ilin]/(NHiADC[mod][view][pln][ch][ilin]/gain_line[mod][view][pln][ch][ilin]); //NPEHi[mod][view][pln][ch][ilin];
	      }
#ifdef DEBUG
	  //if(mod==3){cout<<"Module = "<<mod<<" , Plane = "<<pln<<" ,Channel = "<<ch<<" , View = "<<view<<endl;
	  //  cout<<"ADC Range = "<<ilin*10<<" , Gain = "<<gain_line[mod][view][pln][ch][ilin]<<"  , LoGain = "<< logain_line[mod][view][pln][ch][ilin] <<endl<<endl;
	  //  }
#endif
      }
	if(logain_line[mod][view][pln][ch][7]!=0) Cal=7;
	else if(logain_line[mod][view][pln][ch][8]!=0) Cal=8;
	else if(logain_line[mod][view][pln][ch][6]!=0) Cal=6;
	else if(logain_line[mod][view][pln][ch][5]!=0) Cal=5;
	else if(logain_line[mod][view][pln][ch][4]!=0) Cal=4;
        else if(logain_line[mod][view][pln][ch][3]!=0) Cal=3;
        else if(logain_line[mod][view][pln][ch][2]!=0) Cal=2;
	else{
	  logain_calib[mod][view][pln][ch]=1;
	  //TOREACTIVE IF LINE
	  Gain_Cal_Lo[mod][view][pln][ch]=1;
	  continue;
	}
	logain_calib[mod][view][pln][ch]=logain_line[mod][view][pln][ch][Cal];
	//TOREACTIVE IF LINE
	Gain_Cal_Lo[mod][view][pln][ch]=Lo_Fit[mod][view][pln][ch]->Eval(NLoADC[mod][view][pln][ch][Cal])/(NLoADC[mod][view][pln][ch][Cal]);
#ifdef DEBUG
	//cout<<NLoADC[mod][view][pln][ch][Cal]<<"   "<<loped[mod][view][pln][ch]<<endl;
	//cout<<"Logain Used for Calib = "<<logain_calib[mod][view][pln][ch]<<" , Slope of the CI plot = "<<Gain_Cal_Lo[mod][view][pln][ch]<<endl;
#endif
      }
    }
  }
 }

  for(int ievt = 0; ievt < nevt; ++ievt){
    if(ievt%1000==0)cout<<"Run#"<<run_number<<"\tcalib event:"<<ievt<<endl;
    summary -> Clear();
    wsummary-> Clear();
    tree    -> GetEntry(ievt);
    //la boucle ici est sur les modules, puis sur les cycles, puis suir les hits. Mais c'est bon, on veut le nombre d'ADC du hit, puis savoir dans quel module, plan, view, channel il est, puis corriger
    for(int mod=0; mod<nMod; mod++){
      if(mod==16) continue; //for now, skip the PM
      for(int cyc=0; cyc<23; cyc++){
        int ninghit = summary -> NIngridModHits(mod, cyc);
#ifdef DEBUG
	//cout<<endl;
	//cout<<"Module:"<<mod<<"   , Cycle:"<<cyc<<endl;
#endif
        for(int i=0; i<ninghit; i++){

	  IngridHitSummary *inghitsum;
	  inghitsum   = (IngridHitSummary*) (summary -> GetIngridModHit(i, mod, cyc) );
	  
	  int view = inghitsum -> view;
	  int pln  = inghitsum -> pln;
	  int ch   = inghitsum -> ch;
#ifdef DEBUG
	  //cout<<"Plan = "<<pln<<" , Channel = "<<ch<<" , View = "<<view<<endl;
#endif
	  //	  logain_corrected[mod][view][pln][ch]=LoGain[mod][view][pln][ch]*(Lo_Fit[mod][view][pln][ch]->Eval(inghitsum->loadc - loped[mod][view][pln][ch]))/(Lo_Fit[mod][view][pln][ch]->Eval(LoGain[mod][view][pln][ch]));
	  int LoRange=(int) (inghitsum->adc - ped[mod][view][pln][ch])/10;
	  //cout<<LoRange<<endl;
	  /*int Cal=0;

	  if(logain_line[mod][view][pln][ch][3]!=0) Cal=3;
	  else if(logain_line[mod][view][pln][ch][4]!=0) Cal=4;
	  else if(logain_line[mod][view][pln][ch][2]!=0) Cal=2;
	  else if(logain_line[mod][view][pln][ch][5]!=0) Cal=5;
	  else continue;
	  
	  Cal_ADC=logain_line[mod][view][pln][ch][Cal]
	  Gain_Cal_Lo[mod][view][pln][ch]=Lo_Fit[mod][view][pln][ch]->Eval(Cal)/(logain_line[mod][view][pln][ch][Cal]);*/
	  
	  logain_corrected[mod][view][pln][ch]=logain_calib[mod][view][pln][ch]*(Gain_Cal_Lo[mod][view][pln][ch])/(Lo_Fit[mod][view][pln][ch]->Eval(inghitsum->loadc - loped[mod][view][pln][ch])/(inghitsum->loadc - loped[mod][view][pln][ch]));

	  //gain_corrected[mod][view][pln][ch]=gain_line[mod][view][pln][ch][LoRange];
	  
	  //Gain_Cal_Lo[mod][view][pln][ch]=Lo_Fit[mod][view][pln][ch]->Eval(LoGain[mod][view][pln][ch])/(LoGain[mod][view][pln][ch]);//R1
	  //TOREACTIVEIFPOINTBYPOINT
	  //logain_corrected[mod][view][pln][ch]=LoGain[mod][view][pln][ch]*(Gain_Cal_Lo[mod][view][pln][ch])/(Lo_Fit[mod][view][pln][ch]->Eval(inghitsum->loadc - loped[mod][view][pln][ch])/(inghitsum->loadc - loped[mod][view][pln][ch]));
	  //####We had a First linearity correction (Jan.2013) by B.Quilain#######//
	  //charger la fonction: 
	  //double CI_ADC=Hi_Fit[mod][view][pln][ch]->Eval(inghitsum->adc-ped[mod][view][pln][ch]);
	  //cout<<CI_ADC<<endl;
#ifdef DEBUG
	  //cout<<"Default Low Gain was ="<<logain[mod][view][pln][ch]<<" , Low Gain is = "<<LoGain[mod][view][pln][ch]<<endl;
	  //if(inghitsum->adc>400) cout<<"High Gain was ="<<gain[mod][view][pln][ch]<<"  , Charge Injected is = "<<CI_ADC<<"  ,  ADC equivalent is = "<<inghitsum->adc<<endl;
	  //	  cout<<"High Gain was ="<<gain[mod][view][pln][ch]<<"  , Charge Injected is = "<<CI_ADC<<"  ,  ADC equivalent is = "<<inghitsum->adc<<endl;
#endif
	  //gain_correction[mod][view][pln][ch]=(CI_ADC/(inghitsum->adc-ped[mod][view][pln][ch]))/Gain_Cal[mod][view][pln][ch];//R2/R1
	  //gain_corrected[mod][view][pln][ch]=gain[mod][view][pln][ch]*gain_correction[mod][view][pln][ch];
	  //cout<<"High Gain corrected is now = "<<gain[mod][view][pln][ch]<<endl;
          //double Conversion_CIPE=(CI_ADC/(inghitsum->adc-ped[mod][view][pln][ch]))/gain_corrected[mod][view][pln][ch];
	  //cout<<"Conversion Factor = "<<Conversion_CIPE<<endl;
	  
	  double CI_ADC_Lo=Lo_Fit[mod][view][pln][ch]->Eval(inghitsum->loadc - loped[mod][view][pln][ch]);
          //logain[mod][view][pln][ch]=(1/Conversion_CIPE)*(CI_ADC_Lo/(inghitsum->loadc - loped[mod][view][pln][ch]));
																			  
 //logain_correction[mod][view][pln][ch]=(Lo_Fit[mod][view][pln][ch]->Eval(inghitsum->loadc - loped[mod][view][pln][ch])/(inghitsum->loadc - loped[mod][view][pln][ch]))/Gain_Cal_Lo[mod][view][pln][ch];
																			   //logain_corrected[mod][view][pln][ch]=logain[mod][view][pln][ch]*logain_correction[mod][view][pln][ch];
#ifdef DEBUG
	  //if(CI_ADC>200 || CI_ADC_Lo>200) {
	  //cout<<"High Gain was ="<<gain[mod][view][pln][ch]<<"  , Charge Injected is = "<<CI_ADC<<"  ,  ADC equivalent is = "<<inghitsum->adc - ped[mod][view][pln][ch]<<endl;
	  //cout<<"High Gain corrected is now = "<<gain_corrected[mod][view][pln][ch]<<endl;
	  //cout<<"Conversion Factor = "<<Conversion_CIPE<<endl;
	  //cout<<"Lo Gain was ="<<logain[mod][view][pln][ch]<<"  , Charge Injected is = "<<CI_ADC_Lo<<"  ,  ADC equivalent is = "<<inghitsum->loadc - loped[mod][view][pln][ch]<<endl;
	  //cout<<"LoGain corrected is now = "<<logain_corrected[mod][view][pln][ch]<<endl<<endl;
	  //if(TMath::Abs(logain_corrected[mod][view][pln][ch])<.5){
	  cout<<"Number of ADC in HG Channel = "<<inghitsum->adc - ped[mod][view][pln][ch]<<endl;
	  cout<<"Before Calibration, Amount of p.e is in HG channel = "<<(inghitsum->adc - ped[mod][view][pln][ch])/gain[mod][view][pln][ch]<<" ,  p.e in Lo Gain Channel = "<<(inghitsum->loadc - loped[mod][view][pln][ch])/logain[mod][view][pln][ch]<<endl;
	  cout<<"Amount of p.e is in HG channel = "<<(inghitsum->adc - ped[mod][view][pln][ch])/gain_corrected[mod][view][pln][ch]<<" ,  p.e in Lo Gain Channel = "<<(inghitsum->loadc - loped[mod][view][pln][ch])/logain_corrected[mod][view][pln][ch]<<endl;
	  cout<<"Logain at Calib Point = "<<logain_calib[mod][view][pln][ch]<<" , Slope of the CI at Calib = "<<Gain_Cal_Lo[mod][view][pln][ch]<<endl<<endl;
	  //cout<<"Logain at Calib Point = "<<LoGain[mod][view][pln][ch]<<" , Slope of the CI at Calib = "<<Gain_Cal_Lo[mod][view][pln][ch]<<endl<<endl;
	  //}
	  //}
#endif
	  //comment ajoute t'on le bon fit pour logain channel? On a une charge injectée pour high => on va d'abord collecter CI corrigée, pour le high gain channel. Ensuite, on va voir combien d'ADC low channel y correspondent. Cela nous donnera le logain. Puis on regarde combien as t'on d'ADC dans la logain channel et fait la correction. Attention si quand on passe dun graph à l'autre, on va dans une zone ultra non linéaire! Porrait etre mal fittée
	  /* double Eq_CI_ADC_Lo= Lo_Fit[mod][view][pln][ch]->antecedent(CI_ADC);
	     logain[mod][view][pln][ch]=gain_correction[mod][view][pln][ch]*(Eq_CI_ADC_Lo/inghitsum->adc);//ça c'est bon par contre car on compare low et high adc channels

	  //to be looked on
	  //prendre la charge équivalenter aux adc mesurés
	  double CI_ADC_Lo=Lo_Fit[mod][view][pln][ch]->Eval(inghitsum->loadc);
	  cout<<CI_ADC_Lo<<endl;
	  cout<<"Lo Gain was ="<<logain[mod][view][pln][ch]<<endl;
	  //il faudrait corriger le gain à la charge mesuré en regardant chargeinj/adc mesuré aux 2 points
	  logain_correction[mod][view][pln][ch]=(CI_ADC_Lo/inghitsum->loadc)/(CI_ADC/Eq_CI_ADC_Lo);
	  logain_correction[mod][view][pln][ch]*=logain[mod][view][pln][ch];
	  cout<<"Low Gain corrected is now = "<<logain[mod][view][pln][ch]<<endl;
	  */
	  //ATTENTION: pas fait: channel par channel, je n'ai pas refait le logain qui doit pas etre évalué a 1/10e du gain au début et? est ce que c'est bien le rapport y/x pour avoir le gain?
	  //##### Before conversion from ADC to #p.e.,#######
	  //##### we have to do correction of linearity #####
	  //##### of ADC channel.             ###############
	  //##### Now(Jan.2010), we ignore it ###############
	  //#################################################
	    OldGain->Fill(gain[mod][view][pln][ch]);
	    NewGain->Fill(gain_corrected[mod][view][pln][ch]);
	    double ADC_V=inghitsum->adc - ped[mod][view][pln][ch];
	    double LoADC_V=inghitsum->loadc - loped[mod][view][pln][ch];
	    double weight=gain_corrected[mod][view][pln][ch]/logain_corrected[mod][view][pln][ch];/*/gain[mod][view][pln][ch];*/
	    double Lope= (inghitsum->loadc - loped[mod][view][pln][ch])/logain_corrected[mod][view][pln][ch];
	    double Hipe= (inghitsum->adc - ped[mod][view][pln][ch])/gain_corrected[mod][view][pln][ch];

	    if((gain_corrected[mod][view][pln][ch]!=0) && (logain_corrected[mod][view][pln][ch]!=0) && (Hipe>.1) && (Lope>.1)) {
	      HiLine->Fill(gain_corrected[mod][view][pln][ch], gain_line[mod][view][pln][ch][LoRange]);
	      GNL->Fill(ADC_V,weight);
	      Evt0->Fill((inghitsum->adc - ped[mod][view][pln][ch])/gain[mod][view][pln][ch]);
	      Evt->Fill(Hipe);
	      Evt20->Fill((inghitsum->loadc - loped[mod][view][pln][ch])/gain[mod][view][pln][ch]);
	      Evt2->Fill(Lope);
	      LoHi0->Fill((inghitsum->adc - ped[mod][view][pln][ch])/gain[mod][view][pln][ch],ADC_V);
	      LoHi->Fill(Hipe,ADC_V);
	      Lo0->Fill((inghitsum->loadc - loped[mod][view][pln][ch])/logain[mod][view][pln][ch],LoADC_V);
              Lo->Fill(Lope,LoADC_V);
	      //cout<<Lope<<"  "<<(inghitsum->loadc - loped[mod][view][pln][ch])/logain_corrected[mod][view][pln][ch]<<endl;
	      LoHipe->Fill(Lope,Hipe);
	      LoHipeStd->Fill((inghitsum->loadc - loped[mod][view][pln][ch])/logain[mod][view][pln][ch],(inghitsum->adc - ped[mod][view][pln][ch])/gain[mod][view][pln][ch]);
	    }

	  //##### Conversion from ADC to #p.e. ##############
	  //#################################################
	  //cout<<"Conversion from ADC to #p.e"<<endl;
	  inghitsum -> pe   = 1.0 * ( inghitsum ->  adc - ped[mod][view][pln][ch] ) / gain_corrected[mod][view][pln][ch];

	  //##### Conversion from ADC to #p.e.(Logain)#######
	  //#################################################
	  //cout<<"Conversion from ADC to #p.e.(Logain)"<<endl;
	  inghitsum -> lope = 1.0 * ( inghitsum ->  loadc - loped[mod][view][pln][ch] ) / logain_corrected[mod][view][pln][ch] ;
	  //##### Conversiont from TDC to time[nsec] ########
	  //#################################################
	  //cout<<"Conversion from TDC to time[nsec]"<<endl;
	  long time = 2.5 * ( inghitsum ->  tdc ) - 10.0 * ( summary -> trgtime ); 
	  //##### If we don't have Pulse Per Second, ########
	  //##### time is larger than one second     ########
	  //#################################################
	  while(time>1000000000){
	    time -= 1000000000;
	  }
	  while(cosmic&&time>1000000000-100000000){
	    time -= 1000000000;
	  }

	  //##### We have to do  correction because of ##########
	  //##### difference of cable length b/w  ###############
	  //##### Back end board and front end board  ###########
	  //##### but some VETO channels should be careful to do 
	  //#####################################################
	  //cout<<"Starting cable length corrections"<<endl;
	  float cTimeCorrection;
	  if(!cosmic)
	    cTimeCorrection = cTimeCorrectionBase;
	  if(cosmic)
	    cTimeCorrection = 0.5*cTimeCorrectionBase;
	  switch ( mod ) {
	  case  1:
	    if( pln != 11 ) //#### Because pln 11 at mod 1 is pln 12 at mod 0
	    time -= cTimeCorrection;
	    break;
	  case  2:
	    time -= cTimeCorrection;
	    break;
	  case  3:
	    time -= cTimeCorrection;
	    break;	   
	  case  4:
	    if( pln == 11 ) //#### Because pln 11 at mod 4 is pln 12 at mod 3
	    time -= cTimeCorrection;
	    break;	   
	  case 11:
	    if( pln != 13)  //#### Because pln 13 at mod 11 is pln 14 at mod 10
	    time -= cTimeCorrection;
	    break;
	  case 12:
	    if( pln == 13)  //#### Because pln 13 at mod 12 is pln 14 at mod 11
	    time -= cTimeCorrection;
	    break;
	  defalut:
	    break;
	  }
	  inghitsum -> time = time;

	}//Hit Loop
      }//Cyc
    }//Mod
    wsummary = summary;
    wtree -> Fill();
    if(ievt%1000==0)wtree->AutoSave();
  }

  OldGain->SaveAs("Calib/oldgain.root");
  NewGain->SaveAs("Calib/newgain.root");
  GNL->Sumw2();
  GNL->Divide(Evt);
  GNL->SaveAs("Calib/gnl.root");
  LoHi->Sumw2();
  LoHi->Divide(Evt);
  LoHi->SaveAs("Calib/lohi.root");
  LoHi0->Sumw2();
  LoHi0->Divide(Evt0);
  LoHi0->SaveAs("Calib/lohi0.root");
  Lo->Sumw2();
  Lo->Divide(Evt2);
  Lo->SaveAs("Calib/lo.root");
  Lo0->Sumw2();
  Lo0->Divide(Evt20);
  Lo0->SaveAs("Calib/lo0.root");
  LoHipe->SaveAs("Calib/lohipe.root");
  LoHipeStd->SaveAs("Calib/lohipestd.root");
  HiLine->SaveAs("Calib/hiline.root");
  wtree -> Write();
  wfile -> Write();
  wfile -> Close();
  //app.Run();
}
 
