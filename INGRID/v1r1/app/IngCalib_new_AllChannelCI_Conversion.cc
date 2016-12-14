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
char HGGainName[100],HGNormName[100],HGGainLineName[100],HGNormLineName[100],LGGainLineName[100],LGNormLineName[100],LGCalName[100],LGCalNormName[100];
double NStep=1000;
double CIValue;

const int online_calib_run = 13047;
float     High_ADC_Cal[ nMod ][ nView ][ nPln ][ nCh ];
float     Low_ADC_Cal[ nMod ][ nView ][ nPln ][ nCh ];
float     gain[ nMod ][ nView ][ nPln ][ nCh ];
float     ped[ nMod ][ nView ][ nPln ][ nCh ];
float     logain[ nMod ][ nView ][ nPln ][ nCh ];
float     loped[ nMod ][ nView ][ nPln ][ nCh ];
float     gain_correction[ nMod ][ nView ][ nPln ][ nCh ];
float     Gain_Cal[ nMod ][ nView ][ nPln ][ nCh ];
float     logain_cal[ nMod ][ nView ][ nPln ][ nCh ];
float     SlopeCal[ nMod ][ nView ][ nPln ][ nCh ];
int     NCal[ nMod ][ nView ][ nPln ][ nCh ];
float     logain_cal2[ nMod ][ nView ][ nPln ][ nCh ];
float     SlopeCal2[ nMod ][ nView ][ nPln ][ nCh ];
int     NCal2[ nMod ][ nView ][ nPln ][ nCh ];
float     logain_correction[ nMod ][ nView ][ nPln ][ nCh ];
//float     Gain_Cal_Lo[ nMod ][ nView ][ nPln ][ nCh ];
float     gain_corrected[ nMod ][ nView ][ nPln ][ nCh ];
float     logain_corrected[ nMod ][ nView ][ nPln ][ nCh ];
float     LoGain[ nMod ][ nView ][ nPln ][ nCh ];
float     LoGain2[ nMod ][ nView ][ nPln ][ nCh ];
float value[ nMod ][ nView ][ nPln ][ nCh ];
float gain_line[ nMod ][ nView ][ nPln ][ nCh ][100];
int NPoints[ nMod ][ nView ][ nPln ][ nCh ];
int NPoints2[ nMod ][ nView ][ nPln ][ nCh ];
float NHiADC[ nMod ][ nView ][ nPln ][ nCh ];
float     gain_cal[ nMod ][ nView ][ nPln ][ nCh ];
float NLoADC[ nMod ][ nView ][ nPln ][ nCh ];
float NPEHi[ nMod ][ nView ][ nPln ][ nCh ][100];
float logain_line[ nMod ][ nView ][ nPln ][ nCh ][100];
float logain_calib[ nMod ][ nView ][ nPln ][ nCh ];
float CI2PE[ nMod ][ nView ][ nPln ][ nCh ];

TH1D * TotHG= new TH1D("TotHG","norm",1000,0,1000);
TH1D * TotHGNorm= new TH1D("HGNorm","norm",1000,0,1000);
TH1D * TotHGLine= new TH1D("TotHGLine","norm",100,0,100);
TH1D * TotHGLineNorm= new TH1D("HGLineNorm","norm",100,0,100);
TH1D * TotLGLine= new TH1D("TotLGLine","norm",100,0,100);
TH1D * TotLGLineNorm= new TH1D("LGLineNorm","norm",100,0,100);
TH1D * TotLGCal=new TH1D("LGCal","norm",100,0,100);
TH1D * TotLGCalNorm=new TH1D("LGCalNorm","norm",100,0,100);

TGraph *Hi_Graph;
TGraph *Lo_Graph;

TF1 * Hi_Fit[ nMod ][ nView ][ nPln ][ nCh ];
TF1 * Lo_Fit[ nMod ][ nView ][ nPln ][ nCh ];

TH1D *  HGGain[ nMod ][ nView ][ nPln ][ nCh ];
TH1D *  HGNorm[ nMod ][ nView ][ nPln ][ nCh ];
TH1D *  HGGainLine[ nMod ][ nView ][ nPln ][ nCh ];
TH1D *  HGNormLine[ nMod ][ nView ][ nPln ][ nCh ];
TH1D *  LGGainLine[ nMod ][ nView ][ nPln ][ nCh ];
TH1D *  LGNormLine[ nMod ][ nView ][ nPln ][ nCh ];
TH1D * LGCal[ nMod ][ nView ][ nPln ][ nCh ];
TH1D * LGCalNorm[ nMod ][ nView ][ nPln ][ nCh ];

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
	  logain_cal   [mod][view][pln][ch] = 0;
	  SlopeCal [mod][view][pln][ch] = 0;
	  NCal [mod][view][pln][ch] = 0;
          logain_cal2   [mod][view][pln][ch] = 0;
          SlopeCal2 [mod][view][pln][ch] = 0;
          NCal2 [mod][view][pln][ch] = 0;
	  LoGain   [mod][view][pln][ch] = 0;
	  LoGain2   [mod][view][pln][ch] = 0;
	  value [mod][view][pln][ch] = 0;
	  logain_calib  [mod][view][pln][ch] =  1e5;
	  NPoints[mod][view][pln][ch]=0;
	  NPoints2[mod][view][pln][ch]=0;
	  NHiADC[mod][view][pln][ch]=0;
	  NLoADC[mod][view][pln][ch]=0;
	  gain_cal[mod][view][pln][ch]=0;
	  for(int i=0;i<100;i++){      
	    gain_line[mod][view][pln][ch][i] =0; 
	    //NPoints[mod][view][pln][ch][i]=0;
	    //NHiADC[mod][view][pln][ch][i]=0;
	    //NLoADC[mod][view][pln][ch][i]=0;
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


void GetHisto(){
  
 for(int imod=0; imod<nMod; imod++){
    for(int iview=0; iview<nView; iview++){
      for(int ipln=0; ipln<nPln; ipln++){
	for(int ich=0; ich<nCh ; ich++){
	  if( GoodCh[imod][iview][ipln][ich] ){
#ifdef DEBUG
	    //	    cout<<"new iteration = Module:"<<imod<<"  , Plane:"<<ipln<<"  , Channel:"<<ich<<endl;
#endif	    
	    sprintf(hgraph,"hiGraph_%d",NDChan[imod][iview][ipln][ich]);
	    sprintf(lgraph,"loGraph_%d",NDChan[imod][iview][ipln][ich]);
	    Hi_Graph=(TGraph*) CI_File->Get(hgraph);
	    Lo_Graph=(TGraph*) CI_File->Get(lgraph);
	    
	    Hi_Fit[imod][iview][ipln][ich]=(TF1*)Hi_Graph->GetFunction("fFitHi");
	    Lo_Fit[imod][iview][ipln][ich]=(TF1*)Lo_Graph->GetFunction("fFitLo");
	    double Hi_CI_Cal=Hi_Fit[imod][iview][ipln][ich]->Eval(gain[imod][iview][ipln][ich]);
	    double Lo_CI_Cal=Lo_Fit[imod][iview][ipln][ich]->Eval(logain[imod][iview][ipln][ich]);
#ifdef DEBUG
	    //cout<<"At INGRID CAlib Point, ADC in High Gain Ch is = "<< High_ADC_Cal[imod][iview][ipln][ich] - ped[imod][iview][ipln][ich]<<" , and in Low Gain Ch, ADC = "<<Low_ADC_Cal[imod][iview][ipln][ich] - loped[imod][iview][ipln][ich]<<endl;
	    //cout<<"At INGRID CAlib Point, equivalent CI in High Gain Channel is = "<<Hi_CI_Cal<<" , and Low Gain CI is = "<<Lo_CI_Cal<<endl;
#endif	
	    Gain_Cal[imod][iview][ipln][ich]=Hi_CI_Cal/gain[imod][iview][ipln][ich];/*/(High_ADC_Cal[imod][iview][ipln][ich]-ped[imod][iview][ipln][ich])*/
	    if(Gain_Cal[imod][iview][ipln][ich]!=0) CI2PE[imod][iview][ipln][ich]=1/Hi_CI_Cal;
	    else CI2PE[imod][iview][ipln][ich]=.1;
	    //cout<<"CI2PE factor= "<<CI2PE[imod][iview][ipln][ich]<<endl;
	    //Gain_Cal_Lo[imod][iview][ipln][ich]=Lo_CI_Cal/(logain[imod][iview][ipln][ich]);
#ifdef DEBUG
	    //cout<<"Low Gain Original Value = "<<logain[imod][iview][ipln][ich]<<" ,  Low Gain Value = "<<LoGain[imod][iview][ipln][ich]<<endl;
	    //cout<<Gain_Cal[imod][iview][ipln][ich]<<"   "<<Gain_Cal_Lo[imod][iview][ipln][ich]<<endl<<endl;
#endif
	  }
	}
      }
    }
 }
}

//faire de gain corrected un TGraph que l'on va ensuite fitter au mieux? ou directement une fonction si possible... Comment on fait dans ce cas?
void GetHGLinearization(){
 for(int imod=0; imod<nMod; imod++){
    for(int iview=0; iview<nView; iview++){
      for(int ipln=0; ipln<nPln; ipln++){
	for(int ich=0; ich<nCh ; ich++){
	  if( GoodCh[imod][iview][ipln][ich] ){
	    //if(TMath::Abs(gain[iview][ipln][ich]-10)<7 && (Gain_Cal[mod][view][pln][ch] <1.7)){
	    sprintf(HGGainName,"HGGain%02d%d%02d%02d",imod,iview,ipln,ich);
	    sprintf(HGNormName,"HGNorm%02d%d%02d%02d",imod,iview,ipln,ich);
	    sprintf(HGGainLineName,"HGGainLine%02d%d%02d%02d",imod,iview,ipln,ich);
	    sprintf(HGNormLineName,"HGNormLine%02d%d%02d%02d",imod,iview,ipln,ich);
	    HGGain[imod][iview][ipln][ich] = new TH1D(HGGainName,"Gain",1000,0,1000);
	    HGNorm[imod][iview][ipln][ich] = new TH1D(HGNormName,"Norm",1000,0,1000);
	    HGGainLine[imod][iview][ipln][ich] = new TH1D(HGGainLineName,"Gain",1000,0,1000);
	    HGNormLine[imod][iview][ipln][ich] = new TH1D(HGNormLineName,"Norm",1000,0,1000);
	    sprintf(LGCalName,"LGCal%02d%d%02d%02d",imod,iview,ipln,ich);
	    sprintf(LGCalNormName,"LGCalNorm%02d%d%02d%02d",imod,iview,ipln,ich);
	    LGCal[imod][iview][ipln][ich] = new TH1D(LGCalName,"Gain",1000,0,1000);
	    LGCalNorm[imod][iview][ipln][ich] = new TH1D(LGCalNormName,"Norm",1000,0,1000);
	    sprintf(LGGainLineName,"LGGainLine%02d%d%02d%02d",imod,iview,ipln,ich);
	    sprintf(LGNormLineName,"LGNormLine%02d%d%02d%02d",imod,iview,ipln,ich);
	    LGGainLine[imod][iview][ipln][ich] = new TH1D(LGGainLineName,"Gain",1000,0,1000);
	    LGNormLine[imod][iview][ipln][ich] = new TH1D(LGNormLineName,"Norm",1000,0,1000);
	    
	    NStep=HGGain[imod][iview][ipln][ich]->GetNbinsX();

	    if(TMath::Abs(gain[imod][iview][ipln][ich]-10.)<7. && (Gain_Cal[imod][iview][ipln][ich] <1.7)){
	    for(int iadc=1;iadc<=NStep;iadc++){
	      double ADCNumber =iadc;
	      CIValue=Hi_Fit[imod][iview][ipln][ich]->Eval(ADCNumber);
	      if(CIValue==0) gain_correction[imod][iview][ipln][ich]=Gain_Cal[imod][iview][ipln][ich]/(Hi_Fit[imod][iview][ipln][ich]->Eval(699)/699);
	      else    gain_correction[imod][iview][ipln][ich]=Gain_Cal[imod][iview][ipln][ich]/(CIValue/(ADCNumber));//R1/R2            
	      gain_corrected[imod][iview][ipln][ich]=gain[imod][iview][ipln][ich]* gain_correction[imod][iview][ipln][ich];
	      HGGain[imod][iview][ipln][ich]->Fill(ADCNumber,gain_corrected[imod][iview][ipln][ich]);
	      HGNorm[imod][iview][ipln][ich]->Fill(ADCNumber);
	      if(gain_corrected[imod][iview][ipln][ich]>1000){
		cout<<"ADC = "<<ADCNumber<<" , Original Gain = "<<gain[imod][iview][ipln][ich]<<",Gain corrected = "<<gain_corrected[imod][iview][ipln][ich]<<" , At iadc point, CI = "<<CIValue<<" ,at Cal Point, slope is = "<<Gain_Cal[imod][iview][ipln][ich]<<" ,NDChan = "<<NDChan[imod][iview][ipln][ich]<<endl;   
		//cout<<Hi_Fit[imod][iview][ipln][ich]->Eval(gain[imod][iview][ipln][ich])<<"  "<<Hi_Fit[imod][iview][ipln][ich]->Eval(ADCNumber)<<endl;
	      }//creation of linearized gain (over 1.pe around)
	      int line=(int) ADCNumber/10;
	      HGGainLine[imod][iview][ipln][ich]->Fill(line,gain_corrected[imod][iview][ipln][ich]);
	      HGNormLine[imod][iview][ipln][ich]->Fill(line,1);
	      TotHG->Fill(ADCNumber,gain_corrected[imod][iview][ipln][ich]);
	      TotHGNorm->Fill(ADCNumber);
              TotHGLine->Fill(line,gain_corrected[imod][iview][ipln][ich]);
              TotHGLineNorm->Fill(line);
	    }
	      //    cout<<"Amount of ADC counts = "<<iadc<<" , Original Gain = "<<gain[imod][iview][ipln][ich]<<", Gain corrected is = "<<gain_corrected[imod][iview][ipln][ich]<<endl;
	    }
	  }
	}
      }
    }
 }
 TotHG->Sumw2();
 TotHG->Divide(TotHGNorm);
 TotHG->SaveAs("Calib/HG.root");
 TotHGLine->Sumw2();
 TotHGLine->Divide(TotHGLineNorm);
 TotHGLine->SaveAs("Calib/HGLine.root");

 for(int imod=0; imod<nMod; imod++){
   for(int iview=0; iview<nView; iview++){
     for(int ipln=0; ipln<nPln; ipln++){
       for(int ich=0; ich<nCh ; ich++){
	 if( GoodCh[imod][iview][ipln][ich] ){
	   if(TMath::Abs(gain[imod][iview][ipln][ich]-10.)<7. && (Gain_Cal[imod][iview][ipln][ich] <1.7)){
	   HGGain[imod][iview][ipln][ich]->Sumw2();
	   HGGain[imod][iview][ipln][ich]->Divide(HGNorm[imod][iview][ipln][ich]);
	   HGGainLine[imod][iview][ipln][ich]->Sumw2();
	   HGGainLine[imod][iview][ipln][ich]->Divide(HGNormLine[imod][iview][ipln][ich]);
	   }
	   
	 }
       }
     }
   }
 }
}

void GetLGLinearization(){
 for(int imod=0; imod<nMod; imod++){
    for(int iview=0; iview<nView; iview++){
      for(int ipln=0; ipln<nPln; ipln++){
	for(int ich=0; ich<nCh ; ich++){
	  if( GoodCh[imod][iview][ipln][ich] ){
	    if(TMath::Abs(gain[imod][iview][ipln][ich]-10.)<5. && (Gain_Cal[imod][iview][ipln][ich] <1.5)){
	    
	      int NBins= LGGainLine[imod][iview][ipln][ich]->GetNbinsX();
	    //prendre le(s) bin(s) de calibration: regarder l'histo d'abord. Prendre peut etre les 10 premiers (en gros de 0 a 10 p.e )

	      double LoGainCal=0;
	      double Gain_Cal_Lo=0;
	      double CI_ADC_Lo=0;
	      
	      for(int ilin=0;ilin<NBins;ilin++){
		
		CI_ADC_Lo=Lo_Fit[imod][iview][ipln][ich]->Eval(ilin);
		if(CI_ADC_Lo==0) logain_corrected[imod][iview][ipln][ich]=logain[imod][iview][ipln][ich];//case over CI calibration
		else logain_corrected[imod][iview][ipln][ich]=(ilin/(CI_ADC_Lo*CI2PE[imod][iview][ipln][ich]));
		//cout<<"Logain_corrected = "<<logain_corrected[imod][iview][ipln][ich]<<endl;
		LGGainLine[imod][iview][ipln][ich]->Fill(ilin,logain_corrected[imod][iview][ipln][ich]);
		LGNormLine[imod][iview][ipln][ich]->Fill(ilin,1);
		if(logain_corrected[imod][iview][ipln][ich]!=logain_corrected[imod][iview][ipln][ich])
		  { 
		    cout<<"Module = "<<imod<<" , Plane = "<<ipln<<" , Ch = "<<ich<<" , View = "<<iview<<endl;
		    cout<<"PE="<<ilin<<" , LoGainCal="<<LoGainCal<<" , Gain_Cal_Lo="<<Gain_Cal_Lo<<" , Slope at the point we observed = "<<(CI_ADC_Lo/(ilin))<<" , CIADC="<<CI_ADC_Lo<<" , NCal="<<NCal[imod][iview][ipln][ich]<<endl<<endl;
		  }
		/*		if(ilin<10 && ilin>2){
		  cout<<"PE="<<ilin<<" , LoGainCal="<<LoGainCal<<" , Gain_Cal_Lo="<<Gain_Cal_Lo<<" , Slope at the point we observed = "<<(CI_ADC_Lo/(ilin))<<" , CIADC="<<CI_ADC_Lo<<" , NCal="<<NCal[imod][iview][ipln][ich]<<endl<<endl;
		  }*/
		TotLGLine->Fill(ilin,logain_corrected[imod][iview][ipln][ich]);
		TotLGLineNorm->Fill(ilin);
		  

	      }
	    }
	  }
	}
      }
    }
 }
 
 TotLGLine->Sumw2();
 TotLGLine->Divide(TotLGLineNorm);
 TotLGLine->SaveAs("Calib/LGLine.root");
 
 for(int imod=0; imod<nMod; imod++){
   for(int iview=0; iview<nView; iview++){
     for(int ipln=0; ipln<nPln; ipln++){
       for(int ich=0; ich<nCh ; ich++){
	 if( GoodCh[imod][iview][ipln][ich] ){
	   if(TMath::Abs(gain[imod][iview][ipln][ich]-10.)<7. && (Gain_Cal[imod][iview][ipln][ich] <1.7)){
	     LGGainLine[imod][iview][ipln][ich]->Sumw2();
	     LGGainLine[imod][iview][ipln][ich]->Divide(LGNormLine[imod][iview][ipln][ich]);
	   }
	   
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
  TH1D * LoHi     = new TH1D("LoHi","Lo/Hi p.e in function of hipe",100,0,100);
  TH1D * HiLo     = new TH1D("HiLo","Lo/Hi p.e in fucntion of lope",100,0,100);
  TH1D * NormMean     = new TH1D("NormMean","Lo/Hi p.e in fucntion of Meanp.e",100,0,100);
  TH1D * NormMeanOld     = new TH1D("NormMeanOld","Lo/Hi p.e in function of Meanp.e",100,0,100);
  TH1D * HiLoMean     = new TH1D("HiLoMean","Lo/Hi p.e in fucntion of Meanp.e",100,0,100);
  TH1D * HiLoMeanOld     = new TH1D("HiLoMeanOld","Lo/Hi p.e in function of Meanp.e",100,0,100);
  TH1D * Lo     = new TH1D("Lo","Lo p.e ",100,0,100);
  TH1D * Lo0     = new TH1D("Lo0","Lo p.e ",100,0,100);
  TH1D * LoCal     = new TH1D("LoCal","Lo p.e ",100,0,100);
  TH1D * LoEvt     = new TH1D("LoEvt","Lo p.e ",100,0,100);
  TH2D * HiLine     = new TH2D("HiLine","Lo/Hi p.e ",100,0,30,100,0,30);
  TH2D * LoHipe     = new TH2D("LoHipe","Lo/Hi p.e ",100,0,100,100,0,100);
  TH2D * LoHipeStd     = new TH2D("LoHipeStd","Lo/Hi p.e ",100,0,100,100,0,100);
  TH1D * hGainH = new TH1D("hGainH","High Gain for All Channels",100,8,13);
  TH1D * hSlopeH = new TH1D("hSlopeH","High Gain Channel Slope",100,5,20);
  TH1D * hGainL = new TH1D("hGainL","Low Gain for All Channels",50,0,3);
  TH1D * hSlopeL = new TH1D("hSlopeL","Low Gain Channel Slope",100,0,5);
  TH1D * bGainH = new TH1D("bGainH","High Gain for All Channels",100,8,13);
  TH1D * bSlopeH = new TH1D("bSlopeH","High Gain Channel Slope",100,5,20);
  TH1D * bGainL = new TH1D("bGainL","Low Gain for All Channels",50,0,3);
  TH1D * bSlopeL = new TH1D("bSlopeL","Low Gain Channel Slope",100,0,5);
  TH1D * Hpe10 = new TH1D("Hpe10","p.e in HG Chan",100,0,100);
  TH1D * Lpe10 = new TH1D("Lpe10","p.e in LG Chan",500,0,500);
  TH1D * OldHpe10 = new TH1D("Hpe10Old","p.e in HG Chan, before tuning",100,0,100);
  TH1D * OldLpe10 = new TH1D("Lpe10Old","p.e in LG Chan, before tuning",500,0,500);
  TH1D * PENew = new TH1D("PENew","#p.e",500,0,500);
  TH1D * PEOld = new TH1D("PEOld","#p.e",500,0,500);

  //double CI_ADC=0;
  TH1D * nd     = new TH1D("nd","Lo p.e ",24000,0,24000);
  //float NPoints[100];

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
  cout<<"Initializing Constants"<<endl;
  InitilizeCalibConst(); //### Initialize gain, pedestal, and so on


  //#### Read file before calibration ######
  //########################################
  IngridEventSummary* summary = new IngridEventSummary();

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

  TFile*            wfile = new TFile(Output,"recreate");
  TTree*            wtree = new TTree("tree","tree");
  IngridEventSummary* wsummary = new IngridEventSummary();
  wtree              -> Branch   ("fDefaultReco.","IngridEventSummary",
				  &wsummary,  64000,   99);


  cout<<"Get Calibration Constants"<<endl;
  if(!semioffline)
    GetCalibConst(run_number,sub_run_number);
  else if(semioffline)
    GetCalibConst(online_calib_run,0);


#ifdef DEBUG
  //cout<<" Calib non linearized is done"<<endl;
#endif
  CI_File=new TFile("/home/bquilain/Ingrid_Process/CISimpleOutFile.root");

  Hi_Graph=NULL;
  Lo_Graph=NULL;
  
  cout<<"Get CI Histograms"<<endl;
  GetHisto();
  cout<<"Linearize the gain in High ADC channel"<<endl;
  GetHGLinearization();

#ifdef DEBUG
  //cout<<"Calib linearized is done"<<endl;
#endif

  GetLGLinearization();



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
	cout<<"Module:"<<mod<<"   , Cycle:"<<cyc<<endl;
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
	  int LoRange=(int) (inghitsum->adc - ped[mod][view][pln][ch])/10;

	  //logain_corrected[mod][view][pln][ch]=LoGain[mod][view][pln][ch]*(Gain_Cal_Lo[mod][view][pln][ch])/(Lo_Fit[mod][view][pln][ch]->Eval(inghitsum->loadc - loped[mod][view][pln][ch])/(inghitsum->loadc - loped[mod][view][pln][ch]));
#ifdef DEBUG
	  //cout<<"Default Low Gain was ="<<logain[mod][view][pln][ch]<<" , Low Gain is = "<<LoGain[mod][view][pln][ch]<<endl;
	  //if(inghitsum->adc>400) cout<<"High Gain was ="<<gain[mod][view][pln][ch]<<"  , Charge Injected is = "<<CI_ADC<<"  ,  ADC equivalent is = "<<inghitsum->adc<<endl;
	  //	  cout<<"High Gain was ="<<gain[mod][view][pln][ch]<<"  , Charge Injected is = "<<CI_ADC<<"  ,  ADC equivalent is = "<<inghitsum->adc<<endl;
#endif
	  
	  //double CI_ADC_Lo=Lo_Fit[mod][view][pln][ch]->Eval(inghitsum->loadc - loped[mod][view][pln][ch]);
          //logain[mod][view][pln][ch]=(1/Conversion_CIPE)*(CI_ADC_Lo/(inghitsum->loadc - loped[mod][view][pln][ch]));
	  //double Lope= (inghitsum->loadc - loped[mod][view][pln][ch])/logain_corrected[mod][view][pln][ch];
	  //double Hipe= (inghitsum->adc - ped[mod][view][pln][ch])/gain_corrected[mod][view][pln][ch];
																			  
 //logain_correction[mod][view][pln][ch]=(Lo_Fit[mod][view][pln][ch]->Eval(inghitsum->loadc - loped[mod][view][pln][ch])/(inghitsum->loadc - loped[mod][view][pln][ch]))/Gain_Cal_Lo[mod][view][pln][ch];
																			   //logain_corrected[mod][view][pln][ch]=logain[mod][view][pln][ch]*logain_correction[mod][view][pln][ch];
#ifdef DEBUG
	  //if(CI_ADC>200 || CI_ADC_Lo>200) {
	  //cout<<"High Gain was ="<<gain[mod][view][pln][ch]<<"  , Charge Injected is = "<<CI_ADC<<"  ,  ADC equivalent is = "<<inghitsum->adc - ped[mod][view][pln][ch]<<endl;
	  //cout<<"High Gain corrected is now = "<<gain_corrected[mod][view][pln][ch]<<endl;
	  //cout<<"Conversion Factor = "<<Conversion_CIPE<<endl;
	  //cout<<"Lo Gain was ="<<logain[mod][view][pln][ch]<<"  , Charge Injected is = "<<CI_ADC_Lo<<"  ,  ADC equivalent is = "<<inghitsum->loadc - loped[mod][view][pln][ch]<<endl;
	  //cout<<"LoGain corrected is now = "<<logain_corrected[mod][view][pln][ch]<<endl<<endl;
	  /* if((TMath::Abs(Hipe/Lope)>1.5 || TMath::Abs(Hipe/Lope)<.75) && (Lope>25||Hipe>25) && Lope<1000){
	  cout<<"Number of ADC in HG Channel = "<<inghitsum->adc - ped[mod][view][pln][ch]<<endl;
	  cout<<"Before Calibration, Amount of p.e is in HG channel = "<<(inghitsum->adc - ped[mod][view][pln][ch])/gain[mod][view][pln][ch]<<" ,  p.e in Lo Gain Channel = "<<(inghitsum->loadc - loped[mod][view][pln][ch])/logain[mod][view][pln][ch]<<endl;
	  cout<<"Amount of p.e is in HG channel = "<<(inghitsum->adc - ped[mod][view][pln][ch])/gain_corrected[mod][view][pln][ch]<<" ,  p.e in Lo Gain Channel = "<<(inghitsum->loadc - loped[mod][view][pln][ch])/logain_corrected[mod][view][pln][ch]<<endl;
	  //cout<<"Logain at Calib Point = "<<logain_calib[mod][view][pln][ch]<<" , Slope of the CI at Calib = "<<Gain_Cal_Lo[mod][view][pln][ch]<<endl<<endl;
	  cout<<"High Gain Channel"<<endl;
	  cout<<"Original gain = "<<gain[mod][view][pln][ch]<<" , gain corrected = "<<gain_corrected[mod][view][pln][ch]<<" , Slope of the CI at Calib = "<<Gain_Cal[mod][view][pln][ch]<<"  , CI = "<<CI_ADC<<endl;
	  cout<<"Low Gain Channel"<<endl;
	  if(logain[mod][view][pln][ch]==LoGain[mod][view][pln][ch]) cout<<"logain not well tuned"<<endl;
	  cout<<"Logain at Calib Point = "<<LoGain[mod][view][pln][ch]<<" , Slope of the CI at Calib = "<<Gain_Cal_Lo[mod][view][pln][ch]<<"  , CI = "<<CI_ADC_Lo<<endl<<endl;
	  }*/

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
	  /*OldGain->Fill(gain[mod][view][pln][ch]);
	    NewGain->Fill(gain_corrected[mod][view][pln][ch]);
	    double ADC_V=inghitsum->adc - ped[mod][view][pln][ch];
	    double LoADC_V=inghitsum->loadc - loped[mod][view][pln][ch];
	    double weight=gain_corrected[mod][view][pln][ch]/logain_corrected[mod][view][pln][ch];
	    //double Lope= (inghitsum->loadc - loped[mod][view][pln][ch])/logain_corrected[mod][view][pln][ch];
	    //double Hipe= (inghitsum->adc - ped[mod][view][pln][ch])/gain_corrected[mod][view][pln][ch];

	    if(((mod<16 && pln<11)||(mod==16 && pln<18)) && (gain_corrected[mod][view][pln][ch]!=0) && (logain_corrected[mod][view][pln][ch]!=0) && (Hipe>.1) && (Lope>.1) && (gain[mod][view][pln][ch]<20) && Gain_Cal[mod][view][pln][ch]<1.7) {//éviter les fails = pas top a corriger ensuite en dur dans le code
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
	    }*/
	  double Lope=(inghitsum->loadc - loped[mod][view][pln][ch])/logain[mod][view][pln][ch];
	  double Hipe=(inghitsum->adc - ped[mod][view][pln][ch])/gain[mod][view][pln][ch];
	  if( GoodCh[mod][view][pln][ch] ){//forcement assure car on parcourt les évenements
	    if(TMath::Abs(gain[mod][view][pln][ch]-10.)<5. && (Gain_Cal[mod][view][pln][ch] <1.5)){
	      if((NDChan[mod][view][pln][ch]!=6643)&&(NDChan[mod][view][pln][ch]!=13300)&&(NDChan[mod][view][pln][ch]!=11122)&&(NDChan[mod][view][pln][ch]!=15615)&&(NDChan[mod][view][pln][ch]!=3103)&&(NDChan[mod][view][pln][ch]!=15943)&&(NDChan[mod][view][pln][ch]!=11510)&&(NDChan[mod][view][pln][ch]!=13929)&&(NDChan[mod][view][pln][ch]!=14126)&&(NDChan[mod][view][pln][ch]!=636)&&(NDChan[mod][view][pln][ch]!=19734)){
	      int HLine=(inghitsum->adc - ped[mod][view][pln][ch])/10;
	      double HGCorr=HGGainLine[mod][view][pln][ch]->GetBinContent(HLine);
	      //cout<<HGCorr<<endl;
              int LLine=inghitsum->loadc - loped[mod][view][pln][ch];
              double LGCorr=LGGainLine[mod][view][pln][ch]->GetBinContent(LLine);
	      //cout<<LGCorr<<endl<<endl;
	      Lope= (inghitsum->loadc - loped[mod][view][pln][ch])/LGCorr;                     
	      double HipeOld=(inghitsum->adc - ped[mod][view][pln][ch])/gain[mod][view][pln][ch];                   
	      double LopeOld=(inghitsum->loadc - loped[mod][view][pln][ch])/logain[mod][view][pln][ch];
	      Hipe= (inghitsum->adc - ped[mod][view][pln][ch])/HGCorr;
	      if(Lope<10000 && Hipe<1000) {
		hGainH->Fill(gain[mod][view][pln][ch]);
		hSlopeH->Fill(HGCorr);
		hGainL->Fill(logain_cal[mod][view][pln][ch]);
		hSlopeL->Fill(LGCorr);
	      if(TMath::Abs(Hipe/Lope-1)>1) {
                bGainH->Fill(gain[mod][view][pln][ch]);
                bSlopeH->Fill(HGCorr);
                bGainL->Fill(logain_cal[mod][view][pln][ch]);
                bSlopeL->Fill(LGCorr);

		//cout<<"Module="<<mod<<" ,View="<<view<<" ,Plane="<<pln<<" ,Ch="<<ch<<" ,NDChan="<<NDChan[mod][view][pln][ch]<<endl;
		//cout<<"Hipe="<<Hipe<<" ,HGainCorrected="<<HGCorr<<" ,Lope="<<Lope<<" ,LGainCorrected="<<LGCorr<<endl;
		//cout<<"Original High Gain = "<<gain[mod][view][pln][ch]<<" ,Slope = "<<Gain_Cal[mod][view][pln][ch]<<" , Local = "<<logain_cal[mod][view][pln][ch]<<" , Slope of the LowGain = "<<SlopeCal[mod][view][pln][ch]<<endl<<endl;
		nd->Fill(NDChan[mod][view][pln][ch]);
	      }
	      if((Hipe+Lope)/2<40) PENew->Fill(Hipe);
	      else PENew->Fill(Lope);

	      if((inghitsum->adc - ped[mod][view][pln][ch])/gain[mod][view][pln][ch]<50) PEOld->Fill((inghitsum->adc - ped[mod][view][pln][ch])/gain[mod][view][pln][ch]);
              else PEOld->Fill((inghitsum->loadc - loped[mod][view][pln][ch])/logain[mod][view][pln][ch]);

	      Hpe10->Fill(Hipe);
	      Lpe10->Fill(Lope);
              OldHpe10->Fill((inghitsum->adc - ped[mod][view][pln][ch])/gain[mod][view][pln][ch]);
              OldLpe10->Fill((inghitsum->loadc - loped[mod][view][pln][ch])/logain[mod][view][pln][ch]);
	      Evt->Fill(Hipe);
	      HiLo->Fill(Lope,Hipe/Lope);
	      LoHi->Fill(Hipe,Hipe/Lope);
	      NormMean->Fill((Hipe+Lope)/2);
	      HiLoMean->Fill((Hipe+Lope)/2,Hipe/Lope);
	      NormMeanOld->Fill((HipeOld+LopeOld)/2);
	      HiLoMeanOld->Fill((HipeOld+LopeOld)/2,HipeOld/LopeOld);
	      LoHipe->Fill(Lope,Hipe);
	      LoHipeStd->Fill((inghitsum->loadc - loped[mod][view][pln][ch])/logain[mod][view][pln][ch],(inghitsum->adc - ped[mod][view][pln][ch])/gain[mod][view][pln][ch]);
	      }
	      }
	      else{
		gain_corrected[mod][view][pln][ch]=gain[mod][view][pln][ch];
		logain_corrected[mod][view][pln][ch]=logain[mod][view][pln][ch];
	      }
	    }
	    else{
	      gain_corrected[mod][view][pln][ch]=gain[mod][view][pln][ch];
	      logain_corrected[mod][view][pln][ch]=logain[mod][view][pln][ch];
	    }
	  }	    
	
	    //##### Conversion from ADC to #p.e. ##############
	  //#################################################
	  //cout<<"Conversion from ADC to #p.e"<<endl;
	  inghitsum -> pe   = Hipe;
	  //cout<<inghitsum -> pe<<endl;
	  //##### Conversion from ADC to #p.e.(Logain)#######
	  //#################################################
	  //cout<<"Conversion from ADC to #p.e.(Logain)"<<endl;
	  inghitsum -> lope = Lope;
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
  /*
  OldGain->SaveAs("Calib_Conv/oldgain.root");
  NewGain->SaveAs("Calib_Conv/newgain.root");
  nd->SaveAs("nd280.root");
  GNL->Sumw2();
  GNL->Divide(Evt);
  GNL->SaveAs("Calib_Conv/gnl.root");
  LoCal->Sumw2();
  LoCal->Divide(LoEvt);
  LoCal->SaveAs("Calib_Conv/LoCal.root");
  LoHi->Sumw2();
  LoHi->Divide(Evt);
  LoHi->SaveAs("Calib_Conv/lohi.root");
  LoHi0->Sumw2();
  LoHi0->Divide(Evt0);
  LoHi0->SaveAs("Calib_Conv/lohi0.root");
  Lo->Sumw2();
  Lo->Divide(Evt2);
  Lo->SaveAs("Calib_Conv/lo.root");
  Lo0->Sumw2();
  Lo0->Divide(Evt20);
  Lo0->SaveAs("Calib_Conv/lo0.root");*/
  LoHipe->SaveAs("Calib_Conv/lohipe.root");
  LoHipeStd->SaveAs("Calib_Conv/lohipestd.root");
  HiLine->SaveAs("Calib_Conv/hiline.root");
  nd->SaveAs("Calib_Conv/nd280.root");
  hGainH->Scale(1./hGainH->Integral());
  hGainH->SaveAs("Calib_Conv/hGainH.root");
  hSlopeH->Scale(1./hSlopeH->Integral());
  hSlopeH->SaveAs("Calib_Conv/hSlopeH.root");
  hGainL->Scale(1./hGainL->Integral());
  hGainL->SaveAs("Calib_Conv/hGainL.root");
  hSlopeL->Scale(1./hSlopeL->Integral());
  hSlopeL->SaveAs("Calib_Conv/hSlopeL.root");
  bGainH->Scale(1./bGainH->Integral());
  bGainH->SaveAs("Calib_Conv/bGainH.root");
  bSlopeH->Scale(1./bSlopeH->Integral());
  bSlopeH->SaveAs("Calib_Conv/bSlopeH.root");
  bGainL->Scale(1./bGainL->Integral());
  bGainL->SaveAs("Calib_Conv/bGainL.root");
  bSlopeL->Scale(1./bSlopeL->Integral());
  bSlopeL->SaveAs("Calib_Conv/bSlopeL.root");
  LoHi->Sumw2();                        
  LoHi->Divide(Evt);                                                                                 
  LoHi->SaveAs("Calib_Conv/lohi.root"); 
  HiLo->Sumw2();
  HiLo->Divide(Evt);
  HiLo->SaveAs("Calib_Conv/hilo.root");
  HiLoMean->Sumw2();
  HiLoMean->Divide(NormMean);
  HiLoMean->SaveAs("Calib_Conv/hilomean.root");
  HiLoMeanOld->Sumw2();
  HiLoMeanOld->Divide(NormMeanOld);
  HiLoMeanOld->SaveAs("Calib_Conv/hilomeanold.root");
  Hpe10->SaveAs("Calib_Conv/highperepartition.root");
  Lpe10->SaveAs("Calib_Conv/Lowperepartition.root");
  OldHpe10->SaveAs("Calib_Conv/OldHighperepartition.root");
  OldLpe10->SaveAs("Calib_Conv/OldLowperepartition.root");
  PENew->SaveAs("Calib_Conv/PENew.root");
  PEOld->SaveAs("Calib_Conv/PEOld.root");
  wtree -> Write();
  wfile -> Write();
  wfile -> Close();
  //app.Run();
}
 
