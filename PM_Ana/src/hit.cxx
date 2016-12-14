#include<iostream>
#include<sstream>
#include<fstream>
using namespace std;
#include <iomanip.h>
#include <sys/stat.h>

#include <TROOT.h>
#include <TStyle.h>
#include <TApplication.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TEventList.h>
#include <TBranch.h>
#include <TH1.h>
#include <TGraph.h>
#include <TGaxis.h>
#include <TMarker.h>
#include <TText.h>
#include <TSpectrum.h>
#include <TF1.h>
#include <TMath.h>
#include <TSpectrum.h>
#include <TF1.h>
#include <TMath.h>
#include <TSpectrum.h>
#include "TApplication.h"

#include "hit.hxx"


hit::hit(){
}

Int_t hit::one_hit_every_one_layer(Double_t *pe,Int_t numcyc,Double_t cut,Int_t *hit_channel_x,Int_t *hit_channel_y){
  hit_flag_x=0;
  hit_flag_y=0;
  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t numtfb=0;numtfb<UseNumTFB;numtfb++){
      for(Int_t i=0;i<NumLayer;i++){
	hit_layer[numtfb][i]=0;
      }
      for(Int_t numch=0;numch<LayerNumCh;numch++){
	if(*(pe+nummod*NumTFB*NumCh*NumCyc+numtfb*NumCh*NumCyc+numch*NumCyc+numcyc)>cut){
	  hit_layer[numtfb][0]++;
	}
      }
      for(Int_t numch=LayerNumCh;numch<UseNumCh;numch++){
	if(*(pe+nummod*NumTFB*NumCh*NumCyc+numtfb*NumCh*NumCyc+numch*NumCyc+numcyc)>cut){
	  hit_layer[numtfb][1]++;
	}
      }
    
      if(hit_layer[numtfb][0]==1){
	hit_flag_x++;
      }
      if(hit_layer[numtfb][1]==1){
	hit_flag_y++;
      }
    }//numtfb
  }//nummod
  if(hit_flag_x>UseNumTFB-1||hit_flag_y>UseNumTFB-1){
    for(Int_t nummod=0;nummod<NumMod;nummod++){
      for(Int_t numtfb=0;numtfb<NumTFB;numtfb++){
	for(Int_t numch=0;numch<LayerNumCh;numch++){
	  if(*(pe+nummod*NumTFB*NumCh*NumCyc+numtfb*NumCh*NumCyc+numch*NumCyc+numcyc)>cut){
	    *(hit_channel_x+numtfb)=numch;
	  }
	}
	for(Int_t numch=LayerNumCh;numch<UseNumCh;numch++){
	  if(*(pe+nummod*NumTFB*NumCh*NumCyc+numtfb*NumCh*NumCyc+numch*NumCyc+numcyc)>cut){
	    *(hit_channel_y+numtfb)=numch;
	  }
	}
      }//numtfb
    }//nummod
    return 1;
  }
  else{
    return 0;
  }
}

Int_t hit::one_hit_every_one_layer_x(Double_t *pe,Int_t numcyc,Double_t cut,Int_t *hit_channel_x){
  hit_flag_x=0;
  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t numtfb=0;numtfb<UseNumTFB;numtfb++){
      for(Int_t i=0;i<NumLayer;i++){
	hit_layer[numtfb][i]=0;
      }
      for(Int_t numch=0;numch<LayerNumCh;numch++){
	if(*(pe+nummod*NumTFB*NumCh*NumCyc+numtfb*NumCh*NumCyc+numch*NumCyc+numcyc)>cut){
	  hit_layer[numtfb][0]++;
	}
      }
      for(Int_t numch=LayerNumCh;numch<UseNumCh;numch++){
	if(*(pe+nummod*NumTFB*NumCh*NumCyc+numtfb*NumCh*NumCyc+numch*NumCyc+numcyc)>cut){
	  hit_layer[numtfb][1]++;
	}
      }
    
      if(hit_layer[numtfb][0]==1){
	hit_flag_x++;
      }
    }//numtfb
  }//nummod
  if(hit_flag_x>UseNumTFB-1){
    for(Int_t nummod=0;nummod<NumMod;nummod++){
      for(Int_t numtfb=0;numtfb<NumTFB;numtfb++){
	for(Int_t numch=0;numch<LayerNumCh;numch++){
	  if(*(pe+nummod*NumTFB*NumCh*NumCyc+numtfb*NumCh*NumCyc+numch*NumCyc+numcyc)>cut){
	    *(hit_channel_x+numtfb)=numch;
	  }
	}
      }//numtfb
    }//nummod
    return 1;
  }
  else{
    return 0;
  }
}

Int_t hit::one_hit_every_one_layer_y(Double_t *pe,Int_t numcyc,Double_t cut,Int_t *hit_channel_y){
  hit_flag_y=0;
  for(Int_t nummod=0;nummod<NumMod;nummod++){
    for(Int_t numtfb=0;numtfb<UseNumTFB;numtfb++){
      for(Int_t i=0;i<NumLayer;i++){
	hit_layer[numtfb][i]=0;
      }
      for(Int_t numch=0;numch<LayerNumCh;numch++){
	if(*(pe+nummod*NumTFB*NumCh*NumCyc+numtfb*NumCh*NumCyc+numch*NumCyc+numcyc)>cut){
	  hit_layer[numtfb][0]++;
	}
      }
      for(Int_t numch=LayerNumCh;numch<UseNumCh;numch++){
	if(*(pe+nummod*NumTFB*NumCh*NumCyc+numtfb*NumCh*NumCyc+numch*NumCyc+numcyc)>cut){
	  hit_layer[numtfb][1]++;
	}
      }
    
      if(hit_layer[numtfb][1]==1){
	hit_flag_y++;
      }
    }//numtfb
  }//nummod
  if(hit_flag_y>UseNumTFB-1){
    for(Int_t nummod=0;nummod<NumMod;nummod++){
      for(Int_t numtfb=0;numtfb<NumTFB;numtfb++){
	for(Int_t numch=LayerNumCh;numch<UseNumCh;numch++){
	  if(*(pe+nummod*NumTFB*NumCh*NumCyc+numtfb*NumCh*NumCyc+numch*NumCyc+numcyc)>cut){
	    *(hit_channel_y+numtfb)=numch;
	  }
	}
      }//numtfb
    }//nummod
    return 1;
  }
  else{
    return 0;
  }
}



Int_t hit::three_hit_AND_side_layer(Double_t *pe,Int_t nmod,Int_t ntfb,Int_t nch,Int_t ncyc,Double_t cut){
  hit_flag=0;

  if(nch<24){
  Int_t nch_layer = nch;
  if(ntfb==0){
    Int_t cuttfb_1=1,cuttfb_2=2;
    if((nch_layer)==0){
      Int_t cutch_1=0,cutch_2=1;
      if(*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_1*NumCh*NumCyc+cutch_1*NumCyc+ncyc)>cut||*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_1*NumCh*NumCyc+cutch_2*NumCyc+ncyc)>cut){
	hit_flag++;
      }
      if(*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_2*NumCh*NumCyc+cutch_1*NumCyc+ncyc)>cut||*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_2*NumCh*NumCyc+cutch_2*NumCyc+ncyc)>cut){
	hit_flag++;
      }
    }//(nch_layer-1)==0
    else if((nch_layer)==23){
      Int_t cutch_1=23,cutch_2=22;
      if(*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_1*NumCh*NumCyc+cutch_1*NumCyc+ncyc)>cut||*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_1*NumCh*NumCyc+cutch_2*NumCyc+ncyc)>cut){
	hit_flag++;
      }
      if(*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_2*NumCh*NumCyc+cutch_1*NumCyc+ncyc)>cut||*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_2*NumCh*NumCyc+cutch_2*NumCyc+ncyc)>cut){
	hit_flag++;
      }
    }//(nch_layer-1)==47
    else{
      Int_t cutch_1=nch_layer-1,cutch_2=nch_layer+1;
      if(*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_1*NumCh*NumCyc+cutch_1*NumCyc+ncyc)>cut||*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_1*NumCh*NumCyc+cutch_2*NumCyc+ncyc)>cut){
	hit_flag++;
      }
      if(*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_2*NumCh*NumCyc+cutch_1*NumCyc+ncyc)>cut||*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_2*NumCh*NumCyc+cutch_2*NumCyc+ncyc)>cut){
	hit_flag++;
      }
    }//else
  }//(ntfb-1)==0

  else if(ntfb==10){
    Int_t cuttfb_1=10,cuttfb_2=9;
    if((nch_layer-1)==0){
      Int_t cutch_1=0,cutch_2=1;
      if(*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_1*NumCh*NumCyc+cutch_1*NumCyc+ncyc)>cut||*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_1*NumCh*NumCyc+cutch_2*NumCyc+ncyc)>cut){
	hit_flag++;
      }
      if(*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_2*NumCh*NumCyc+cutch_1*NumCyc+ncyc)>cut||*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_2*NumCh*NumCyc+cutch_2*NumCyc+ncyc)>cut){
	hit_flag++;
      }
    }//(nch_layer-1)==0
    else if((nch_layer)==23){
      Int_t cutch_1=23,cutch_2=22;
      if(*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_1*NumCh*NumCyc+cutch_1*NumCyc+ncyc)>cut||*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_1*NumCh*NumCyc+cutch_2*NumCyc+ncyc)>cut){
	hit_flag++;
      }
      if(*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_2*NumCh*NumCyc+cutch_1*NumCyc+ncyc)>cut||*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_2*NumCh*NumCyc+cutch_2*NumCyc+ncyc)>cut){
	hit_flag++;
      }
    }//(nch_layer-1)==23
    else{
      Int_t cutch_1=nch_layer-1,cutch_2=nch_layer+1;
      if(*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_1*NumCh*NumCyc+cutch_1*NumCyc+ncyc)>cut||*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_1*NumCh*NumCyc+cutch_2*NumCyc+ncyc)>cut){
	hit_flag++;
      }
      if(*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_2*NumCh*NumCyc+cutch_1*NumCyc+ncyc)>cut||*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_2*NumCh*NumCyc+cutch_2*NumCyc+ncyc)>cut){
	hit_flag++;
      }
    }//else
  }
  else{
    Int_t cuttfb_1=ntfb-1,cuttfb_2=ntfb+1;
    if((nch_layer)==0){
      Int_t cutch_1=0,cutch_2=1;
      if(*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_1*NumCh*NumCyc+cutch_1*NumCyc+ncyc)>cut||*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_1*NumCh*NumCyc+cutch_2*NumCyc+ncyc)>cut){
	hit_flag++;
      }
      if(*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_2*NumCh*NumCyc+cutch_1*NumCyc+ncyc)>cut||*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_2*NumCh*NumCyc+cutch_2*NumCyc+ncyc)>cut){
	hit_flag++;
      }
    }//(nch_layer-1)==0
    else if((nch_layer)==23){
      Int_t cutch_1=23,cutch_2=22;
      if(*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_1*NumCh*NumCyc+cutch_1*NumCyc+ncyc)>cut||*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_1*NumCh*NumCyc+cutch_2*NumCyc+ncyc)>cut){
	hit_flag++;
      }
      if(*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_2*NumCh*NumCyc+cutch_1*NumCyc+ncyc)>cut||*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_2*NumCh*NumCyc+cutch_2*NumCyc+ncyc)>cut){
	hit_flag++;
      }
    }//(nch_layer-1)==23
    else{
      Int_t cutch_1=nch_layer-1,cutch_2=nch_layer+1;
      if(*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_1*NumCh*NumCyc+cutch_1*NumCyc+ncyc)>cut||*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_1*NumCh*NumCyc+cutch_2*NumCyc+ncyc)>cut){
	hit_flag++;
      }
      if(*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_2*NumCh*NumCyc+cutch_1*NumCyc+ncyc)>cut||*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_2*NumCh*NumCyc+cutch_2*NumCyc+ncyc)>cut){
	hit_flag++;
      }
    }
  }
  }

  if(nch>=24){//Y layer
  Int_t nch_layer = nch-24;
  if(ntfb==0){
    Int_t cuttfb_1=1,cuttfb_2=2;
    if((nch_layer)==0){
      Int_t cutch_1=24,cutch_2=25;
      if(*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_1*NumCh*NumCyc+cutch_1*NumCyc+ncyc)>cut||*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_1*NumCh*NumCyc+cutch_2*NumCyc+ncyc)>cut){
	hit_flag++;
      }
      if(*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_2*NumCh*NumCyc+cutch_1*NumCyc+ncyc)>cut||*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_2*NumCh*NumCyc+cutch_2*NumCyc+ncyc)>cut){
	hit_flag++;
      }
    }//(nch_layer-1)==0
    else if((nch_layer)==23){
      Int_t cutch_1=47,cutch_2=46;
      if(*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_1*NumCh*NumCyc+cutch_1*NumCyc+ncyc)>cut||*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_1*NumCh*NumCyc+cutch_2*NumCyc+ncyc)>cut){
	hit_flag++;
      }
      if(*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_2*NumCh*NumCyc+cutch_1*NumCyc+ncyc)>cut||*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_2*NumCh*NumCyc+cutch_2*NumCyc+ncyc)>cut){
	hit_flag++;
      }
    }//(nch_layer-1)==47
    else{
      Int_t cutch_1=nch_layer-1+24,cutch_2=nch_layer+1+24;
      if(*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_1*NumCh*NumCyc+cutch_1*NumCyc+ncyc)>cut||*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_1*NumCh*NumCyc+cutch_2*NumCyc+ncyc)>cut){
	hit_flag++;
      }
      if(*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_2*NumCh*NumCyc+cutch_1*NumCyc+ncyc)>cut||*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_2*NumCh*NumCyc+cutch_2*NumCyc+ncyc)>cut){
	hit_flag++;
      }
    }//else
  }//(ntfb-1)==0

  else if(ntfb==10){
    Int_t cuttfb_1=10,cuttfb_2=9;
    if((nch_layer)==0){
      Int_t cutch_1=24,cutch_2=25;
      if(*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_1*NumCh*NumCyc+cutch_1*NumCyc+ncyc)>cut||*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_1*NumCh*NumCyc+cutch_2*NumCyc+ncyc)>cut){
	hit_flag++;
      }
      if(*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_2*NumCh*NumCyc+cutch_1*NumCyc+ncyc)>cut||*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_2*NumCh*NumCyc+cutch_2*NumCyc+ncyc)>cut){
	hit_flag++;
      }
    }//(nch_layer-1)==0
    else if((nch_layer)==23){
      Int_t cutch_1=47,cutch_2=46;
      if(*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_1*NumCh*NumCyc+cutch_1*NumCyc+ncyc)>cut||*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_1*NumCh*NumCyc+cutch_2*NumCyc+ncyc)>cut){
	hit_flag++;
      }
      if(*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_2*NumCh*NumCyc+cutch_1*NumCyc+ncyc)>cut||*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_2*NumCh*NumCyc+cutch_2*NumCyc+ncyc)>cut){
	hit_flag++;
      }
    }//(nch_layer-1)==23
    else{
      Int_t cutch_1=nch_layer-1+24,cutch_2=nch_layer+1+24;
      if(*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_1*NumCh*NumCyc+cutch_1*NumCyc+ncyc)>cut||*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_1*NumCh*NumCyc+cutch_2*NumCyc+ncyc)>cut){
	hit_flag++;
      }
      if(*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_2*NumCh*NumCyc+cutch_1*NumCyc+ncyc)>cut||*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_2*NumCh*NumCyc+cutch_2*NumCyc+ncyc)>cut){
	hit_flag++;
      }
    }//else
  }
  else{
    Int_t cuttfb_1=ntfb-1,cuttfb_2=ntfb+1;
    if((nch_layer)==0){
      Int_t cutch_1=24,cutch_2=25;
      if(*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_1*NumCh*NumCyc+cutch_1*NumCyc+ncyc)>cut||*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_1*NumCh*NumCyc+cutch_2*NumCyc+ncyc)>cut){
	hit_flag++;
      }
      if(*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_2*NumCh*NumCyc+cutch_1*NumCyc+ncyc)>cut||*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_2*NumCh*NumCyc+cutch_2*NumCyc+ncyc)>cut){
	hit_flag++;
      }
    }//(nch_layer-1)==0
    else if((nch_layer)==23){
      Int_t cutch_1=47,cutch_2=46;
      if(*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_1*NumCh*NumCyc+cutch_1*NumCyc+ncyc)>cut||*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_1*NumCh*NumCyc+cutch_2*NumCyc+ncyc)>cut){
	hit_flag++;
      }
      if(*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_2*NumCh*NumCyc+cutch_1*NumCyc+ncyc)>cut||*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_2*NumCh*NumCyc+cutch_2*NumCyc+ncyc)>cut){
	hit_flag++;
      }
    }//(nch_layer-1)==23
    else{
      Int_t cutch_1=nch_layer-1+24,cutch_2=nch_layer+1+24;
      if(*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_1*NumCh*NumCyc+cutch_1*NumCyc+ncyc)>cut||*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_1*NumCh*NumCyc+cutch_2*NumCyc+ncyc)>cut){
	hit_flag++;
      }
      if(*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_2*NumCh*NumCyc+cutch_1*NumCyc+ncyc)>cut||*(pe+nmod*NumTFB*NumCh*NumCyc+cuttfb_2*NumCh*NumCyc+cutch_2*NumCyc+ncyc)>cut){
	hit_flag++;
      }
    }
  }
  }

  if(hit_flag==2){
    return 1;
  }
  else{
    return 0;
  }

}

Int_t hit::exist_Tdc(Long_t *Tdc,Int_t nmod,Int_t ntfb,Int_t nch,Int_t ncyc){
  Long_t tdc = *(Tdc+nmod*NumTFB*NumCh*NumCyc+ntfb*NumCh*NumCyc+nch*NumCyc+ncyc);

  if(tdc<Tdc_max&&tdc>Tdc_min){
    return 1;
  }
  else{
    return 0;
  }
}
