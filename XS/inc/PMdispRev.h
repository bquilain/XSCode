#ifndef __PMDRAW_HXX__
#define __PMDRAW_HXX__

#define dist 240  

#include <TArc.h>
#include <Hit.h>
#include<iostream>
#include<sstream>
#include<fstream>
using namespace std;
#include <iomanip>
#include <sys/stat.h>
#include <cmath>

void xhit(int x,int y,double r, vector <int> trackid){
  double X,Y,R;
  if(x==0)X=5;
  else X=46*x+9;
  if(x==0)Y=y*50+25;
  else{
  if(y<8)Y=y*50+25;
  else if(y<24)Y=412.5+(y-8)*25;
  else Y=(y-8)*50+25;
  }
  if(r<2)R=0;
  else R=sqrt(r-2)*3;
  TArc *arc=new TArc(X,Y,R);
  if(trackid.size()!=0){
    arc->SetFillColor(trackid.front());
    arc->SetLineColor(trackid.back());
  }
  else{
    arc->SetFillColor(kBlue);
    arc->SetLineColor(0);
  }
  arc->Draw("SAME");

};


void yhit(int x,int y,double r, vector <int> trackid){
  double X,Y,R;
  if(x==0)X=28;
  else X=46*x+32;
  if(x==0)Y=y*50+25;
  else{
  if(y<8)Y=y*50+25;
  else if(y<24)Y=412.5+(y-8)*25;
  else Y=(y-8)*50+25;
  }
  if(r<2)R=0;
  else R=sqrt(r-2)*3;
  TArc *arc=new TArc(X,Y,R);
  if(trackid.size()!=0){
    arc->SetFillColor(trackid.front());
    arc->SetLineColor(trackid.back());
  } 
  else{
    arc->SetFillColor(kBlue);
    arc->SetLineColor(0);
  }
  arc->Draw("SAME");
};


void vhit(int x,int y,double r, vector <int> trackid){
  double X,Y,R;
  if(x==0)Y=-55;
  else Y=1255;
  //X=9.5+y*50;
  X=5+9.5+y*50;
  if(r<2)R=0;
  else R=sqrt(r-2)*3;
  TArc *arc=new TArc(X,Y,R);
  if(trackid.size()!=0){
    arc->SetFillColor(trackid.front());
    arc->SetLineColor(trackid.back());
  } 
  else{
    arc->SetFillColor(kBlue);
    arc->SetLineColor(0);
  }
  arc->Draw("SAME");
};





void xihit(int x,int y,double r, vector <int> trackid){
  double X,Y,R;
  //X=105*x+791+dist;
  X=5+105*x+1074.5;
  Y=y*50+25;
  if(r<2)R=0;
  else R=sqrt(r-2)*3;
  TArc *arc=new TArc(X,Y,R);
  if(trackid.size()!=0){
    arc->SetFillColor(trackid.front());
    arc->SetLineColor(trackid.back());
  } 
  else{
    arc->SetFillColor(kBlue);
    arc->SetLineColor(0);
  }
  arc->Draw("SAME");

};



void yihit(int x,int y,double r, vector <int> trackid){
  double X,Y,R;
  //X=105*x+801+dist;
  X=5+105*x+1074.5+10;
  Y=y*50+25;
  if(r<2)R=0;
  else R=sqrt(r-2)*3;
  TArc *arc=new TArc(X,Y,R);
  if(trackid.size()!=0){
    arc->SetFillColor(trackid.front());
    arc->SetLineColor(trackid.back());
  } 
  else{
    arc->SetFillColor(kBlue);
    arc->SetLineColor(0);
  }
  arc->Draw("SAME");
};


void vihit(int x,int y,double r, vector <int> trackid){
  double X,Y,R;
  if(x==0)Y=-55-20;
  else Y=1255+20;
  //X=9.5+y*50+1020;
  X=5+9.5+y*50+1075;
  if(r<2)R=0;
  else R=sqrt(r-2)*3;
  TArc *arc=new TArc(X,Y,R);
  if(trackid.size()!=0){
    arc->SetFillColor(trackid.front());
    arc->SetLineColor(trackid.back());
  } 
  else{
    arc->SetFillColor(kBlue);
    arc->SetLineColor(0);
  }
  arc->Draw("SAME");
};


void sci_ing(double x,double y,double x1,double y1){
  TBox *b1=new TBox(x,y,x1,y1);
  b1->SetLineColor(kGreen);
  b1->SetLineWidth(2);
  b1->SetFillStyle(0);
  b1->Draw("SAME");
};

void sci_par(double x,double y,double x1,double y1){
  TBox *b1=new TBox(x,y,x1,y1);
  b1->SetLineColor(kYellow);
  b1->SetLineWidth(2);
  b1->SetFillStyle(0);
  b1->Draw("SAME");
};

void sci_sci(double x,double y,double x1,double y1){
  TBox *b1=new TBox(x,y,x1,y1);
  b1->SetLineColor(kGreen+2);
  b1->SetLineWidth(2);
  b1->SetFillStyle(0);
  b1->Draw("SAME");
};

void sci_veto(double x,double y,double x1,double y1){
  TBox *b1=new TBox(x,y,x1,y1);
  b1->SetLineColor(kBlue);
  b1->SetLineWidth(2);
  b1->SetFillStyle(0);
  b1->Draw("SAME");
};


void iron(double x,double y,double x1,double y1){
  TBox *b1=new TBox(x,y,x1,y1);
  b1->SetFillColor(17);
  b1->SetLineWidth(2);
  b1->Draw("SAME");
};




void drawx(void){
  
  int k,j,i;

  for(j=0;j<24;j++)sci_ing(0,j*50,10,j*50+50);
  sci_par(23,0,33,1200);
  for(k=0;k<17;k++){
    for(j=0;j<8;j++)sci_ing(50+46*k,j*50,60+46*k,j*50+50);
    for(j=16;j<24;j++)sci_ing(50+46*k,j*50,60+46*k,j*50+50);
    for(j=0;j<16;j++)sci_sci(48.5+46*k,400+j*25,61.5+46*k,j*25+425);
  sci_par(73+46*k,0,83+46*k,1200);
  }
  for(j=0;j<17;j++)sci_veto(-15.5+50*j,-60,34.5+50*j,-50);
  for(j=0;j<17;j++)sci_veto(-15.5+50*j,1250,34.5+50*j,1260);



  for(k=0;k<11;k++){
    //for(j=0;j<24;j++)sci_ing(105*k+786+dist,j*50,105*k+796+dist,j*50+50);
    for(j=0;j<24;j++)sci_ing(105*k+1074.5,j*50,105*k+1084.5,j*50+50);
  }
  for(k=0;k<9;k++)iron(105*k+816+287,0,105*k+881+287,1200);

  //for(j=0;j<22;j++)sci_veto(-15.5+50*j+1020,-60-20,34.5+50*j+1020,-50-20);
  for(j=0;j<22;j++)sci_veto(-15.5+50*j+1075,-60-20,34.5+50*j+1075,-50-20);
  //for(j=0;j<22;j++)sci_veto(-15.5+50*j+1020,1250+20,34.5+50*j+1020,1260+20);
  for(j=0;j<22;j++)sci_veto(-15.5+50*j+1075,1250+20,34.5+50*j+1075,1260+20);


  
};


void drawy(void){
  
  int k,j,i;
  for(j=0;j<24;j++)sci_ing(23,j*50,33,j*50+50);
  sci_par(0,0,10,1200);

  for(k=0;k<17;k++){
    for(j=0;j<8;j++)sci_ing(73+46*k,j*50,83+46*k,j*50+50);
    for(j=16;j<24;j++)sci_ing(73+46*k,j*50,83+46*k,j*50+50);
    for(j=0;j<16;j++)sci_sci(71.5+46*k,400+j*25,84.5+46*k,j*25+425);
  sci_par(50+46*k,0,60+46*k,1200);
  }
  for(j=0;j<17;j++)sci_veto(-15.5+50*j,-60,34.5+50*j,-50);
  for(j=0;j<17;j++)sci_veto(-15.5+50*j,1250,34.5+50*j,1260);
 


  for(k=0;k<11;k++){
    //for(j=0;j<24;j++)sci_ing(105*k+796+dist,j*50,105*k+806+dist,j*50+50);
    for(j=0;j<24;j++)sci_ing(105*k+1084.5,j*50,105*k+1094.5,j*50+50);
  }
  for(k=0;k<9;k++)iron(105*k+816+287,0,105*k+881+287,1200);

  //for(j=0;j<22;j++)sci_veto(-15.5+50*j+1020,-60-20,34.5+50*j+1020,-50-20);
  for(j=0;j<22;j++)sci_veto(-15.5+50*j+1075,-60-20,34.5+50*j+1075,-50-20);
  //for(j=0;j<22;j++)sci_veto(-15.5+50*j+1020,1250+20,34.5+50*j+1020,1260+20);
  for(j=0;j<22;j++)sci_veto(-15.5+50*j+1075,1250+20,34.5+50*j+1075,1260+20);
};


void tline(int mod,int xup,int xdown,float b,float slope){
  TLine *l1=new TLine(xup,b+xup*slope, xdown, b+xdown*slope);
  l1->SetLineWidth(5);
  l1->Draw("SAME");
};


/*
void cline(int mod, int view, int pln, float up, float down, int value, int dist1){
  TLine *l1=new TLine(zposi(mod,view,pln),up, zposi(mod,view,pln+1+dist1), down);
  //l1->SetLineWidth(sqrt(sqrt(value)*8));
  l1->SetLineWidth(2);
  //l1->SetLineColor(28);
  l1->SetLineColor(kViolet-9+value);
  l1->Draw("SAME");
};
*/


TCanvas *c1;
bool firstdraw=true;
void EvtDisp(vector <Hit3D> Vec){
  if(firstdraw){
    c1 = new TCanvas("c1","c1",2050*0.8/2*0.8,600*1.03*2*0.8);    
    firstdraw=false;
  }
  //*****Draw module*****
  //TH1F *h=new TH1F("","Side view",2050*10,0,1875);
  TH1F *h=new TH1F("","Side view",2050*10,0,1945);
  //TH1F *h=new TH1F("","Side view",820*10,0,820);
  h->SetMinimum(0);
  h->SetMaximum(1200);
  h->GetXaxis()->SetLabelSize(0);
  h->GetYaxis()->SetLabelSize(0);
  h->SetStats(0);
  h->SetNdivisions(0);
  //TH1F *v=new TH1F("","Top view",2050*10,0,1875);
  TH1F *v=new TH1F("","Top view",2050*10,0,1945);
  //TH1F *v=new TH1F("","Top view",820*10,0,820);
  v->SetMinimum(0);
  v->SetMaximum(1200);
  v->GetXaxis()->SetLabelSize(0);
  v->GetYaxis()->SetLabelSize(0);
  v->SetStats(0);
  v->SetNdivisions(0);
  
  
  c1->Divide(1,2);
  c1->cd(1);
  h->Draw();
  drawx();
  c1->cd(2);
  v->Draw();
  drawy();
  
  
  //*****Start drawing*****
 
  for(int VIEW=0;VIEW<2;VIEW++){	
    c1->cd(VIEW+1);
    
    //*****Draw hits*****
    for(Long64_t i=0;i<Vec.size();i++){
      if(Vec[i].view!=VIEW) continue;
	
      if(Vec[i].mod==16){
	if(Vec[i].pln<18 && Vec[i].pln>0 && Vec[i].ch>7 && Vec[i].ch<24) Vec[i].pecorr/=2;
	//cout<<"Mod="<<Vec[i].mod<<", pln="<<Vec[i].pln<<", view="<<Vec[i].view<<", channel="<<Vec[i].ch<<", pe="<<Vec[i].pecorr/2<<endl;
	if(Vec[i].view==0){
	  if(Vec[i].pln<18) xhit(Vec[i].pln,Vec[i].ch,(double)Vec[i].pecorr, Vec[i].RecTrk);//place le hit en x
	  else vhit(Vec[i].pln-18,Vec[i].ch,(double)Vec[i].pecorr,Vec[i].RecTrk);
	}
	else{
	  if(Vec[i].pln<18) yhit(Vec[i].pln,Vec[i].ch,(double)Vec[i].pecorr,Vec[i].RecTrk);
	  else vhit(Vec[i].pln-21,Vec[i].ch,(double)Vec[i].pecorr,Vec[i].RecTrk);
	}
      }

      else if(Vec[i].mod==3){
	if(Vec[i].view==0){
	  if(Vec[i].pln<11)xihit(Vec[i].pln,Vec[i].ch,(double)Vec[i].pecorr,Vec[i].RecTrk);
	  else vihit(Vec[i].pln-13,Vec[i].ch,(double)Vec[i].pecorr,Vec[i].RecTrk);
	}
	else{
	  if(Vec[i].pln<11)yihit(Vec[i].pln,Vec[i].ch,(double)Vec[i].pecorr,Vec[i].RecTrk);
	  else vihit(Vec[i].pln-11,Vec[i].ch,(double)Vec[i].pecorr,Vec[i].RecTrk);
	}
      }
      else continue;
    }
  

    //*****Draw cells*****
    /*
    for(PLN=0;PLN<plnmax(mod)-1;PLN++){
      for(DIST=0;DIST<2;DIST++){
	if(PLN==plnmax(mod)-2&&DIST==1)continue;
	for(CELL=0;CELL<ncell[VIEW][PLN][DIST];CELL++){
	  cline(mod,VIEW,PLN,clcenter[VIEW][PLN][cellu[VIEW][PLN][CELL][DIST]],clcenter[VIEW][PLN+1+DIST][celld[VIEW][PLN][CELL][DIST]],value[VIEW][PLN][CELL][DIST],DIST);
	}
      }
    }
    */

    //*****Draw tracks*****
    /*       for(TRA=0;TRA<ntrack[VIEW];TRA++){
      if(ttrack[VIEW][TRA])
	tline(mod,XX[VIEW][TRA][ntracell[VIEW][TRA]],XX[VIEW][TRA][0],par[VIEW][TRA][0],par[VIEW][TRA][1]);
    }
    */
    /*    TLine *l1=new TLine(Zi[VIEW],b[VIEW]+Zi[VIEW]*a[VIEW], Zf[VIEW],b[VIEW]+Zf[VIEW]*a[VIEW]);
    l1->SetLineWidth(5);
    l1->Draw("SAME");    
    */ 
  }

  c1->cd(0);
  c1->Update();
  
  printf("Type \'s\' to save the event display.\n");
  printf("Type \'q\' to quit.\n");
  
  char ans[8];
  fgets(ans,8,stdin);
  
  if( *ans == 's') c1->Print("~/Evtdisp.pdf");
  if( *ans == 'q') exit(0);
  
  c1->Clear();
  
  


};


#endif
