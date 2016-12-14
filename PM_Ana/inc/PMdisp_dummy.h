#ifndef __PMDRAW_HXX__
#define __PMDRAW_HXX__

#define dist 240  

#include <TArc.h>

void xhit(int x,int y,double r){
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
  arc->SetFillColor(kRed);
  arc->SetLineColor(kRed);
  arc->Draw("SAME");

};


void yhit(int x,int y,double r){
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
  arc->SetFillColor(kRed);
  arc->SetLineColor(kRed);
  arc->Draw("SAME");
};


void vhit(int x,int y,double r){
  double X,Y,R;
  if(x==0)Y=-55;
  else Y=1255;
  //X=9.5+y*50;
  X=5+9.5+y*50;
  if(r<2)R=0;
  else R=sqrt(r-2)*3;
  TArc *arc=new TArc(X,Y,R);
  arc->SetFillColor(kRed);
  arc->SetLineColor(kRed);
  arc->Draw("SAME");
};





void xihit(int x,int y,double r){
  double X,Y,R;
  //X=105*x+791+dist;
  X=5+105*x+1074.5-1074.5;
  Y=y*50+25;
  if(r<2)R=0;
  else R=sqrt(r-2)*3;
  TArc *arc=new TArc(X,Y,R);
  arc->SetFillColor(kRed);
  arc->SetLineColor(kRed);
  arc->Draw("SAME");

};



void yihit(int x,int y,double r){
  double X,Y,R;
  //X=105*x+801+dist;
  X=5+105*x+1074.5+10-1074.5;
  Y=y*50+25;
  if(r<2)R=0;
  else R=sqrt(r-2)*3;
  TArc *arc=new TArc(X,Y,R);
  arc->SetFillColor(kRed);
  arc->SetLineColor(kRed);
  arc->Draw("SAME");
};


void vihit(int x,int y,double r){
  double X,Y,R;
  if(x==0)Y=-55-20;
  else Y=1255+20;
  //X=9.5+y*50+1020;
  X=5+9.5+y*50+1075-1074.5;
  if(r<2)R=0;
  else R=sqrt(r-2)*3;
  TArc *arc=new TArc(X,Y,R);
  arc->SetFillColor(kRed);
  arc->SetLineColor(kRed);
  arc->Draw("SAME");
};


void mod_ing(double x,double y){
  TBox *b1=new TBox(x-60,y-60,x+60,y+60);
  b1->SetLineColor(kBlack);
  b1->SetLineWidth(1);
  b1->SetFillStyle(0);
  b1->Draw("SAME");
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

  /*
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
  */


  for(k=0;k<11;k++){
    //for(j=0;j<24;j++)sci_ing(105*k+786+dist,j*50,105*k+796+dist,j*50+50);
    for(j=0;j<24;j++)sci_ing(105*k+1074.5-1074.5,j*50,105*k+1084.5-1074.5,j*50+50);
  }
  for(k=0;k<9;k++)iron(105*k+816+287-1074.5,0,105*k+881+287-1074.5,1200);

  //for(j=0;j<22;j++)sci_veto(-15.5+50*j+1020,-60-20,34.5+50*j+1020,-50-20);
  for(j=0;j<22;j++)sci_veto(-15.5+50*j+1075-1074.5,-60-20,34.5+50*j+1075-1074.5,-50-20);
  //for(j=0;j<22;j++)sci_veto(-15.5+50*j+1020,1250+20,34.5+50*j+1020,1260+20);
  for(j=0;j<22;j++)sci_veto(-15.5+50*j+1075-1074.5,1250+20,34.5+50*j+1075-1074.5,1260+20);


  
};


void drawy(void){

  int k,j,i;
  /*
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
  */


  for(k=0;k<11;k++){
    //for(j=0;j<24;j++)sci_ing(105*k+796+dist,j*50,105*k+806+dist,j*50+50);
    for(j=0;j<24;j++)sci_ing(105*k+1084.5-1074.5,j*50,105*k+1094.5-1074.5,j*50+50);
  }
  for(k=0;k<9;k++)iron(105*k+816+287-1074.5,0,105*k+881+287-1074.5,1200);

  //for(j=0;j<22;j++)sci_veto(-15.5+50*j+1020,-60-20,34.5+50*j+1020,-50-20);
  for(j=0;j<22;j++)sci_veto(-15.5+50*j+1075-1074.5,-60-20,34.5+50*j+1075-1074.5,-50-20);
  //for(j=0;j<22;j++)sci_veto(-15.5+50*j+1020,1250+20,34.5+50*j+1020,1260+20);
  for(j=0;j<22;j++)sci_veto(-15.5+50*j+1075-1074.5,1250+20,34.5+50*j+1075-1074.5,1260+20);
};


void tline(int mod,int xup,int xdown,float seppen,float katamuki){
  TLine *l1=new TLine(xup-1074.5,seppen+xup*katamuki, xdown-1074.5, seppen+xdown*katamuki);
  l1->SetLineWidth(5);
  l1->Draw("SAME");
};



void cline(int mod, int view, int pln, float up, float down, int value, int dist1){
  TLine *l1=new TLine(zposi(mod,view,pln),up, zposi(mod,view,pln+1+dist1), down);
  //l1->SetLineWidth(sqrt(sqrt(value)*8));
  l1->SetLineWidth(2);
  //l1->SetLineColor(28);
  l1->SetLineColor(kViolet-9+value);
  l1->Draw("SAME");
};



void Timing(int i, int bmax){
  Int_t beam =2830;
  Int_t bwspill =  581;
  //TBox *b1=new TBox(2780+i*(577.14+3),0,2860+i*(577.14+3),1);              
  TBox *b1=new TBox(beam+bwspill*i-53,0,beam+bwspill*i+53,bmax);
  b1->SetFillColor(kCyan);
  b1->SetLineWidth(2);
  b1->SetFillStyle(3354);
  b1->Draw("SAME");
 TBox *b2=new TBox(beam+bwspill*i-53,0,beam+bwspill*i+53,bmax);
  b2->SetLineColor(kCyan);
  b2->SetLineWidth(1);
  b2->SetFillStyle(0);
  b2->Draw("SAME");
};



TCanvas *c1;
bool firstdraw=true;
void EvtDisp(int mod){





  if(firstdraw){
    c1 = new TCanvas("c1","c1",2050*1.4/2*0.88*1.08*1.006,600*1.03*2*0.8*1.05*1.006);    
    firstdraw=false;
  }
  //*****Draw module*****
  //TH1F *h=new TH1F("","Side view",2050*10,0,1875);
  TH1F *h=new TH1F("","Side view",930*10,0,1057);
  //TH1F *h=new TH1F("","Side view",820*10,0,820);
  h->SetMinimum(0);
  h->SetMaximum(1200);
  h->GetXaxis()->SetLabelSize(0);
  h->GetYaxis()->SetLabelSize(0);
  h->GetYaxis()->SetAxisColor(0);
  h->GetXaxis()->SetAxisColor(0);
  h->SetStats(0);
  h->SetNdivisions(0);
  //TH1F *v=new TH1F("","Top view",2050*10,0,1875);
  TH1F *v=new TH1F("","Top view",930*10,0,1057);
  //TH1F *v=new TH1F("","Top view",820*10,0,820);
  v->SetMinimum(0);
  v->SetMaximum(1200);
  v->GetXaxis()->SetLabelSize(0);
  v->GetYaxis()->SetLabelSize(0);
  v->GetYaxis()->SetAxisColor(0);
  v->GetXaxis()->SetAxisColor(0);
  v->SetStats(0);
  v->SetNdivisions(0);
  
  
  c1->Divide(1,2);
  c1->cd(1)->Divide(3,1);
  c1->cd(2)->Divide(2,1);

  c1->cd(1)->cd(1);
  h->Draw();
  drawx();
  c1->cd(1)->cd(2);
  v->Draw();
  drawy();
  
  
  //*****Start drawing*****
  for(int VIEW=0;VIEW<2;VIEW++){	
    c1->cd(1)->cd(VIEW+1);
    
    //*****Draw hits*****
    for(Long64_t i=0;i<hit[VIEW];i++){
      if(mod==16){
	if(VIEW==0){
	  if(pln[VIEW][i]<18)xhit(pln[VIEW][i],ch[VIEW][i],(double)pe[VIEW][i]);
	  else vhit(pln[VIEW][i]-18,ch[VIEW][i],(double)pe[VIEW][i]);
	}
	else{
	  if(pln[VIEW][i]<18)yhit(pln[VIEW][i],ch[VIEW][i],(double)pe[VIEW][i]);
	  else vhit(pln[VIEW][i]-21,ch[VIEW][i],(double)pe[VIEW][i]);
	}
      }
      else{
	if(VIEW==0){
	  if(pln[VIEW][i]<11)xihit(pln[VIEW][i],ch[VIEW][i],(double)pe[VIEW][i]);
	  else vihit(pln[VIEW][i]-13,ch[VIEW][i],(double)pe[VIEW][i]);
	}
	else{
	  if(pln[VIEW][i]<11)yihit(pln[VIEW][i],ch[VIEW][i],(double)pe[VIEW][i]);
	  else vihit(pln[VIEW][i]-11,ch[VIEW][i],(double)pe[VIEW][i]);
	}
      }      
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
    for(TRA=0;TRA<ntrack[VIEW];TRA++){
      if(ttrack[VIEW][TRA])
	tline(mod,XX[VIEW][TRA][ntracell[VIEW][TRA]],XX[VIEW][TRA][0],par[VIEW][TRA][0],par[VIEW][TRA][1]);
    }
    
    
    
    
  }



  c1->cd(1)->cd(3);
  TH1F *m=new TH1F("","Upstream view",600,-600,600);
  //TH1F *h=new TH1F("","Side view",820*10,0,820);
  m->SetMinimum(0);
  m->SetMaximum(1300);
  m->GetXaxis()->SetLabelSize(0);
  m->GetYaxis()->SetLabelSize(0);
  m->GetYaxis()->SetAxisColor(0);
  m->GetXaxis()->SetAxisColor(0);
  m->SetStats(0);
  m->SetNdivisions(0);
  m->Draw();

  for(int ib=0;ib<7;ib++)
    mod_ing(-(ib-3)*150,600);
  for(int ib=0;ib<7;ib++)
    mod_ing(0,(ib-3)*150+600);

  if(mod<7){
    float xm=-(mod-3)*150;
    TBox *b2=new TBox(xm-60,600-60,xm+60,600+60);
    b2->SetLineColor(kBlack);
    b2->SetFillColor(kRed);
    b2->SetLineWidth(1);
    b2->SetFillStyle(1);
    b2->Draw("SAME");
  }
  
  else if(mod<14){
    float ym=(mod-10)*150;
    TBox *b2=new TBox(-60,ym+600-60,60,ym+600+60);
    b2->SetLineColor(kBlack);
    b2->SetFillColor(kRed);
    b2->SetLineWidth(1);
    b2->SetFillStyle(1);
    b2->Draw("SAME");

  }


  TLatex *tx0=new TLatex(150*3-30,600-30,"0");
  tx0->SetTextSize(0.08);
  tx0->Draw("SAME");
  TLatex *tx1=new TLatex(150*2-30,600-30,"1");
  tx1->SetTextSize(0.08);
  tx1->Draw("SAME");
  TLatex *tx2=new TLatex(150*1-30,600-30,"2");
  tx2->SetTextSize(0.08);
  tx2->Draw("SAME");
  TLatex *tx6=new TLatex(-150*3-30,600-30,"6");
  tx6->SetTextSize(0.08);
  tx6->Draw("SAME");
  TLatex *tx5=new TLatex(-150*2-30,600-30,"5");
  tx5->SetTextSize(0.08);
  tx5->Draw("SAME");
  TLatex *tx4=new TLatex(-150*1-30,600-30,"4");
  tx4->SetTextSize(0.08);
  tx4->Draw("SAME");
  TLatex *tx7=new TLatex(-30,-150*3+600-30,"7");
  tx7->SetTextSize(0.08);
  tx7->Draw("SAME");
  TLatex *tx8=new TLatex(-30,-150*2+600-30,"8");
  tx8->SetTextSize(0.08);
  tx8->Draw("SAME");
  TLatex *tx9=new TLatex(-30,-150*1+600-30,"9");
  tx9->SetTextSize(0.08);
  tx9->Draw("SAME");
  TLatex *tx13=new TLatex(-58,150*3+600-30,"13");
  tx13->SetTextSize(0.08);
  tx13->Draw("SAME");
  TLatex *tx12=new TLatex(-58,150*2+600-30,"12");
  tx12->SetTextSize(0.08);
  tx12->Draw("SAME");
  TLatex *tx11=new TLatex(-55,150*1+600-30,"11");
  tx11->SetTextSize(0.08);
  tx11->Draw("SAME");

  if(mod!=10&&mod!=3){
    TLatex *tx10=new TLatex(-30+8-9.5,+600-30+28,"10");
    tx10->SetTextSize(0.06);
    tx10->Draw("SAME");
    TLatex *tx3=new TLatex(-30-28+2,+600-30-28,"3");
    tx3->SetTextSize(0.06);
    tx3->Draw("SAME");
  }
  else if(mod==3){
    TLatex *tx3=new TLatex(-30,600-30,"3");
    tx3->SetTextSize(0.08);
    tx3->Draw("SAME");
  }
  else if(mod==10){
    TLatex *tx10=new TLatex(-58,600-30,"10");
    tx10->SetTextSize(0.08);
    tx10->Draw("SAME");
  }


  TLatex *tit=new TLatex(-150*2-30,1170,"Hit Module");
  tit->SetTextSize(0.08);
  tit->SetTextColor(kRed);
  tit->Draw("SAME");


  TLatex *ts0=new TLatex(-150*3-142,220,"0~6:");
  ts0->SetTextSize(0.042);
  ts0->Draw("SAME");
  TLatex *ts1=new TLatex(-150*3-142,170,"Horizontal Module");
  ts1->SetTextSize(0.042);
  ts1->Draw("SAME");
  TLatex *ts2=new TLatex(-150*3-142,220-150,"7~13:");
  ts2->SetTextSize(0.042);
  ts2->Draw("SAME");
  TLatex *ts3=new TLatex(-150*3-142,170-150,"Vertical Module");
  ts3->SetTextSize(0.042);
  ts3->Draw("SAME");



  c1->cd(2)->cd(1);
  TH1F *hp=new TH1F("","P.E. distribution",60,0,120);
  hp->SetFillColor(kRed);
  hp->SetFillStyle(3002);
  hp->GetXaxis()->SetLabelSize(0.05);
  hp->GetYaxis()->SetLabelSize(0.05);
  hp->GetXaxis()->SetTitleSize(0.05);
  hp->GetYaxis()->SetTitleSize(0.05);
  hp->GetXaxis()->SetTitle("Number of P.E.");
  hp->GetYaxis()->SetTitle("Number of hits");
  hp->SetStats(0);


  TH1F *ht=new TH1F("","Timing distribution",600,2000,8000);
  ht->SetFillColor(kBlack);
  //ht->SetFillStyle(1);
  ht->GetXaxis()->SetLabelSize(0.05);
  ht->GetYaxis()->SetLabelSize(0.05);
  ht->GetXaxis()->SetTitleSize(0.05);
  ht->GetYaxis()->SetTitleSize(0.05);
  ht->GetXaxis()->SetTitle("Hit timing [nsec]");
  ht->GetYaxis()->SetTitle("Number of hits");
  ht->SetStats(0);


  for(int VIEW=0;VIEW<2;VIEW++){	
    for(Long64_t i=0;i<hit[VIEW];i++){
      hp->Fill(pe[VIEW][i]);
      ht->Fill(timin[VIEW][i]);
      //cout<<timin[VIEW][i]<<endl;
    }
  }
  cout<<"Timing     : "<<ht->GetMean()<<" ns"<<endl;
  cout<<"Mean P.E.  : "<<hp->GetMean()<<endl;

  hp->Draw();

  c1->cd(2)->cd(2);

  int binmax=ht->GetBinContent(ht->GetMaximumBin());
  TH1F *ht2=new TH1F("","Timing distribution",600,2000,8000);
  ht2->SetFillColor(kBlack);
  ht2->GetXaxis()->SetLabelSize(0.05);
  ht2->GetYaxis()->SetLabelSize(0.05);
  ht2->GetXaxis()->SetTitleSize(0.05);
  ht2->GetYaxis()->SetTitleSize(0.05);
  ht2->GetXaxis()->SetTitle("Hit timing [nsec]");
  ht2->GetYaxis()->SetTitle("Number of hits");
  ht2->SetStats(0);
  ht2->SetMinimum(0);
  ht2->SetMaximum(binmax+1);
  ht2->Draw();

  for(Long64_t index=0;index<8;index++){
    Timing(index,binmax+1);
  }

  ht->Draw("SAME");

  TLatex *tt=new TLatex(4500,(binmax+1)*1.02,"Expected timing");
  tt->SetTextSize(0.06);
  tt->SetTextColor(kCyan);
  tt->Draw("SAME");



  c1->cd(0);
  c1->Update();
  
  //printf("Type \'s\' to save the event display.\n");
  //printf("Type \'q\' to quit.\n");
  
  //char ans[8];
  //fgets(ans,8,stdin);
  
  printf("\a");fflush(stdout);


  usleep(10*1000);


  //if( *ans == 's') c1->Print("Evtdisp.pdf");
  //if( *ans == 'q') exit(0);
  
  c1->Clear();
  
  


};





#endif
