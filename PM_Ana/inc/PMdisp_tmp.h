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
  X=5+105*x+1074.5;
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
  X=5+105*x+1074.5+10;
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
  X=5+9.5+y*50+1075;
  if(r<2)R=0;
  else R=sqrt(r-2)*3;
  TArc *arc=new TArc(X,Y,R);
  arc->SetFillColor(kRed);
  arc->SetLineColor(kRed);
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
  b1->SetLineColor(kGreen);
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
  //for(j=0;j<24;j++)sci_ing(23,j*50,33,j*50+50);
  //sci_par(0,0,10,1200);

  for(k=8;k<17;k++){
    //for(j=0;j<8;j++)sci_ing(73+46*k,j*50,83+46*k,j*50+50);
    //for(j=16;j<24;j++)sci_ing(73+46*k,j*50,83+46*k,j*50+50);
    for(j=0;j<16;j++)sci_sci(71.5+46*k,400+j*25,84.5+46*k,j*25+425);
    //sci_par(50+46*k,0,60+46*k,1200);
  }
  /*
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
  */

};


void tline(int mod,int xup,int xdown,float seppen,float katamuki){
  TLine *l1=new TLine(xup,seppen+xup*katamuki, xdown, seppen+xdown*katamuki);
  l1->SetLineWidth(6);
  l1->Draw("SAME");
};



void cline(int mod, int view, int pln, float up, float down, int value, int dist1){
  TLine *l1=new TLine(zposi(mod,view,pln),up, zposi(mod,view,pln+1+dist1), down);
  //l1->SetLineWidth(sqrt(sqrt(value)*8));
  l1->SetLineWidth(3);
  //l1->SetLineColor(28);
  /*
  if(value==0)
    l1->SetLineColor(kBlue-9);
  else
  l1->SetLineColor(kBlue);
  */
  if(value==0)l1->SetLineColor(kCyan-9);
  if(value==1)l1->SetLineColor(kCyan);
  if(value==2)l1->SetLineColor(kAzure+2);
  if(value==3)l1->SetLineColor(kBlue);
  if(value==4)l1->SetLineColor(kViolet+4);
  if(value==5)l1->SetLineColor(kMagenta);
  if(value>=6)l1->SetLineColor(kMagenta+3);


  l1->SetLineColor(kBlue-9+value*2);
  /*
  if(value==0)
    l1->SetLineColor(kBlue-9);
  else
    l1->SetLineColor(kBlue-9+2);
  */
  l1->Draw("SAME");
};





void cluster(int mod,int view,int pln,int chi,int chf){
  TBox *b1=new TBox(zposi(mod,view,pln)-15,xyposi(mod,pln,chi)-15,zposi(mod,view,pln)+15,xyposi(mod,pln,chf)+15);
  //b1->SetLineColor(kViolet+1);
  b1->SetLineColor(kGray+2);
  b1->SetLineWidth(2);
  b1->SetFillStyle(0);
  b1->Draw("SAME");
}



TCanvas *c1;
bool firstdraw=true;
void EvtDisp(int mod){
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
  TH1F *v=new TH1F("","Top view",2050*10,0,830);
  //TH1F *v=new TH1F("","Top view",820*10,0,820);
  v->SetMinimum(0);
  v->SetMaximum(1200);
  v->GetXaxis()->SetLabelSize(0);
  v->GetYaxis()->SetLabelSize(0);
  v->SetStats(0);
  v->SetNdivisions(0);
  
  
  //c1->Divide(1,2);
  //c1->cd(1);
  //h->Draw();
  //drawx();
  //c1->cd(2);
  v->Draw();
  drawy();
  
  
  //*****Start drawing*****
  for(int VIEW=1;VIEW<2;VIEW++){	
    //c1->cd(VIEW+1);
    
    //*****Draw hits*****
    for(Long64_t i=0;i<hit[VIEW];i++){
      if(mod==16){
	if(VIEW==0){
	  if(pln[VIEW][i]<18)xhit(pln[VIEW][i],ch[VIEW][i],(double)pe[VIEW][i]);
	  else vhit(pln[VIEW][i]-18,ch[VIEW][i],(double)pe[VIEW][i]);
	}
	else{
	  if(pln[VIEW][i]<18)yhit(pln[VIEW][i],ch[VIEW][i],(double)pe[VIEW][i]/5);
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

        //*****Draw clusters*****                                               
        for(PLN=0;PLN<plnmax(mod);PLN++){
          for(CL=0;CL<ncl[VIEW][PLN];CL++){
            cluster(mod,VIEW,PLN,clchi[VIEW][PLN][CL],clchf[VIEW][PLN][CL]);                                                                       
          }
        }



    //*****Draw cells*****
 
    for(PLN=0;PLN<plnmax(mod)-1;PLN++){
      for(DIST=0;DIST<2;DIST++){
	if(PLN==plnmax(mod)-2&&DIST==1)continue;
	for(CELL=0;CELL<ncell[VIEW][PLN][DIST];CELL++){
	  //cline(mod,VIEW,PLN,clcenter[VIEW][PLN][cellu[VIEW][PLN][CELL][DIST]],clcenter[VIEW][PLN+1+DIST][celld[VIEW][PLN][CELL][DIST]],value[VIEW][PLN][CELL][DIST],DIST);
	}
      }
    }
 

    //*****Draw tracks*****
    for(TRA=0;TRA<ntrack[VIEW];TRA++){
      //if(ttrack[VIEW][TRA])
	tline(mod,XX[VIEW][TRA][ntracell[VIEW][TRA]],XX[VIEW][TRA][0],par[VIEW][TRA][0],par[VIEW][TRA][1]);
    }
    
    
    
    
  }
  c1->cd(0);
  c1->Update();
  
  printf("Type \'s\' to save the event display.\n");
  printf("Type \'q\' to quit.\n");
  
  char ans[8];
  fgets(ans,8,stdin);
  
  if( *ans == 's') c1->Print("tracking.eps");
  if( *ans == 'q') exit(0);
  
  c1->Clear();
  
  


};


#endif
