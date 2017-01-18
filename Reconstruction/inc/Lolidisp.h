#ifndef __LOLIDRAW_HXX__
#define __LOLIDRAW_HXX__

#define dist 240  

#include <TArc.h>
#include "INGRID_Dimension.cxx"
#include "LoliAna.hxx"

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


void xloliparticle(){
//starting point
  double X,Y,R;
  X = ( simpar_for_disp->ipos[2] + 143.3 )*10.+409.5-240;
  Y = simpar_for_disp->ipos[1] * 10. + 600;
  std::cout << X << " " << Y << "\n";
  R=10.;
  TArc *arc=new TArc(X,Y,R);
  arc->SetFillColor(6);
  arc->SetLineColor(6);
  arc->Draw("SAME");
 
//end point
  double X2,Y2,R2;
  X2 = ( simpar_for_disp->fpos[2] + 143.3 )*10.+409.5-240;
  Y2 = simpar_for_disp->fpos[1] * 10. + 600;
  std::cout << X2 << " " << Y2 << "\n";
  R2=10.;
  TArc *arc2=new TArc(X2,Y2,R2);
  arc2->SetFillColor(6);
  arc2->SetLineColor(6);
  arc2->Draw("SAME");

//line
  TLine *l1=new TLine(X,Y,X2,Y2);
  l1->SetLineWidth(2);
  l1->SetLineColor(6);
  l1->Draw("SAME");

};

void yloliparticle(){
//starting point
  double X,Y,R;
  X = ( simpar_for_disp->ipos[2] + 143.3 )*10.+409.5-240.;
  Y = simpar_for_disp->ipos[0] * 10. + 600;

  std::cout << X << " " << Y << "\n";
  R=10.;
  TArc *arc=new TArc(X,Y,R);
  arc->SetFillColor(6);
  arc->SetLineColor(6);
  arc->Draw("SAME");
 
//end point
  double X2,Y2,R2;
  X2 = ( simpar_for_disp->fpos[2] + 143.3 )*10.+409.5-240.;
  Y2 = simpar_for_disp->fpos[0] * 10. + 600;
  std::cout << X2 << " " << Y2 << "\n";
  R2=10.;
  TArc *arc2=new TArc(X2,Y2,R2);
  arc2->SetFillColor(6);
  arc2->SetLineColor(6);
  arc2->Draw("SAME");

//line
  TLine *l1=new TLine(X,Y,X2,Y2);
  l1->SetLineWidth(2);
  l1->SetLineColor(6);
  l1->Draw("SAME");

};





void xlolihit(){
  
  int mod,view,pln,ch;
  double pe,X,Y,R,pe_cross;

  INGRID_Dimension *fdim = new INGRID_Dimension();
  for(int i=0; i<(int)allhit_for_disp.size(); i++){
	mod =allhit_for_disp[i].mod;
	view=allhit_for_disp[i].view;
	pln =allhit_for_disp[i].pln;
	ch  =allhit_for_disp[i].ch;
	pe  =allhit_for_disp[i].pe;
	pe_cross  =allhit_for_disp[i].pe_cross;
	if(view==1)continue;
	if(mod==15){
  		fdim -> get_pos_loli_xy( mod, view, pln, ch, &Y, &X);
		X=10.*X+409.5;//cm->mm
		Y=10.*Y+600;//cm->mm
	}
	else if(mod==3){
  		fdim -> get_posXY( mod, view, pln, ch, &Y, &X);
		X=10.*X + 1074.5;//cm->mm
		Y=10.*Y + 25;//cm->mm
	}
	else continue;
	std::cout << "hit " << "mod" << mod << " view" << view << " pln" << pln << " ch" << ch << " X" << X << " Y" << Y << " pe" << pe << " pe_cross" << pe_cross << " id" << allhit_for_disp[i].id << endl;
  	//if(pe<1.5)R=0;
	//else R=sqrt(pe-1)*1;
	//R=sqrt(pe)*1;
	R=sqrt(pe)*2;
	TArc *arc=new TArc(X,Y,R);
	arc->SetFillColor(kRed);
	arc->SetLineColor(kRed);
	arc->Draw("SAME");

	//crosstalk
  	//if(pe_cross<1.5)R=0;
	//else R=sqrt(pe_cross-1)*1;
	//R=sqrt(pe_cross)*2;
	//TArc *arc_cross=new TArc(X,Y,R);
	//arc_cross->SetFillColor(kBlue);
	//arc_cross->SetLineColor(kBlue);
	//arc_cross->Draw("SAME");

//	//pe+crosstalk
//	R=sqrt(pe+pe_cross)*2;
//	TArc *arc=new TArc(X,Y,R);
//	arc->SetFillColor(kRed);
//	arc->SetLineColor(kRed);
//	arc->Draw("SAME");
  }
  delete fdim;


};

void ylolihit(){
  
  int mod,view,pln,ch;
  double pe,X,Y,R,pe_cross;

  INGRID_Dimension *fdim = new INGRID_Dimension();
  for(int i=0; i<(int)allhit_for_disp.size(); i++){
	mod =allhit_for_disp[i].mod;
	view=allhit_for_disp[i].view;
	pln =allhit_for_disp[i].pln;
	ch  =allhit_for_disp[i].ch;
	pe  =allhit_for_disp[i].pe;
	pe_cross  =allhit_for_disp[i].pe_cross;
	if(view==0)continue;
	if(mod==15){
	  	fdim -> get_pos_loli_xy( mod, view, pln, ch, &Y, &X);
		X=10.*X+409.5;//cm->mm
		Y=10.*Y+600;//cm->mm
	}
	else if(mod==3){
  		fdim -> get_posXY( mod, view, pln, ch, &Y, &X);
		X=10.*X + 1074.5;//cm->mm
		Y=10.*Y + 25;//cm->mm
	}
	else continue;
	std::cout << "hit " << "mod" << mod << " view" << view << " pln" << pln << " ch" << ch << " X" << X << " Y" << Y << " pe" << pe << " pe_cross" << pe_cross << "\n";
	//if(pe<1.5)R=0;
	//else R=sqrt(pe-1)*2;
	R=sqrt(pe)*2;
	TArc *arc=new TArc(X,Y,R);
	arc->SetFillColor(kRed);
	arc->SetLineColor(kRed);
	arc->Draw("SAME");

	//crosstalk
  	//if(pe_cross<1.5)R=0;
	//else R=sqrt(pe_cross-1)*1;
	//R=sqrt(pe_cross)*2;
	//TArc *arc_cross=new TArc(X,Y,R);
	//arc_cross->SetFillColor(kBlue);
	//arc_cross->SetLineColor(kBlue);
	//arc_cross->Draw("SAME");

//	//pe+crosstalk
//	R=sqrt(pe+pe_cross)*2;
//	TArc *arc=new TArc(X,Y,R);
//	arc->SetFillColor(kRed);
//	arc->SetLineColor(kRed);
//	arc->Draw("SAME");


  }
  delete fdim;


};



void xloli_reconhit(){
  
  int mod,view,pln,ch;
  double pe,X,Y,R,pe_cross;
  int id;

  INGRID_Dimension *fdim = new INGRID_Dimension();
  for(int i=0; i<alltrack.size(); i++){
    for(int j=0; j<alltrack[i].hitid.size(); j++){
        id=-1;
        for(int k=0; k<allhit_for_disp.size(); k++){
  	  if(alltrack[i].hitid[j] == allhit_for_disp[k].id){
		id=k;
		continue;
	  }
        }
        if(id==-1)std::cout << "hit id is not found" << std::endl;
	mod =allhit_for_disp[id].mod;
	view=allhit_for_disp[id].view;
	pln =allhit_for_disp[id].pln;
	ch  =allhit_for_disp[id].ch;
	pe  =allhit_for_disp[id].pe;
	pe_cross  =allhit_for_disp[id].pe_cross;
	if(view==1)continue;
	if(mod==15){
  		fdim -> get_pos_loli_xy( mod, view, pln, ch, &Y, &X);
		X=10.*X+409.5;//cm->mm
		Y=10.*Y+600;//cm->mm
	}
	else if(mod==3){
  		fdim -> get_posXY( mod, view, pln, ch, &Y, &X);
		X=10.*X + 1074.5;//cm->mm
		Y=10.*Y + 25;//cm->mm
	}
	else continue;
	std::cout << "reconstructed hit " << "mod" << mod << " view" << view << " pln" << pln << " ch" << ch << " X" << X << " Y" << Y << " pe" << pe << " pe_cross" << pe_cross << "\n";
  	//if(pe<1.5)R=0;
	//else R=sqrt(pe-1)*1;
	//R=sqrt(pe)*1;
	R=sqrt(pe)*2;
	TArc *arc=new TArc(X,Y,R);
	arc->SetFillColor(kBlue);
	arc->SetLineColor(kBlue);
	arc->Draw("SAME");

	//crosstalk
  	//if(pe_cross<1.5)R=0;
	//else R=sqrt(pe_cross-1)*1;
	//R=sqrt(pe_cross)*2;
	//TArc *arc_cross=new TArc(X,Y,R);
	//arc_cross->SetFillColor(kBlue);
	//arc_cross->SetLineColor(kBlue);
	//arc_cross->Draw("SAME");

//	//pe+crosstalk
//	R=sqrt(pe+pe_cross)*2;
//	TArc *arc=new TArc(X,Y,R);
//	arc->SetFillColor(kRed);
//	arc->SetLineColor(kRed);
//	arc->Draw("SAME");
  }//j
  }//i
  delete fdim;


};



void yloli_reconhit(){
  
  int mod,view,pln,ch;
  double pe,X,Y,R,pe_cross;
  int id;

  INGRID_Dimension *fdim = new INGRID_Dimension();
  for(int i=0; i<alltrack.size(); i++){
    for(int j=0; j<alltrack[i].hitid.size(); j++){
        id=-1;
        for(int k=0; k<allhit_for_disp.size(); k++){
  	  if(alltrack[i].hitid[j] == allhit_for_disp[k].id){
		id=k;
		continue;
	  }
        }
	mod =allhit_for_disp[id].mod;
	view=allhit_for_disp[id].view;
	pln =allhit_for_disp[id].pln;
	ch  =allhit_for_disp[id].ch;
	pe  =allhit_for_disp[id].pe;
	pe_cross  =allhit_for_disp[id].pe_cross;
	if(view==0)continue;
	if(mod==15){
	  	fdim -> get_pos_loli_xy( mod, view, pln, ch, &Y, &X);
		X=10.*X+409.5;//cm->mm
		Y=10.*Y+600;//cm->mm
	}
	else if(mod==3){
  		fdim -> get_posXY( mod, view, pln, ch, &Y, &X);
		X=10.*X + 1074.5;//cm->mm
		Y=10.*Y + 25;//cm->mm
	}
	else continue;
	std::cout << "reconstructed hit " << "mod" << mod << " view" << view << " pln" << pln << " ch" << ch << " X" << X << " Y" << Y << " pe" << pe << " pe_cross" << pe_cross << "\n";
	//if(pe<1.5)R=0;
	//else R=sqrt(pe-1)*2;
	R=sqrt(pe)*2;
	TArc *arc=new TArc(X,Y,R);
	arc->SetFillColor(kBlue);
	arc->SetLineColor(kBlue);
	arc->Draw("SAME");

	//crosstalk
  	//if(pe_cross<1.5)R=0;
	//else R=sqrt(pe_cross-1)*1;
	//R=sqrt(pe_cross)*2;
	//TArc *arc_cross=new TArc(X,Y,R);
	//arc_cross->SetFillColor(kBlue);
	//arc_cross->SetLineColor(kBlue);
	//arc_cross->Draw("SAME");

//	//pe+crosstalk
//	R=sqrt(pe+pe_cross)*2;
//	TArc *arc=new TArc(X,Y,R);
//	arc->SetFillColor(kRed);
//	arc->SetLineColor(kRed);
//	arc->Draw("SAME");

  }//j
  }//i
  delete fdim;


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

void water(double x,double y,double x1,double y1){
  TBox *b1=new TBox(x,y,x1,y1);
  b1->SetLineColor(kBlue);
  b1->SetLineWidth(2);
  b1->SetFillStyle(0);
  b1->Draw("SAME");
};



void drawx(void){
  
  int k,j;

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
  
  int k,j;
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



void drawlolix(void){

	int mod,view,pln,ch,axis;
	double posixy,posiz,thick,width;

	//mod=15; view=0; axis=0;
	mod=15; view=0; axis=1;
	for(pln=0;pln<plnmax(mod,view,pln,axis);pln++){
	for(ch=0; ch <chmax(mod,view,pln,axis);ch++){
		posiz =zposi(mod,view,pln,ch,axis);
		posixy=xyposi(mod,view,pln,ch,axis);
		thick =scithick(mod,view,pln,0,axis);
		width =sciwidth(mod,view,pln,ch,axis);
    		//sci_ing(posiz-thick,posixy-width,posiz+thick,posixy+width);
    		sci_ing(posixy-width,posiz-thick,posixy+width,posiz+thick);
	}
	}

	mod=3; view=0; axis=0;
	for(pln=0;pln<plnmax(mod,view,pln,axis);pln++){
	for(ch=0; ch <chmax(mod,view,pln,axis);ch++){
		posiz =zposi(mod,view,pln,ch,axis);
		posixy=xyposi(mod,view,pln,ch,axis);
		thick =scithick(mod,view,pln,0,axis);
		width =sciwidth(mod,view,pln,ch,axis);
    		sci_ing(posiz-thick,posixy-width,posiz+thick,posixy+width);
	}
	}
	//water
	water(409.5 -233,-30,409 +233,1230);
	//iron
  	for(int k=0;k<9;k++)iron(105*k+816+287,0,105*k+881+287,1200);

};

void drawloliy(void){

	int mod,view,pln,ch,axis;
	double posixy,posiz,thick,width;

	//mod=15; view=1; axis=0;
	mod=15; view=1; axis=1;
	for(pln=0;pln<plnmax(mod,view,pln,axis);pln++){
	for(ch=0; ch <chmax(mod,view,pln,axis);ch++){
		posiz =zposi(mod,view,pln,ch,axis);
		posixy=xyposi(mod,view,pln,ch,axis);
		thick =scithick(mod,view,pln,0,axis);
		width =sciwidth(mod,view,pln,ch,axis);
    		//sci_ing(posiz-thick,posixy-width,posiz+thick,posixy+width);
    		sci_ing(posixy-width,posiz-thick,posixy+width,posiz+thick);
	}
	}

	mod=3; view=1; axis=0;
	for(pln=0;pln<plnmax(mod,view,pln,axis);pln++){
	for(ch=0; ch <chmax(mod,view,pln,axis);ch++){
		posiz =zposi(mod,view,pln,ch,axis);
		posixy=xyposi(mod,view,pln,ch,axis);
		thick =scithick(mod,view,pln,0,axis);
		width =sciwidth(mod,view,pln,ch,axis);
    		sci_ing(posiz-thick,posixy-width,posiz+thick,posixy+width);
	}
	}
	//water
	water(409.5 -233,-30,409 +233,1230);
	//iron
  	for(int k=0;k<9;k++)iron(105*k+816+287,0,105*k+881+287,1200);

};





void tline(int mod,int xup,int xdown,float seppen,float katamuki){
  TLine *l1=new TLine(xup,seppen+xup*katamuki, xdown, seppen+xdown*katamuki);
  l1->SetLineWidth(5);
  l1->Draw("SAME");
};


void tline_ana(int mod,float xup,float xdown,float seppen,float katamuki){
  TLine *l1=new TLine(xup,seppen, xdown, seppen+(xdown-xup)*katamuki);
  l1->SetLineWidth(5);
  l1->Draw("SAME");
};


void cline(int mod, int view, int pln, float up, float down, int value, int dist1, int axis){
  //TLine *l1=new TLine(zposi(mod,view,pln),up, zposi(mod,view,pln+1+dist1), down);
  //TLine *l1=new TLine(zposi(mod,view,pln,0,axis),up, zposi(mod,view,pln+1+dist1,0,axis), down);
  TLine *l1=new TLine(up,zposi(mod,view,pln,0,axis), down,zposi(mod,view,pln+1+dist1,0,axis));
  //TLine *l1=new TLine(zposi(mod,view,pln,0,0),up, zposi(mod,view,pln+1+dist1,0,0), down);
  //l1->SetLineWidth(sqrt(sqrt(value)*8));
  l1->SetLineWidth(2);
  //l1->SetLineColor(28);
  l1->SetLineColor(kViolet-9+value);
  l1->Draw("SAME");
};

void xloliline(){
	for(int i=0;i<(int)alltrack.size();i++){
		if(alltrack[i].view==0){
			TLine* l1 = new TLine(alltrack[i].iz, alltrack[i].ixy, alltrack[i].fz, alltrack[i].fxy);
  			l1->SetLineWidth(2);
	  		l1->Draw("SAME");
			std::cout << "track: iz " << alltrack[i].iz << ", ixy " <<  alltrack[i].ixy << ", fz " <<  alltrack[i].fz << ", fxy " <<  alltrack[i].fxy << ", slope  " << alltrack[i].slope << ", ang " << alltrack[i].ang << ", intcpt " << alltrack[i].intcpt << "\n";
		}
	}
}

void yloliline(){
	for(int i=0;i<(int)alltrack.size();i++){
		if(alltrack[i].view==1){
			TLine* l1 = new TLine(alltrack[i].iz, alltrack[i].ixy, alltrack[i].fz, alltrack[i].fxy);
  			l1->SetLineWidth(2);
	  		l1->Draw("SAME");
			std::cout << "track: iz " << alltrack[i].iz << ", ixy " <<  alltrack[i].ixy << ", fz " <<  alltrack[i].fz << ", fxy " <<  alltrack[i].fxy << ", slope " << alltrack[i].slope << ", ang " << alltrack[i].ang << ", intcpt " << alltrack[i].intcpt << "\n";
		}
	}
}

TCanvas *c1;
bool firstdraw=true;

void EvtDisp(int mod, int ievt=0, bool batch=false){
  if(batch) gROOT->SetBatch();
  cout << "Start EvtDisp" << endl;
  if(firstdraw){
    c1 = new TCanvas("c1","c1",(int)(2050*0.8/2*0.8),(int)(600*1.03*2*0.8));
    firstdraw=false;
  }
  //*****Draw module*****
  TH1F *h=new TH1F("","Side view",2050*10,0,1945);
  h->SetMinimum(0);
  h->SetMaximum(1200);
  h->GetXaxis()->SetLabelSize(0);
  h->GetYaxis()->SetLabelSize(0);
  h->SetStats(0);
  h->SetNdivisions(0);
  TH1F *v=new TH1F("","Top view",2050*10,0,1945);
  v->SetMinimum(0);
  v->SetMaximum(1200);
  v->GetXaxis()->SetLabelSize(0);
  v->GetYaxis()->SetLabelSize(0);
  v->SetStats(0);
  v->SetNdivisions(0);
  
  
  c1->Divide(1,2);
  c1->cd(1);
  h->Draw();
  //cout << "after h->Draw()" << endl;
  if(mod==15 || mod==3){
  	drawlolix();
  }
  else if(mod==16 || mod<15){
  	drawx();
  }
  c1->cd(2);
  v->Draw();
  if(mod==15 || mod==3){
  	drawloliy();
  }
  else if(mod==16 || mod<15){
  	drawy();
  }
  
  //*****Start drawing*****
  for(int VIEW=0;VIEW<2;VIEW++){	
    c1->cd(VIEW+1);
   
    //*****Draw Vertex and Muon***** //simulation
    //if(mod==15){
    //      if(VIEW==0)xloliparticle();
    //      if(VIEW==1)yloliparticle();
    //}
 
    //*****Draw hits*****
    if(mod==15 || mod==3){
          if(VIEW==0)xlolihit();
          else if(VIEW==1)ylolihit();
    }
    else{

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

    }

//    //*****Draw cells*****
//    
//    for(PLN=0;PLN<plnmax(mod,VIEW,PLN,1)-1;PLN++){
//    //for(PLN=0;PLN<plnmax(mod,VIEW,PLN,0)-1;PLN++){
//      for(DIST=0;DIST<6;DIST++){
//	//if(PLN==plnmax(mod)-2&&DIST==1)continue;
//	if(PLN==plnmax(mod,VIEW,PLN,1)-DIST)continue;
//	//if(PLN==plnmax(mod,VIEW,PLN,0)-2&&DIST==1)continue;
//	for(CELL=0;CELL<ncell[VIEW][PLN][DIST];CELL++){
//	  cline(mod,VIEW,PLN,clcenter[VIEW][PLN][cellu[VIEW][PLN][CELL][DIST]],clcenter[VIEW][PLN+1+DIST][celld[VIEW][PLN][CELL][DIST]],value[VIEW][PLN][CELL][DIST],DIST,1);
//	}
//      }
//    }
    

//    //*****Draw tracks*****
//    //std::cout << ntrack[VIEW] << "\n";
//    for(TRA=0;TRA<ntrack[VIEW];TRA++){
//      if(ttrack[VIEW][TRA]){
//	std::cout << "XX " << XX[VIEW][TRA][ntracell[VIEW][TRA]] << " " << XX[VIEW][TRA][0] << " " << par[VIEW][TRA][0] << " " << par[VIEW][TRA][1] << "\n";
//	tline(mod,XX[VIEW][TRA][ntracell[VIEW][TRA]],XX[VIEW][TRA][0],par[VIEW][TRA][0],par[VIEW][TRA][1]);
//      }
//    }
//    
	if(VIEW==0){   
		xloliline();
		xloli_reconhit();
	}
	else if(VIEW==1){   
		yloliline();
		yloli_reconhit();
	}
	
	
    
  }
  c1->cd(0);
  c1->Update();
  
  if(batch){
    //c1->Print(Form("../data/image/EvtDisp_d160720c/EvtDisp_%06d.png",ievt)); //for developing by hosomi
  }
  else{
    printf("Type \'s\' to save the event display.\n");
    printf("Type \'q\' to quit.\n");
    
    char ans[8];
    fgets(ans,8,stdin);
    
    if( *ans == 's') c1->Print("~/Evtdisp.pdf");
    if( *ans == 'q') exit(0);
    
  }
  //c1->Print(Form("/export/scraid3/data/taichiro/cosmic_na_2016_4_26/picture_2016_5_4/test%d_hitbundle03_%d.jpg",ievt,alltrack.size()));
  //int nxtrack=0;
  //int nytrack=0;
  //for(int i=0;i<alltrack.size();i++){
  //      if(alltrack[i].view==0){nxtrack++;}
  //      if(alltrack[i].view==1){nytrack++;}
  //}
  //c1->Print(Form("/export/scraid3/data/taichiro/cosmic_na_2016_4_26/picture_2016_5_7/test_%d_%d_%d.jpg",ievt,nxtrack,nytrack));

  c1->Clear();

};

void EvtDisp_Analoli(){
  if(firstdraw){
    c1 = new TCanvas("c1","c1",(int)(2050*0.8/2*0.8),(int)(600*1.03*2*0.8));
    firstdraw=false;
  }
  //*****Draw module*****
  TH1F *h=new TH1F("","Side view",2050*10,0,1945);
  h->SetMinimum(0);
  h->SetMaximum(1200);
  h->GetXaxis()->SetLabelSize(0);
  h->GetYaxis()->SetLabelSize(0);
  h->SetStats(0);
  h->SetNdivisions(0);
  TH1F *v=new TH1F("","Top view",2050*10,0,1945);
  v->SetMinimum(0);
  v->SetMaximum(1200);
  v->GetXaxis()->SetLabelSize(0);
  v->GetYaxis()->SetLabelSize(0);
  v->SetStats(0);
  v->SetNdivisions(0);
  
  
  c1->Divide(1,2);
  c1->cd(1);
  h->Draw();
  drawlolix();
  c1->cd(2);
  v->Draw();
  drawloliy();
  
  //*****Start drawing*****
  for(int VIEW=0;VIEW<2;VIEW++){	
    c1->cd(VIEW+1);
  
    //*****Draw Vertex and Muon*****
    if(VIEW==0)xloliparticle();
    if(VIEW==1)yloliparticle();
   

    //*****Draw hits*****
    if(VIEW==0)xlolihit();
    else if(VIEW==1)ylolihit();


//    //*****Draw cells*****
//    
//    for(PLN=0;PLN<plnmax(mod,VIEW,PLN,1)-1;PLN++){
//    //for(PLN=0;PLN<plnmax(mod,VIEW,PLN,0)-1;PLN++){
//      for(DIST=0;DIST<6;DIST++){
//	//if(PLN==plnmax(mod)-2&&DIST==1)continue;
//	if(PLN==plnmax(mod,VIEW,PLN,1)-DIST)continue;
//	//if(PLN==plnmax(mod,VIEW,PLN,0)-2&&DIST==1)continue;
//	for(CELL=0;CELL<ncell[VIEW][PLN][DIST];CELL++){
//	  cline(mod,VIEW,PLN,clcenter[VIEW][PLN][cellu[VIEW][PLN][CELL][DIST]],clcenter[VIEW][PLN+1+DIST][celld[VIEW][PLN][CELL][DIST]],value[VIEW][PLN][CELL][DIST],DIST,1);
//	}
//      }
//    }
    

//    //*****Draw tracks*****
//    //std::cout << ntrack[VIEW] << "\n";
//    for(TRA=0;TRA<ntrack[VIEW];TRA++){
//      if(ttrack[VIEW][TRA]){
//	std::cout << "XX " << XX[VIEW][TRA][ntracell[VIEW][TRA]] << " " << XX[VIEW][TRA][0] << " " << par[VIEW][TRA][0] << " " << par[VIEW][TRA][1] << "\n";
//	tline(mod,XX[VIEW][TRA][ntracell[VIEW][TRA]],XX[VIEW][TRA][0],par[VIEW][TRA][0],par[VIEW][TRA][1]);
//      }
//    }
    
//	if(VIEW==0){   
//		xloliline_ana();
//	}
//	else if(VIEW==1){   
//		yloliline_ana();
//	}
	
	
//    //*****Draw 3D tracks*****

      for(int i=0;i<(int)pmtrack.size();i++){
      for(int j=0;j<(int)pmtrack[i].trk.size();j++){
	int mod  = 15;
	int view = VIEW;

	if(view==0){
		int ixpln = pmtrack[i].trk[j].startxpln;
		int fxpln = pmtrack[i].trk[j].endxpln;
		if(pmtrack[i].trk[j].ing_trk){
			fxpln=pmtrack[i].trk[j].ing_endpln;
			mod=3;
		}
		int ix    = (int)pmtrack[i].trk[j].x;
		double slopex = tan(3.14/180.*pmtrack[i].trk[j].thetax);
		tline_ana(15,zposi(15,view,ixpln,0,0),zposi(mod,view,fxpln,0,0),ix,slopex);
		std::cout << j << " " << mod << " " << view << " " << ixpln << " " << fxpln << " " << ix  << " " << slopex << " " << pmtrack[i].trk[j].thetax << "\n";
	}
	else if(view==1){
		int iypln = pmtrack[i].trk[j].startypln;
		int fypln = pmtrack[i].trk[j].endypln;
		if(pmtrack[i].trk[j].ing_trk){
			fypln=pmtrack[i].trk[j].ing_endpln;
			mod=3;
		}
		int iy    = (int)pmtrack[i].trk[j].y;
		double slopey = tan(3.14/180.*pmtrack[i].trk[j].thetay);
		tline_ana(15,zposi(15,view,iypln,0,0),zposi(mod,view,fypln,0,0),iy,slopey);
		std::cout << j << " " << mod << " " << view << " " << iypln << " " << fypln << " " << iy  << " " << slopey << " " << pmtrack[i].trk[j].thetay << "\n";
	}
      }
      }


    
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
