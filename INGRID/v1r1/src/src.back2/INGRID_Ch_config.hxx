#ifndef _INGRID_Ch_config_H
#define _INGRID_Ch_config_H

#include<iostream>
#include<sstream>
#include<fstream>
using namespace std;

#define nRMM 5

class INGRID_Ch_config{
private:
  int ReadoutId[300][2];
  //ReadoutId[][0]:RMM 0~4
  //ReadoutId[][1]:TFB 0~47

  int PhysId[300][2];
  bool tpl_flag[300];
  //PhysId[][0]:Module 0~14
  //PhysId[][1]:TPL 0~10 VETO 0~3
  //tpl[]=true:tracking plane
  //tpl[]=false:VETO plane

  bool horizontal;
  bool vertical;
  bool offdiag;
  bool full_veto;

  int index;
public:
  INGRID_Ch_config(){
    //make a list of RMM and TFB
    index=0;
    for(int numrmm=0; numrmm<nRMM; numrmm++){
      for(int numtfb=0;numtfb<48;numtfb++){
	if((numrmm==1||numrmm==3)&&numtfb==43)
	  break;
	else if( numrmm==4 && numtfb==26 )
	  break;
	ReadoutId[index][0]=numrmm;
	ReadoutId[index][1]=numtfb;
	index++;
      }//numtfb
    }//numrmm
    
    //make a list of nummod and numplane
    index=0;
    for(int nummod=0;nummod<16;nummod++){
      for(int numplane=0;numplane<13;numplane++){
	if(nummod<7){
	  PhysId[index][0]=6-nummod;
	}
	else if(nummod>=7&&nummod<16){
	  PhysId[index][0]=nummod;
	}


	PhysId[index][1]=numplane;
	if(numplane==11||numplane==12){tpl_flag[index]=false;}
	else{tpl_flag[index]=true;}
	index++;
      }//numplane
    }//nummod
  };
  ~INGRID_Ch_config(){};


  bool channel_configuration(int *rmm,int *tfb,int *trip,int *trip_ch,int *mod,int *plane,int *ch, bool *tpl, bool *veto);
  //rmm=0~4, tfb=0~47, trip=0~3, trip_ch=0~15
  //mod=0~13, plane=0~10(tpl), 0~4(veto), ch=0~47,0~21
  bool channel_configuration(int *rmm,int *tfb,int *trip,int *trip_ch,int *mod,int *view,int *plane,int *ch);
  //new channel_configuration for new data structure
  //rmm=0~4, tfb=0~47, trip=0~3, trip_ch=0~15
  //mod=0~13, plane=0~10(tpl), 0~4(veto), ch=0~47,0~21
  bool active_veto(int *mod,int *plane);
  bool active_pln(int mod, int view, int plane);
  bool get_global_pln(int *mod,int *plane, bool *tpl, bool *veto, int *global_pln);
  bool get_global_ch(int *mod, int *plane,int *ch, bool *tpl, bool *veto, int *global_ch);
  bool get_global_ch(int *mod, int *view, int *plane,int *ch, int *global_ch);
  //for new data structure
  bool get_global_ch(int mod, int view, int plane,int ch, int *global_ch);

};
#endif

