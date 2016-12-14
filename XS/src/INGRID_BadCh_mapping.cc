#ifndef _INGRID_BadCh_mapping_C
#define _INGRID_BadCh_mapping_C

#include<iostream>
#include<sstream>
#include<fstream>
#include"INGRID_BadCh_mapping.h"
using namespace std;
static const int TPL = 0;
static const int VETO = 1;

bool INGRID_BadCh_mapping::badchannel(int mod, int view, int pln,int ch){
  int number_of_bad=17;
  int bad[17][4]={
    // mod, view,  pln,   ch
    {    0,    0,   13,   12}, //cable damage channel
    {    1,    1,    1,   14}, //cable damage channel
    {    5,    1,    2,    0}, //cable damage channel
    {    5,    1,    4,    2}, //high gain channel
    {    5,    1,    4,   18}, //cable damage channel
    {    5,    1,    7,    9}, //cable damage channel
    {    7,    0,    5,    0}, //pedestal channel
    {    9,    1,    0,   12}, //cable damage chnnale
    {    0,    1,    8,   14}, //pedestal shift channel
    {   11,    0,    5,   13}, //pedestal shift channel

    {    2,    0,   13,    2}, //added 2010/11/1
    {    5,    1,    4,   20}, //added 2010/11/1
    {    7,    0,    6,    9}, //added 2010/11/1
    {    5,    0,    4,   23}, //added 2010/11/1 HighGain
    {   15,    1,    8,   14}, //added 2010/11/1 cable damaged
    {   15,    1,   10,   22}, //added 2010/11/1 cable damaged??
    {    9,    0,    4,   11}  //added 2010/11/10 ??
  };
 
  for(int i=0;i<number_of_bad;i++){
    if(bad[i][0]==mod){
      if(bad[i][1]==view){
	if(bad[i][2]==pln){
	  if(bad[i][3]==ch){
	    return true;
	  }//channel
	}//planen
      }//TPL or VETO
    }//mod
  }//i
  return false;  

}

bool INGRID_BadCh_mapping::badchannel(int *mod,int *plane,int *ch,bool *tpl,bool *veto){
  int number_of_bad=8;
  int bad[8][4]={
      //{module, TPLorVETO, Plane, Channel}
    //### channel 0 ~ 23 -> view = 1, channel 24 ~ 48 -> view = 0
      0, VETO, 2, 12,
      1, TPL, 1, 14,
      5, TPL, 2, 0,
      5, TPL, 4, 18,
      5, TPL, 4,  2,
      7, TPL, 5, 24,
      9, TPL, 0, 12,
      5, TPL, 7, 9
  };
 
  for(int i=0;i<number_of_bad;i++){
    if(bad[i][0]==*mod){
      if((*tpl&&bad[i][1]==IDTPL)||(*veto&&bad[i][1]==IDVETO)){
	if(bad[i][2]==*plane){
	  if(bad[i][3]==*ch){
	    return true;
	  }//channel
	}//planen
      }//TPL or VETO
    }//mod
  }//i
  return false;  
}


#endif
