#ifndef _INGRID_Geo_C
#define _INGRID_Geo_C
#include<iostream>
#include<sstream>
#include"INGRID_Geometry.hxx"
using namespace std;

bool INGRID_Geometry::fGet_Geometry(int *mod, int *plane, int *ch, bool *tpl, bool *veto, double *posx, double *posy, double *posz){
  if(tpl){
    *posz = 1.0 * (*plane) * (ThickTPL+ThickIron);
    if(*mod!=3){
      if(*ch<24){
	*posy = PosYModule[*mod] + 1.0 * WidthScinti * (*ch + 0.5);
      }//Y layer
      if(*ch>=24){
	*posx = PosXModule[*mod] + 1.0 * WidthScinti * (23 - *ch + 0.5);
      }//X layer
    }//mod!=3
    if(*mod!=3){
      if(*ch<24){
	*posx = PosXModule[*mod] + 1.0 * WidthScinti * (23 - *ch + 0.5);
      }//X layer
      if(*ch>=24){
	*posy = PosYModule[*mod] + 1.0 * WidthScinti * (*ch + 0.5);
      }//Y layer
    }//mod!=3
    return true;
  }
  if(veto){
    return true;
  }

}
#endif
