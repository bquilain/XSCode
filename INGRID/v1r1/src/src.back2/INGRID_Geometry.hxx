#ifndef _INGRID_Geo_H
#define _INGRID_Geo_H
#include<iostream>
#include<sstream>
using namespace std;

#define XLayer 0
#define YLayer 1
#define RVETO 0
#define LVETO 1
#define BVETO 2
#define UVETO 3
const static double PosXModule[14]={0,0,0,0,0,0,0,0,0,0,0,0,0,0};
const static double PosYModule[14]={0,0,0,0,0,0,0,0,0,0,0,0,0,0};
static double WidthScinti=5.0;//cm
const static double LengthScinti=120.0;//cm
const static double ThickScinti=1.0;//cm
const static double ThickTPL=3.2;//cm
const static double ThickIron=6.4;//cm

class INGRID_Geometry{
public:
  bool fGet_Geometry(int *mod, int *plane, int *ch, bool *tpl, bool *veto, double *posx, double *posy, double *posz);

};
#endif
