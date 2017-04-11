#ifndef Hit_h
#define Hit_h
#include<iostream>
#include<sstream>
#include<fstream>
using namespace std;
#include <iomanip>
#include <sys/stat.h>
#include <cmath>


class Hit2D{
 public:
  int mod;
  int view;
  double pe;
  double xy;
  double z;

  void clear(){
    mod=1e-5;
    view=1e-5;
    pe=0;
    xy=0;
    z=0;
  }
};


class Hit3D{
 public:
  int mod;
  int pln;
  int view;//la view du hit d'origin                                                                                                                                       
  int ch;
  double pe;
  double pecorr;
  double truepe; //if <0, dn hit                                                                                                                                           
  int pdg;
  double x;
  double y;
  double z;
  double time;
  int inttype;
  int used;
  bool InTrk;
  int trackid;
  double dist_plastic;
  double dist_iron;
  vector <int> RecTrk;

  void clear(){
    pln=1e-5;
    mod=1e-5;
    view=1e-5;
    ch=0;
    pe=0;
    pecorr=0;
    truepe=0;
    pdg=-1;
    x=0;
    y=0;
    z=0;
    time=1e-5;
    inttype=0;
    used=0;
    InTrk=false;
    trackid=0;
    RecTrk.clear();
    dist_plastic=0;
    dist_iron=0;
  }

  friend bool operator < (const Hit3D left, const Hit3D right)
  {
    return (left.z < right.z ? true : false);
  }
  friend bool operator == (Hit3D A, Hit3D B)
  {
    return (((A.mod == B.mod)&&(A.pln == B.pln)&&(A.ch == B.ch)&&(A.view == B.view)&&(A.time == B.time)) ? true : false);
  }
};


class HitTemp{
 public:
  //int x;
  //double y;
  //double z;
  int view;
  int ch;
  int pln;
  int trk;
  int hit;
  int mod;

  void clear(){
    ch=0;
    pln=0;
    view=0;
    trk=0;
    hit=0;
    mod=-1;
  }
  friend bool operator == (HitTemp A, HitTemp B)
  {
    return (((A.view == B.view)&&(A.pln == B.pln)&&(A.ch == B.ch)&&(A.mod == B.mod)) ? true : false);
  }
};


class RecTrack{
 public:
  double Mom;
  double IronDist;
  double PlasticDist;
  double MuCL;
  int Sample;
  double Angle;
  double TrueAngle;
  int TruePart;

  void clear(){
  Mom=0;
  IronDist=0;
  PlasticDist=0;
  MuCL=0;
  Sample=-1;
  Angle=0;
  TrueAngle=0;
  TruePart=0;
  }
  friend bool operator < (const RecTrack left, const RecTrack right)
  {
    return (left.MuCL < right.MuCL ? true : false);
  }
};

#endif
