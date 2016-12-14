#ifndef Corrections_h
#define Corrections_h
#include<iostream>
#include<sstream>
#include<fstream>
using namespace std;
#include <iomanip>
#include <vector>
#include <sys/stat.h>
#include <cmath>
#include "Hit.h"

class Corrections{
 public:
  //Method
  double GetMCCorrections(double NbEvent,int Module);
  vector <Hit3D> GetFiberAttenuation(vector <Hit3D> Vec);
  vector <Hit3D> GetDXCorrection(vector <Hit3D> Vec, double dx);
  vector <Hit3D> RemoveFiberAttenuation(vector <Hit3D> Vec);
  void GetHitPlaneCorrection(vector <Hit3D> Vec);
};

#endif
