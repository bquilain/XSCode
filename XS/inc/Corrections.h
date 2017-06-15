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
  Corrections(bool=true);
  ~Corrections();
  void SetDetector(bool);
  bool GetDetector();

  double GetMCCorrections(double NbEvent,int Module);
  vector <Hit3D> GetFiberAttenuation(vector <Hit3D> Vec);
  vector <Hit3D> GetDXCorrection(vector <Hit3D> Vec, double dx);
  void GetDXCorrectionWM(vector <Hit3D> & Vec,double angle3D, double thetax, double thetay);
  vector <Hit3D> RemoveFiberAttenuation(vector <Hit3D> Vec);
  void GetHitPlaneCorrection(vector <Hit3D> Vec);
  double dzWM(double angle, double theta, bool grid);

 private:
  bool _isPM;
  INGRID_Dimension* IngDimCor;
};

#endif
