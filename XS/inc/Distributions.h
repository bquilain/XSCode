#ifndef Distributions_h
#define Distributions_h
#include<iostream>
#include<sstream>
#include<fstream>
using namespace std;
#include <iomanip>
#include <sys/stat.h>
#include <cmath>
#include "Hit.h"
#include "setup.h"

class Distributions{
 public:

  Double_t Cumulative_PMIng(Double_t *pe, Double_t *par); Double_t Cumulative_PMSci(Double_t *pe, Double_t *par); Double_t Cumulative_Ing(Double_t *pe, Double_t *par);
  Double_t fPDF_PMIng(Double_t *pe, Double_t *par); Double_t fPDF_PMSci(Double_t *pe, Double_t *par); Double_t fPDF_Ing(Double_t *pe, Double_t *par);
};
#endif
