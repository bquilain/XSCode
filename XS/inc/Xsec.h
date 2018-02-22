#ifndef Xsec_h
#define Xsec_h
#include<iostream>
#include<sstream>
#include<fstream>
using namespace std;
#include <iomanip>
#include <sys/stat.h>
#include <cmath>
#include "Hit.h"
#include "setup.h"

class Xsec{
 public:
  TRandom3 * r;
  Xsec(bool=true);
  ~Xsec();
  
  void Initialize();
  void SetDetector(bool);
  bool GetDetector();

  int DetermineFSI(int IsSand,int IsAnti,int IsNuE,int IsBkgH,int IsBkgV, int IsSciBkg,IngridEventSummary * evt);
  void DetermineNuType(bool& IsSand,bool& IsAnti,bool& IsNuE,bool& IsBkgH,bool& IsBkgV, bool & IsSciBkg,int nutype, int intmode, double* pos);
    
  void LoadMuCLDistributions();
  
  double GetMuCL(TF1 * CL_Ing,TF1 * CL_PMIng,TF1 * CL_PMSci,vector <Hit3D> Vec, double dx, int TrackSample, bool SystematicPE, int RandomPE);
  double GetMuCL_Plan(TF1 * CL_Ing,TF1 * CL_PMIng,TF1 * CL_PMSci,vector <Hit3D> Vec, double dx, int TrackSample, bool SystematicPE, int RandomPE);
  
  void LoadInputFiles_OnlyUnfoldedData(char * InputName, double **DataUnfoldedEvents, double ** DataTrueEvents);
    
  void LoadInputFiles_OnlySelectedData(char * fDataName,double ** DataReconstructedEvents);

  void LoadInputFiles_OnlySelectedDataSB(char * fDataName,char * fDataNameSB,double ** DataReconstructedEvents);
    
  void LoadInputFiles(char * fDataName,char * fMCName, double **** MCReconstructedEvents_TrueSignal, double ** DataReconstructedEvents,double ** MCReconstructedEvents, double ** MCReconstructedBkgEvents, double ** Efficiency, double *NumberOfPOT,double **** MCReconstructedEvents_TrueSignal=NULL,double ** DataEfficiency=NULL,bool FakeData=false);

  void LoadInputFilesSB(char * fDataName,char * fMCName, char * fDataNameSB,char * fMCNameSB, double **** MCReconstructedEvents_TrueSignal, double ** DataReconstructedEvents,double ** MCReconstructedEvents, double ** MCReconstructedBkgEvents, double ** Efficiency, double *NumberOfPOT,double **** MCReconstructedEvents_TrueSignal=NULL,double ** DataEfficiency=NULL,bool FakeData=false);
  //void LoadInputFilesSB(char * fDataName,char * fMCName, char * fDataNameSB,char * fMCNameSB, double **** MCReconstructedEvents_TrueSignal, double ** DataReconstructedEvents,double ** MCReconstructedEvents, double ** MCReconstructedBkgEvents, double ** Efficiency, double *NumberOfPOT);

  void LoadNeutrinoFlux(TH1D * NeutrinoFlux);
  
  void BuildLikelihood(double **** vLikelihood, double ** vPriorMC, double **** MCReconstructedEvents_TrueSignal);

  void ProjectOnTruePhaseSpace_OnlySignal(double ** TrueEventsDistribution_SignalOnly, double **** ReconstructedEvents_TrueSignal);
    
  void BuildUnfolding(double **** vUnfolding, double **** vLikelihood, double ** vPrior);

  void ApplyUnfoldingBkgSubstraction(double ** vPosteriorEvents, double **** vUnfolding, double ** DataReconstructedEvents, double ** MCReconstructedBkgEvents);

  void ApplyUnfolding(double ** vPosteriorEvents, double ** vPosteriorEvents_SignalOnly, double **** vUnfolding, double ** DataReconstructedEvents);

  void SetPrior(double ** vPriorNormalised, double ** vPrior, double ** vInitialPriorMC, double ** vInitialPrior,double ** vPosterior, bool PriorMC, double IterationStep);

  void GenerateStatisticalFluctuations(double ** DataReconstructedEvents);

  void GenerateDetectorMCFluctuations(double **** MCReconstructedEvents_TrueSignal,double **** RelativeSigma);
  void GenerateXSModelMCFluctuations(double ** Efficiency, double **** MCReconstructedEvents_TrueSignal, double **** RelativeSigma);

  void GenerateDataFluctuations(double ** DataReconstructedEvents, double **** DataReconstructedEvents_TrueSignal,double **** RelativeSigma);
    
 private:
  bool _isPM;
};
#endif
