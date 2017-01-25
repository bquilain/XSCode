#ifndef Reconstruction_h
#define Reconstruction_h
#include<iostream>
#include<sstream>
#include<fstream>
using namespace std;
#include <iomanip>
#include <sys/stat.h>
#include <cmath>
#include "Hit.h"

class Reconstruction{
 public:

  vector <Hit3D> ApplyPEError(vector <Hit3D> Vec, double angle);
  vector <Hit3D> SearchIngridHit(vector <Hit3D> Vec, vector <Hit3D> VecAll, double thetaX, double thetaY, double TrackSample);
  //vector <HitTemp> EraseDoubleHits(IngridBasicReconSummary * recon, int itrk, vector <HitTemp> HitV);
  //vector <Hit3D> Hit2DMatching(IngridEventSummary* evt, IngridBasicReconSummary * recon,vector <HitTemp> HitV,vector <Hit3D> Vec, bool MC);
  vector <double> GetTrackAngle(vector <Hit3D> Vec);
  //vector <Hit3D> SeveralHitsPlane(IngridBasicReconSummary * recon,vector <HitTemp> HitV,vector <Hit3D> Vec);
  //vector <Hit3D> Hit2DMatchingCluster(IngridEventSummary* evt, IngridBasicReconSummary * recon,vector <HitTemp> HitV,vector <Hit3D> Vec);
  //vector <HitTemp> EraseDoubleHitsAllTracks(IngridBasicReconSummary * recon, vector <HitTemp> HitV);
  //vector <Hit3D> Hit2DMatchingAllTracks(IngridBasicReconSummary * recon);
  vector <HitTemp> EraseDoubleHitsPM(PMAnaSummary * recon, int itrk, vector <HitTemp> HitV);
  vector <Hit3D> Hit2DMatchingPM(IngridEventSummary* evt, PMAnaSummary * recon,vector <HitTemp> HitV,vector <Hit3D> Vec, bool MC);
  vector <double> GetTrackAnglePM(vector <Hit3D> Vec, double anglex, double angley, int TrackSample);
  vector <Hit3D> SeveralHitsPlanePM(PMAnaSummary * recon,vector <HitTemp> HitV,vector <Hit3D> Vec);
  vector <Hit3D> Hit2DMatchingClusterPM(IngridEventSummary* evt, PMAnaSummary * recon);
  vector <Hit3D> ClusterPM(IngridEventSummary* evt, PMAnaSummary * recon, int nTracks);
  vector <HitTemp> EraseDoubleHitsAllTracksPM(PMAnaSummary * recon, vector <HitTemp> HitV);
  vector <Hit3D> Hit2DMatchingAllTracksPM(PMAnaSummary * recon, bool MC);
  vector <double> GetAllTracksLength(vector <Hit3D> Vec);
  double GetTrackLength(vector <Hit3D> Vec);
  double GetTrackEnergyPerDistance(vector <Hit3D> Vec, double dx); 
  vector <double> GetVertex(vector <Hit3D> Vec);
  vector <int> TrackComposition(vector <Hit3D> Vec, double Vx, double Vy, double Vz);
  vector <double> GetKinematic(double ang1, double thetax1, double thetay1,double ang2, double thetax2, double thetay2);
  double GetBeamAngle(double ang1, double thetax1, double thetay1);
  vector <Hit3D> CountSharedHits(vector <Hit3D> Vec, vector< vector<Hit3D> > VecDouble, int Trk);
  vector <Hit3D> IsInTrk(vector <Hit3D> VecCluster, vector<Hit3D> VecAllTracks);
  vector <double> TrackRelativeLength(double posx, double posy, double posz, double tx, double ty, double TrkLength);
  vector <double> TheoreticalTrack(int startplnx, double startchx, double thetax, int startplny, double startchy, double thetay);
  bool HasGeomTrack(int mod, int startplnx, int startchx, double thetax, int startplny, int startchy, double thetay);
  vector <double> IngridTrack(int mod, int startplnx, int startchx, double thetax, int startplny, int startchy, double thetay);
  vector <double> TrackPenetration(int Mod, int pln_iniX, double ch_iniX, double thetax,int pln_iniY, double ch_iniY, double thetay, int pln_finX, double ch_finX, int pln_finY, double ch_finY, double dx_Ing);
  vector <double> TrackPenetrationPM(int pln_iniX, double ch_iniX, double thetax, int pln_iniY, double ch_iniY, double thetay,int pln_finX, double ch_finX, int pln_finY, double ch_finY, int IngMod, int pln_ini_Ing, int pln_fin_Ing, double dx_Ing, bool Geom);
  vector <double> ConvertTruePM(float ipos[4], float fpos[4]);
  int SelectTrackSample(bool pm_stop, bool Geom, bool has_ingrid, bool ingrid_stop, int ing_last_pln);
  double GetFSI(IngridEventSummary * evt);
  bool InPMFV(IngridSimVertexSummary * simver);
  vector <double> GetTrueMuonInformation(IngridEventSummary * evt);
  vector <double> GetTruePionInformation(IngridEventSummary * evt);//ML
  bool IsFV(int mod, double posx, double posy, double posz);
  int GetTrackParticle(IngridEventSummary * evt, PMAnaSummary * recon, int itrk, double TrkLength);
  bool IsINGRID(int mod, int pln, int ch);
  vector <double> GetMatchingPMINGRID(vector <Hit3D> Vec);
  double GetINGRIDTrackWidth(vector <Hit3D> Vec);
  vector <double> GetLastINGRIDChannel(vector <Hit3D> Vec, double TrackSample);
  void GetSelectionPM(bool * SelectionFV, bool * SelectionOV, PMAnaSummary * recon, bool MC);    

};
#endif
