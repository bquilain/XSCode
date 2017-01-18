#ifndef _INGRID_Dimension_H
#define _INGRID_Dimension_H

#include<iostream>
#include<sstream>
#include<fstream>
using namespace std;

#include"IngridConstants.h"

#define CHMAX 80
#define MOD_WAGASCI 15

#ifndef __INGRID_CONSTANTS__
const static double IronThick       =   6.5   ;     //cm
const static double ScintiWidth     =   5.0   ;     //cm
const static double ScintiThick     =   1.0   ;     //cm
const static double PlnThick        =   4.2   ;     //cm
const static double VetoOffsetZX    =   -2.0  ;     //cm
const static double VetoStartZ      =   +1.0  ;     //cm
const static double VetoOffsetZY    =   -1.0  ;     //cm
const static double VetoOffsetRight = -10.575-2.5 ; //cm
const static double VetoOffsetLeft  = 130.9  -2.5 ; //cm
const static double VetoOffsetBottom=  -5.9   ;     //cm
const static double VetoOffsetUp    = 128.4   ;     //cm
#endif

const static double ScibarWidth     =   2.5   ;     //cm
const static double PlnThick_PM     =   4.6   ;     //cm
const static double PlnThick_front_PM=   5.0   ;     //cm
const static double PlnDist_PM      =   2.3   ;     //cm
const static double VetoOffsetZX_PM    =   -1.55  ;     //cm
const static double VetoOffsetZY_PM    =   -1.55  ;     //cm
const static double VetoOffsetRight_PM = -6 ; //cm
const static double VetoOffsetLeft_PM  = 125 ; //cm
const static double VetoOffsetBottom_PM=  -6   ;     //cm
const static double VetoOffsetUp_PM    = 125   ;     //cm
const static double VetoStartZ_PM      =   -0.55  ;     //cm

//general
const int    loli_plnnum	  = 8;
const int    loli_chnum		  = 40;
const int    loli_gridchnum	  = 20;
const double loli_watersurface_z  = -23.3;//cm upstream water surface  now half of waterZ length (46.6cm)
const double loli_firstdistance_z = 0.5;  //cm distance between water target and scintillator
const double loli_firstoffset_z   = loli_watersurface_z + loli_firstdistance_z;//cm distance between water target and scintillator
const double loli_scinti_thick 	  = 0.3;  //cm thickness of scinti
const double loli_scinti_width 	  = 2.5;  //cm width of scinti
const double loli_offsetxy        = - loli_scinti_width * loli_chnum/2.;
const double loli_gap 		  = 5.7;  //cm gap between hlayer and hlayer
const double loli_offset_hv 	  = 2.85; //cm gap between hlayer and vlayer
//grid layer
const double loli_cutgap	  = 5.0;   //cm gap of grid cut
const double loli_cutwidth	  = 0.35;  //cm width of grid cut
const double loli_cutthick	  = 1.3;   //cm thickness of grid cut
const double loli_offsetxy_grid	  = - loli_scinti_width * loli_chnum/2. +2.335;//cm distance between first grid and edge of scintillator 

class INGRID_Dimension{
private:

  double VETOOffsetZ;
public:
  INGRID_Dimension(){

  };
  ~INGRID_Dimension(){};
  bool get_pos(int mod, int pln, int ch, bool tpl, bool veto, double *posxy, double *posz);

  //########## For New Data Structure #############
  //###############################################
  bool get_posXY(int mod, int view, int pln, int ch,
		 double *posxy, double *posz);
  bool get_posVeto(int mod, int view, int pln, int ch, 
		   double *posxy, double *posz);//for new data structure (not complete)
  //###############################################


  bool get_expch(int mod, int pln, int *ch, bool tpl, bool veto, double a, double b);
  bool get_expchXY(int mod, int view, int pln, int *ch, double a, double b);

  //added for prototype of WAGASCI
  bool get_pos_loli(int mod, int view, int pln, int ch, int grid,
			double *posx, double *posy, double *posz);
  bool get_pos_loli(int mod, int view, int pln, int ch,
			double *posx, double *posy, double *posz);
  bool get_grid_loli(int mod, int view, int pln, int ch,
			int *grid, int *gridch);
  bool get_pos_loli_xy(int mod, int view, int pln, int ch,
			 double* posxy, double* posz);
  bool get_pos_loli_xy(int mod, int view, int pln, int gridch, int grid,
			 double* posxy, double* posz);

  //for Lolirecon
  bool get_reconplnch_loli(int mod, int view, int pln, int ch, int axis,
			 int* reconpln, int* reconch);

  bool get_plnch_fromrecon_loli(int mod, int view, int reconpln, int reconch, int axis,
			 int* pln, int* gridch, int* grid);

  int get_chmax_loli(int mod, int view, int pln, int axis);

  int get_plnmax_loli(int mod, int view, int pln, int axis);

  bool get_posi_lolirecon(int mod, int view, int reconpln, int reconch, int axis,
			 double* posxy, double* posz);

  float get_sciwidth(int mod, int view, int reconpln, int reconch, int axis);

  float get_scithick(int mod, int view, int reconpln, int reconch, int axis);

  void get_loli_gridcell_id(int mod, int view, int pln, int ch, double posx, double posy, double posz,
				int* gridcell_id_x1, int* gridcell_id_x2, int* gridcell_id_y1, int* gridcell_id_y2);

  void get_loli_scintiid_from_cellid(int mod, int view, int pln, int gridcell_id, int* cross_n, int* ch);
};
#endif

