#ifndef _SETUP_H
#define _SETUP_H
#include<iostream>
#include<sstream>


#define NumMod     17
//#define NumMod     16
#define NumTPL     11
#define NumVETO    4
//#define NumTFB     11
#define NumTFB     24
#define UseNumTFB  13
#define NumCh      64
#define NumTPLCh   48 
#define NumVETOCh  22
#define UseNumCh   48
#define NumCyc     23
#define UseNumTFB  11
#define VetoNumTFB 2
#define UseNumCh   48
//#define LayerNumCh 24
#define LayerNumCh 32
#define NumLayer   2
#define GATE       800 //nsec


//for Proton Module
#define NumTPL_PM     18
#define NumVETO_PM    4
#define NumTFB_PM     22
#define UseNumTFB_PM  21
#define NumCh_PM      64
#define NumTPLCh_PM   64 
#define NumVETOCh_PM  17
#define UseNumCh_PM   64
#define UseNumTFB_PM  18
#define VetoNumTFB_PM 3
#define UseNumCh_PM   64
#define LayerNumCh_PM 32
#define NumLayer_PM   2

/*
#define NumTFB     24
#define UseNumTFB  13
#define NumCh      64
#define NumTPLCh   64
#define NumVETOCh  22
#define UseNumCh   64
#define NumCyc     23
#define UseNumTFB  18
#define VetoNumTFB 3
#define UseNumCh   48
//#define LayerNumCh 24                                                                                                                                                     
#define LayerNumCh 32
#define NumLayer   2
#define GATE       800 
*/

//definition of the name of folder
/*static const char *data_file_folder  = "/home/ingrid/data/daqdata";
static const char *calib_file_folder = "/home/daq/data/calib_file";
static const char *dailycheck_folder = "/home/daq/data/dailycheck";
static const char *gain_and_pedestal_folder = "/home/daq/data/gain_and_pedestal_file";
*/
static const long Tdc_max=15000000;
static const long Tdc_min=10;

static const double CUT=2.5;
static const int cosmic_event_cycle_1=18;
static const int cosmic_event_cycle_2=19;

//static const char *data_file_folder  = "/home/bquilain/Ingrid_Process/daqdata";                    
static const char *data_file_folder  = "/export/scraid2/data/ingrid/daqdata";
                         
static const char *calib_file_folder = "/home/bquilain/Ingrid_Process/calib_file";                     
                        
static const char *dailycheck_folder = "/home/bquilain/Ingrid_Process/dailycheck";
                                             
static const char *gain_and_pedestal_folder = "/home/bquilain/Ingrid_Process/gain_and_pedestal_file";
bool isPM(int rmm, int tfb){
  if( rmm==4 && tfb>=26)
    return true;
  else
    return false;
}

bool isPM(int mod){
  if( mod==16 )
    return true;
  else
    return false;
}



#endif
