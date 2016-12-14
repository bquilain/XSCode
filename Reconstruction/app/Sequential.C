#include<iomanip>
#include<iostream>
#include <cstdio>
#include <sstream>
#include <fstream>
using namespace std;

int main(int argc, char **argv){

  int NFiles=500;
  char Name1[500], Name2[500];
  char Command1[500], Command2[500], Command3[500], Command4[500], Command11[500];
  char CommandR1[500], CommandR10[500], CommandR2[500], CommandR20[500],CommandR3[500], CommandR30[500], CommandR4[500], CommandR40[500];
  for(int i=1;i<=NFiles;i++){
    sprintf(Name1,"AnaBeam%d.sh",i);
    //sprintf(Command1,"./PMrecon_HitInfo -f ../../Ing_MC_Beam%d_NewGainNewMIP9_wNoise.root -o ../../Ing_MC_Beam%d_NewGainNewMIP9_wNoise_recon.root",i,i);
    //sprintf(Command1,"PM_Ana/app/IngAddNoiseMC -f Ingrid_Process/Ing_%d_MIP5.root -o Ingrid_Process/Ing_MC_RealBeam%d_NewGainNewMIP5_wNoise.root",i,i);
    //sprintf(Command1,"./PMrecon_HitInfo -r 14000 -s 0 -f ../../Ing_MC_Sand%d_NewGainNewMIP9_Form4_wNoise.root -o ../../Ing_MC_Sand%d_NewGainNewMIP9_Form4_wNoise_recon.root",i,i);
    //sprintf(Command11,"./IngAnaVertex_HitInfo -r 14000 -s 0 -f ../../Ing_MC_Sand%d_NewGainNewMIP9_Form4_wNoise_recon.root -o ../../Ing_MC_Sand%d_NewGainNewMIP9_Form4_wNoise_ana.root",i,i);
    //sprintf(Command1,"./PMrecon_HitInfo -r 14000 -s 0 -f ../../Ing_MC_Beam%d_MIP9_NuE_HighStat.root -o ../../Ing_MC_Beam%d_MIP9_NuE_HighStat_recon.root",i,i);
    //sprintf(Command11,"./IngAnaVertex_HitInfo -r 14000 -s 0 -f ../../Ing_MC_Beam%d_MIP9_NuE_HighStat_recon.root -o ../../Ing_MC_Beam%d_MIP9_NuE_HighStat_ana.root",i,i);
   
    //sprintf(Command1,"./PMrecon_HitInfo -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/Ing_MC_Beam%dNuE_BirksOKLargeSquareSci_wNoise.root -o /export/scraid2/data/bquilain/MCfiles/Ing_MC_Beam%dNuE_BirksOKLargeSquareSci_wNoise_recon.root",i,i);
    //sprintf(Command11,"./IngAnaVertex_HitInfo -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/Ing_MC_Beam%dNuE_BirksOKLargeSquareSci_wNoise_recon.root -o /export/scraid2/data/bquilain/MCfiles/Ing_MC_Beam%dNuE_BirksOKLargeSquareSci_wNoise_ana.root",i,i);





    //Here are the files I'm using finally (sand and beam)

    //sprintf(Command1,"./PMreconRev -r 14000 -s %d -f /export/scraid2/data/bquilain/MCfiles/Ing_MC_Sand%d_KillNeutrons5MeVNewMIP_wNoise.root -o Test2.root",i,i,i);
    //sprintf(Command1,"./PMreconRev -r 14000 -s %d -f /export/scraid2/data/bquilain/MCfiles/Ing_MC_Sand%d_KillNeutrons5MeVNewMIP_wNoise.root -o Test2.root",i,i,i);
    //sprintf(Command1,"./PMrecon_HitInfo -r 14000 -s %d -f /export/scraid2/data/bquilain/MCfiles/Ing_MC_Sand%d_KillNeutrons5MeVNewMIP_wNoise.root -o Test2.root",i,i,i);   
    /*
    sprintf(Command1,"/export/scraid3/data/bquilain/Programs/Programs/INGRID_analysis/app/PMreconRev_KS_woXTalk -r 14000 -s %d -f /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_Run0_wNoise.root -o /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_Run0_wNoise_recon.root",i,i,i);
    sprintf(Command11,"/export/scraid3/data/bquilain/Programs/Programs/INGRID_analysis/app/IngAnaVertex_HitInfo -r 14000 -s %d -f /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_Run0_wNoise_recon.root -o /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_Run0_wNoise_ana.root",i,i,i);

    sprintf(CommandR1,"/export/scraid3/data/bquilain/Programs/Programs/INGRID_analysis/app/PMreconRev_KS_woXTalk -r 14000 -s %d -f /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_Run1_wNoise.root -o /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_Run1_wNoise_recon.root",i,i,i);
    sprintf(CommandR10,"/export/scraid3/data/bquilain/Programs/Programs/INGRID_analysis/app/IngAnaVertex_HitInfo -r 14000 -s %d -f /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_Run1_wNoise_recon.root -o /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_Run1_wNoise_ana.root",i,i,i);

    sprintf(CommandR2,"/export/scraid3/data/bquilain/Programs/Programs/INGRID_analysis/app/PMreconRev_KS_woXTalk -r 14000 -s %d -f /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_Run2_wNoise.root -o /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_Run2_wNoise_recon.root",i,i,i);
    sprintf(CommandR20,"/export/scraid3/data/bquilain/Programs/Programs/INGRID_analysis/app/IngAnaVertex_HitInfo -r 14000 -s %d -f /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_Run2_wNoise_recon.root -o /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_Run2_wNoise_ana.root",i,i,i);

    sprintf(CommandR3,"/export/scraid3/data/bquilain/Programs/Programs/INGRID_analysis/app/PMreconRev_KS_woXTalk -r 14000 -s %d -f /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_Run3_wNoise.root -o /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_Run3_wNoise_recon.root",i,i,i);
    sprintf(CommandR30,"/export/scraid3/data/bquilain/Programs/Programs/INGRID_analysis/app/IngAnaVertex_HitInfo -r 14000 -s %d -f /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_Run3_wNoise_recon.root -o /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_Run3_wNoise_ana.root",i,i,i);

    sprintf(CommandR4,"/export/scraid3/data/bquilain/Programs/Programs/INGRID_analysis/app/PMreconRev_KS_woXTalk -r 14000 -s %d -f /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_Run4_wNoise.root -o /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_Run4_wNoise_recon.root",i,i,i);
    sprintf(CommandR40,"/export/scraid3/data/bquilain/Programs/Programs/INGRID_analysis/app/IngAnaVertex_HitInfo -r 14000 -s %d -f /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_Run4_wNoise_recon.root -o /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_Run4_wNoise_ana.root",i,i,i);
*/
    /*
    sprintf(Command1,"/export/scraid3/data/bquilain/Programs/Programs/INGRID_analysis/app/PMreconRev_KS_woXTalk -r 14000 -s %d -f /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46__HorizShiftMinus0.2cm_wNoise.root -o /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46__HorizShiftMinus0.2cm_wNoise_recon.root",i,i,i);
    sprintf(Command11,"/export/scraid3/data/bquilain/Programs/Programs/INGRID_analysis/app/IngAnaVertex_HitInfo -r 14000 -s %d -f /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46__HorizShiftMinus0.2cm_wNoise_recon.root -o /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46__HorizShiftMinus0.2cm_wNoise_ana.root",i,i,i);

    sprintf(CommandR1,"/export/scraid3/data/bquilain/Programs/Programs/INGRID_analysis/app/PMreconRev_KS_woXTalk -r 14000 -s %d -f /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46__HorizShiftMinus0.4cm_wNoise.root -o /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46__HorizShiftMinus0.4cm_wNoise_recon.root",i,i,i);
    sprintf(CommandR10,"/export/scraid3/data/bquilain/Programs/Programs/INGRID_analysis/app/IngAnaVertex_HitInfo -r 14000 -s %d -f /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46__HorizShiftMinus0.4cm_wNoise_recon.root -o /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46__HorizShiftMinus0.4cm_wNoise_ana.root",i,i,i);

    sprintf(CommandR2,"/export/scraid3/data/bquilain/Programs/Programs/INGRID_analysis/app/PMreconRev_KS_woXTalk -r 14000 -s %d -f /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46__HorizShiftMinus0.6cm_wNoise.root -o /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46__HorizShiftMinus0.6cm_wNoise_recon.root",i,i,i);
    sprintf(CommandR20,"/export/scraid3/data/bquilain/Programs/Programs/INGRID_analysis/app/IngAnaVertex_HitInfo -r 14000 -s %d -f /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46__HorizShiftMinus0.6cm_wNoise_recon.root -o /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46__HorizShiftMinus0.6cm_wNoise_ana.root",i,i,i);

    sprintf(CommandR3,"/export/scraid3/data/bquilain/Programs/Programs/INGRID_analysis/app/PMreconRev_KS_woXTalk -r 14000 -s %d -f /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46__HorizShiftMinus0.8cm_wNoise.root -o /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46__HorizShiftMinus0.8cm_wNoise_recon.root",i,i,i);
    sprintf(CommandR30,"/export/scraid3/data/bquilain/Programs/Programs/INGRID_analysis/app/IngAnaVertex_HitInfo -r 14000 -s %d -f /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46__HorizShiftMinus0.8cm_wNoise_recon.root -o /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46__HorizShiftMinus0.8cm_wNoise_ana.root",i,i,i);

    sprintf(CommandR4,"/export/scraid3/data/bquilain/Programs/Programs/INGRID_analysis/app/PMreconRev_KS_woXTalk -r 14000 -s %d -f /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_HShit1.0_wNoise.root -o /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_HShit1.0_wNoise_recon.root",i,i,i);
    sprintf(CommandR40,"/export/scraid3/data/bquilain/Programs/Programs/INGRID_analysis/app/IngAnaVertex_HitInfo -r 14000 -s %d -f /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_HShit1.0_wNoise_recon.root -o /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_HShit1.0_wNoise_ana.root",i,i,i);
    */
 sprintf(Command1,"/export/scraid3/data/bquilain/Programs/Programs/INGRID_analysis/app/PMreconRev_KS_woXTalk -r 14000 -s %d -f /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_HorizontalShift0.5cm_wNoise.root -o /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_HorizontalShift0.5cm_wNoise_recon.root",i,i,i);
    sprintf(Command11,"/export/scraid3/data/bquilain/Programs/Programs/INGRID_analysis/app/IngAnaVertex_HitInfo -r 14000 -s %d -f /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_HorizontalShift0.5cm_wNoise_recon.root -o /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_HorizontalShift0.5cm_wNoise_ana.root",i,i,i);

    sprintf(CommandR1,"/export/scraid3/data/bquilain/Programs/Programs/INGRID_analysis/app/PMreconRev_KS_woXTalk -r 14000 -s %d -f /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_HorizontalShift1.5cm_wNoise.root -o /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_HorizontalShift1.5cm_wNoise_recon.root",i,i,i);
    sprintf(CommandR10,"/export/scraid3/data/bquilain/Programs/Programs/INGRID_analysis/app/IngAnaVertex_HitInfo -r 14000 -s %d -f /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_HorizontalShift1.5cm_wNoise_recon.root -o /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_HorizontalShift1.5cm_wNoise_ana.root",i,i,i);

    sprintf(CommandR2,"/export/scraid3/data/bquilain/Programs/Programs/INGRID_analysis/app/PMreconRev_KS_woXTalk -r 14000 -s %d -f /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_HorizontalShift3.0cm_wNoise.root -o /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_HorizontalShift3.0cm_wNoise_recon.root",i,i,i);
    sprintf(CommandR20,"/export/scraid3/data/bquilain/Programs/Programs/INGRID_analysis/app/IngAnaVertex_HitInfo -r 14000 -s %d -f /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_HorizontalShift3.0cm_wNoise_recon.root -o /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_HorizontalShift3.0cm_wNoise_ana.root",i,i,i);

    sprintf(CommandR3,"/export/scraid3/data/bquilain/Programs/Programs/INGRID_analysis/app/PMreconRev_KS_woXTalk -r 14000 -s %d -f /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_HorizontalShift6.0cm_wNoise.root -o /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_HorizontalShift6.0cm_wNoise_recon.root",i,i,i);
    sprintf(CommandR30,"/export/scraid3/data/bquilain/Programs/Programs/INGRID_analysis/app/IngAnaVertex_HitInfo -r 14000 -s %d -f /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_HorizontalShift6.0cm_wNoise_recon.root -o /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_HorizontalShift6.0cm_wNoise_ana.root",i,i,i);

    sprintf(CommandR4,"/export/scraid3/data/bquilain/Programs/Programs/INGRID_analysis/app/PMreconRev_KS_woXTalk -r 14000 -s %d -f /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_HorizontalShift9.0cm_wNoise.root -o /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_HorizontalShift9.0cm_wNoise_recon.root",i,i,i);
    sprintf(CommandR40,"/export/scraid3/data/bquilain/Programs/Programs/INGRID_analysis/app/IngAnaVertex_HitInfo -r 14000 -s %d -f /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_HorizontalShift9.0cm_wNoise_recon.root -o /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_HorizontalShift9.0cm_wNoise_ana.root",i,i,i);

    /*
    sprintf(CommandR4,"/export/scraid3/data/bquilain/Programs/Programs/INGRID_analysis/app/PMreconRev_KS_woXTalk -r 14000 -s %d -f /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_HShit1.0_wNoise.root -o /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_HShit1.0_wNoise_recon.root",i,i,i);
    sprintf(CommandR40,"/export/scraid3/data/bquilain/Programs/Programs/INGRID_analysis/app/IngAnaVertex_HitInfo -r 14000 -s %d -f /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_HShit1.0_wNoise_recon.root -o /export/scraid3/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_HShit1.0_wNoise_ana.root",i,i,i);
    */


    /*  sprintf(Command1,"/export/scraid2/data/bquilain/Programs/INGRID_analysis/app/PMreconRev_KS_woXTalk -r 14000 -s %d -f /export/scraid2/data/bquilain/MCfiles/Ing_MC_Beam25_BirksCorrectedMIP46_Run1_wNoise.root -o /export/scraid2/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_Run1_wNoise_recon.root",i,i,i);
    sprintf(Command11,"/export/scraid2/data/bquilain/Programs/INGRID_analysis/app/IngAnaVertex_HitInfo -r 14000 -s %d -f /export/scraid2/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_Run1_wNoise_recon.root -o /export/scraid2/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_Run1_wNoise_ana.root",i,i,i);
*/    
//sprintf(Command1,"/export/scraid2/data/bquilain/Programs/INGRID_analysis/app/PMreconRev_KS_woXTalk -r 14000 -s %d -f /export/scraid2/data/bquilain/MCfiles/Ing_MC_Sand%d_BirksCorrectedMIP46_ReWeight_wNoise.root -o /export/scraid2/data/bquilain/MCfiles/Ing_MC_Sand%d_BirksCorrectedMIP46_ReWeight_wNoise_recon.root",i,i,i);
    //sprintf(Command11,"/export/scraid2/data/bquilain/Programs/INGRID_analysis/app/IngAnaVertex_HitInfo -r 14000 -s %d -f /export/scraid2/data/bquilain/MCfiles/Ing_MC_Sand%d_BirksCorrectedMIP46_ReWeight_wNoise_recon.root -o /export/scraid2/data/bquilain/MCfiles/Ing_MC_Sand%d_BirksCorrectedMIP46_ReWeight_wNoise_ana.root",i,i,i);


    //sprintf(Command1,"./PMrecon_HitInfo -r 14000 -s %d -f /export/scraid2/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_wNoise.root -o /export/scraid2/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_wNoise_recon.root",i,i,i);
    //sprintf(Command11,"./IngAnaVertex_HitInfo -r 14000 -s %d -f /export/scraid2/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_wNoise_recon.root -o /export/scraid2/data/bquilain/MCfiles/Ing_MC_Beam%d_BirksCorrectedMIP46_wNoise_ana.root",i,i,i);



    
    //    sprintf(Command1,"./PMrecon_HitInfo -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/Ingrid/Ing_MC_Beam%d_AllPlasticPlane_wNoise.root -o /export/scraid2/data/bquilain/MCfiles/Ingrid/Ing_MC_Beam%d_AllPlasticPlane_wNoise_recon.root",i,i);
    //sprintf(Command11,"./IngAnaVertex_HitInfo -r 14000 -s 0 -f /export/scraid2/data/bquilain/MCfiles/Ingrid/Ing_MC_Beam%d_AllPlasticPlane_wNoise_recon.root -o /export/scraid2/data/bquilain/MCfiles/Ingrid/Ing_MC_Beam%d_AllPlasticPlane_wNoise_ana.root",i,i);
   //sprintf(Command1,"./PMrecon_HitInfo -r 14000 -s 0 -f ../../Ing_MC_Beam%d_.5mmcut_v2.root -o ../../Ing_MC_Beam%d_.5mmcut_v2_recon.root",i,i);
    //sprintf(Command11,"./IngAnaVertex_HitInfo -r 14000 -s 0 -f ../../Ing_MC_Beam%d_.5mmcut_v2_recon.root -o ../../Ing_MC_Beam%d_.5mmcut_v2_ana.root",i,i);
    //sprintf(Command1,"./IngAnaVertex_HitInfo -r 14000 -s 0 -f ../../Ing_MC_Beam%d_NewGainNewMIP9_wNoise_recon.root -o ../../Ing_MC_Beam%d_NewGainNewMIP9_wNoise_ana.root",i,i);

    //sprintf(Command1,"./PMrecon_HitInfo -r 14000 -s 0 -f ../../PM_MC_RealBeam%d_NewGainNewMIP6.root -o ../../PM_MC_RealBeam%d_NewGainNewMIP6_recon.root",i,i);
    //sprintf(Command1,"./IngAnaVertex_HitInfo -r 14000 -s 0 -f ../../PM_MC_RealBeam%d_NewGainNewMIP6_recon.root -o ../../PM_MC_RealBeam%d_NewGainNewMIP6_ana.root",i,i);
    //sprintf(Command1,"./IngAnaVertex_HitInfo -r 14000 -s 0 -f ../../Ing_MC_AllMergedSand_NewGainNewMIP5_wNoise_%d_recon.root -o ../../Ing_MC_AllMergedSand_NewGainNewMIP5_wNoise_%d_ana.root",i,i);
    ofstream Script(Name1);
    if(Script)
      {
	Script<<"#!/bin/bash +x"<<endl;
	  //      Script<<Command1<<endl;
	  //Script<<Command11<<endl;
	/*	Script<<CommandR1<<endl;
	Script<<CommandR10<<endl;
	Script<<CommandR2<<endl;
	Script<<CommandR20<<endl;
	Script<<CommandR3<<endl;
	Script<<CommandR30<<endl;
	*/	Script<<CommandR4<<endl;
	Script<<CommandR40<<endl;
      }


    sprintf(Command2,"Executable = AnaBeam%d.sh",i);
    sprintf(Command3,"Output = condor_AnaBeamlog%d.txt",i);
    sprintf(Command4,"Error = condor_AnaBeamerr%d.txt",i);
    sprintf(Name2,"condorAnaBeam%d.sh",i);
    cout<<Name2<<endl;
    ofstream Condor(Name2);
    if(Condor){
      //cout<<"hello"<<endl;
      Condor<<Command2<<endl
	    <<"Universe        = vanilla"<<endl
	    <<"Rank            = kflops"<<endl
	    <<"Requirements    = CpuIsBusy == FALSE && (machine != \"scbn00\" && machine != \"scbn01.hepnet.scphys.kyoto-u.ac.jp\" )"<<endl
	    <<"Getenv          = True"<<endl
	    <<"Arguments      =  $(Process)"<<endl
	    <<Command3<<endl
	    <<Command4<<endl
	    <<"Notification    = Never"<<endl
	    <<"QUEUE 1"<<endl;
    }
    
  }
  
}
