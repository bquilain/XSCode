export INSTALLREPOSITORY=/home/mlicciardi/T2K/work/XSCode
export MCINPUTSTORAGE_WM=/export/scraid4/data2/taichiro/neut
export MCINPUTSTORAGE=${MCINPUTSTORAGE_WM}
export DATAINPUTSTORAGE=/export/scraid2/data/bquilain
# need DATAINPUTSTORAGE_WM
export MCOUTPUTSTORAGE=/export/scraid4/data2/mlicciardi/XSFiles/Prod0/PM
export MCOUTPUTSTORAGE_WM=/export/scraid4/data2/mlicciardi/XSFiles/Prod0/WM
export MCOUTPUTSTORAGE_dev=/export/scraid4/data2/mlicciardi/XSFiles/dev/PM
export MCOUTPUTSTORAGE_dev_WM=/export/scraid4/data2/mlicciardi/XSFiles/dev/WM
export DATAOUTPUTSTORAGE=/export/scraid2/data/bquilain
export INGRIDSOFTWARE=/home/bquilain/T2KSoftware2
#export INGRIDSOFTWARE=/home/mlicciardi/T2K/work/basesoft

#Set ND280 software package PATH
export ND280ROOT=/home/mlicciardi/T2K/work/basesoft/ROOT/v5r22p00n00
export ND280GEANT=/home/mlicciardi/T2K/work/basesoft/GEANT/v9r2p01n00
export ND280CLHEP=/home/mlicciardi/T2K/work/basesoft/CLHEP/v2r0p3

export INGRIDSOFTWARE_INGRID=${INGRIDSOFTWARE}/INGRID
export INGRIDSOFTWARE_INGRIDFORMAT=${INGRIDSOFTWARE}/ingrid_format

export DATAOUTPUTSTORAGE_DST=$DATAOUTPUTSTORAGE/dst
export DATAOUTPUTSTORAGE_CALIB=$DATAOUTPUTSTORAGE/calib_root_file
export DATAOUTPUTSTORAGE_RECONSTRUCTED=$DATAOUTPUTSTORAGE/DataNew
export DATAOUTPUTSTORAGE_ROOT=$DATAOUTPUTSTORAGE/PM

#source $INSTALLREPOSITORY/Run_At_Start.bashrc //ML perso i don't need that
#Set Env various by cmt
source $ND280GEANT/cmt/setup.sh
#source $ND280ROOT/cmt/setup.sh  commented out to keep ROOTSYS at root5.34 while using MC (see GNUmakefile)
source $ND280CLHEP/cmt/setup.sh

export TMVASYS=/home/mlicciardi/T2K/work/basesoft/TMVA-v4.2.0

#Set MY Env
#export G4WORKDIR=$INSTALLREPOSITORY/MC/ingmc_IngridRevReWeightFinal
export G4WORKDIR=$INSTALLREPOSITORY/MC
export G4BIN=$G4WORKDIR/bin
export G4TMP=$G4WORKDIR/tmp
#export G4NEUTRONHPDATA=~/T2KSoftware2/GEANT/v9r2p01n00/Linux-x86_64/G4NDL4.2
