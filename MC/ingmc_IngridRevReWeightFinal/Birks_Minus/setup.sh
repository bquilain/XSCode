#Set ND280 software package PATH
export ND280ROOT=/export/scbn03/data1/kikawa/T2K/work/ROOT/v5r22p00n00
export ND280GEANT=/export/scbn03/data1/kikawa/T2K/work/GEANT/v9r2p01n00
export ND280CLHEP=/export/scbn03/data1/kikawa/T2K/work/CLHEP/v2r0p3


#Set Env various by cmt
source $ND280GEANT/cmt/setup.sh
source $ND280ROOT/cmt/setup.sh
source $ND280CLHEP/cmt/setup.sh

#Set MY Env
export G4WORKDIR=`pwd`
export G4BIN=$G4WORKDIR/bin
export G4TMP=$G4WORKDIR/tmp
