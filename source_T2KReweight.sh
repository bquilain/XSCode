#!/bin/bash
# ML 2016/12/15 not modified to fit in my own environment yet

export T2KREWEIGHTREPOSITORY=/home/bquilain/T2KSoftware2/T2KReWeight/v1r23

export CVSROOT=:ext:anoncvs@repo.t2k.org:/home/trt2kmgr/T2KRepository
export CVS_RSH=ssh 
unset CVS_SERVER
export CERN_DIR=/home/bquilain/T2KSoftware2/CERNLIB/2005
export CERN_ROOT=$CERN_DIR
export CERN_LIB=$CERN_DIR/lib
export CERN_INC=$CERN_DIR/include
export CERN=/home/bquilain/T2KSoftware2/CERNLIB
export CERN_LEVEL=2005
export PATH=$CERN_ROOT/bin:$PATH
export LD_LIBRARY_PATH=$CERN_ROOT/lib:$LD_LIBRARY_PATH;

export ROOTSYS=/home/bquilain/T2KSoftware2/ROOT/v5r34p09n03/Linux-x86_64
export PATH=$ROOTSYS/bin:$PATH
export LD_LIBRARY_PATH=$ROOTSYS/lib/root:$LD_LIBRARY_PATH
export MANPATH=${ROOTSYS}/share/man/:${MANPATH}
export NEUT_ROOT=/home/bquilain/T2KSoftware2/neut_5.3.3_v1r21
export PATH=$NEUT_ROOT/src/neutsmpl/bin:$PATH
export LD_LIBRARY_PATH=$NEUT_ROOT/src/reweight:$LD_LIBRARY_PATH;

export NIWG=/home/bquilain/T2KSoftware2/NIWGReWeight
export LD_LIBRARY_PATH=${NIWG}:$LD_LIBRARY_PATH;
export NIWGREWEIGHT_INPUTS=${NIWG}/inputs

export JREWEIGHT=/home/bquilain/T2KSoftware2/JReWeight
export JREWEIGHT_INPUTS=$JREWEIGHT/inputs
export LD_LIBRARY_PATH=$JREWEIGHT:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$OAANALYSISLIBS:$LD_LIBRARY_PATH
export JNUBEAM=$JREWEIGHT

export T2KREWEIGHT=/home/bquilain/T2KSoftware2/T2KReWeight/v1r23
export PATH=$T2KREWEIGHT/bin:$PATH
export LD_LIBRARY_PATH=$T2KREWEIGHT/lib:$LD_LIBRARY_PATH
