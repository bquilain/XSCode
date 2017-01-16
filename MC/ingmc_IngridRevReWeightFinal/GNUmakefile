###
CPPFLAGS += -Df2cFortran

###
CERN_ROOT = /cern/pro
CERNINC = $(CERN_ROOT)/include
CERNLIB = -L$(CERN_ROOT)/lib 
CERNLIBS = -lmathlib -lpacklib -lgraflib -lgrafX11 -lpacklib -lkernlib -lmathlib
CPPFLAGS += -I$(CERNINC)
EXTRALIBS += $(CERNLIBS)

###
MYTEMP = -L/usr/X11R6/lib -L/usr/lib/gcc/x86_64-redhat-linux/3.4.6/ -lg2c
EXTRALIBS += $(MYTEMP)

DATALIBDIR = $(INSTALLREPOSITORY)/libs/libsPM/lib_5.22
CPPFLAGS += -I$(DATALIBDIR)

MYLIB = -lm $(DATALIBDIR)/INGRIDEVENTSUMMARY.so $(DATALIBDIR)/IngridHitSummary.so $(DATALIBDIR)/IngridSimHitSummary.so $(DATALIBDIR)/IngridSimVertexSummary.so $(DATALIBDIR)/IngridSimParticleSummary.so $(DATALIBDIR)/BeamInfoSummary.so $(DATALIBDIR)/Ingrid1stReducSummary.so $(DATALIBDIR)/IngridBasicReconSummary.so $(DATALIBDIR)/IngridTrackSummary.so $(DATALIBDIR)/NeutInfoSummary.so $(DATALIBDIR)/PMReconSummary.so $(DATALIBDIR)/PMAnaSummary.so
EXTRALIBS += -L$(MYLIBDIR) $(MYLIB)

### following to keep ROOTSYS to 5.34 while using the MC
ROOTCINT=YES
ND280ROOTSYS=$(ND280ROOT)/Linux-x86_64/root
ROOTLIBS =  $(shell $(ND280ROOTSYS)/bin/root-config --glibs)
CPPFLAGS += $(shell $(ND280ROOTSYS)/bin/root-config --cflags)
EXTRALIBS += $(ROOTLIBS)

###
name := IngMC
G4TARGET := $(name)
G4EXLIB := true

.PHONY: all
#all: lib bin
all: bin

include $(G4INSTALL)/config/binmake.gmk
#-DCERNLIB_TYPE -Dextname -Dlynx -I$(CERNINC)

LDFLAGS += $(CERNLIB)
#LDLIBS += $(CERNLIBS)

CXXFLAGS_WITHOUT_O := $(filter-out -O% , $(CXXFLAGS)) 
CXXFLAGS_WITHOUT_O := $(filter-out +O% , $(CXXFLAGS_WITHOUT_O))
