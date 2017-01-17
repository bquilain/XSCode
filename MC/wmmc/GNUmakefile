###
CPPFLAGS += -Df2cFortran

###
CERN_ROOT = /cern/pro
#changed due to Suzuki-san's mail
#CERN_ROOT = /home/t2k/kentsuzu/cern/cernlib/2005
CERNINC = $(CERN_ROOT)/include
CERNLIB = -L$(CERN_ROOT)/lib 
CERNLIBS = -lmathlib -lpacklib -lgraflib -lgrafX11 -lpacklib -lkernlib -lmathlib
CPPFLAGS += -I$(CERNINC)
EXTRALIBS += $(CERNLIBS)

###
MYTEMP = -L/usr/X11R6/lib -L/usr/lib/gcc/x86_64-redhat-linux/3.4.6/ -lg2c
EXTRALIBS += $(MYTEMP)

DATALIBDIR = $(INSTALLREPOSITORY)/libs/libsWM
DATALIBSO = $(DATALIBDIR)/lib_5.22
CPPFLAGS += -I$(DATALIBDIR)/src

MYLIB = -lm $(DATALIBSO)/INGRIDEVENTSUMMARY.so $(DATALIBSO)/IngridHitSummary.so $(DATALIBSO)/IngridSimHitSummary.so $(DATALIBSO)/IngridSimVertexSummary.so $(DATALIBSO)/IngridSimParticleSummary.so $(DATALIBSO)/BeamInfoSummary.so $(DATALIBSO)/Ingrid1stReducSummary.so $(DATALIBSO)/IngridBasicReconSummary.so $(DATALIBSO)/IngridTrackSummary.so $(DATALIBSO)/NeutInfoSummary.so $(DATALIBSO)/PMReconSummary.so $(DATALIBSO)/PMAnaSummary.so

EXTRALIBS += -L$(MYLIBDIR) $(MYLIB)

### modif ML 2016/10/10
ROOTCINT=YES
ND280ROOTSYS=$(ND280ROOT)/Linux-x86_64/root
ROOTLIBS =  $(shell $(ND280ROOTSYS)/bin/root-config --glibs)
CPPFLAGS += $(shell $(ND280ROOTSYS)/bin/root-config --cflags)
EXTRALIBS += $(ROOTLIBS)

###
name := IngWMMC
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
