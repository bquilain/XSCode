###############################################################################

#include $(ROOTSYS)/test/Makefile.arch
#include $(ROOTSYS)/etc/Makefile.arch

SrcSuf	= cxx
ObjSuf	= o
ExeSuf	=
OutPutOpt = -o 

ROOTGLIBS = $(shell root-config --glibs)
CXXFLAGS = $(shell root-config --cflags) -O2 -Wall -fPIC -D_REENTRANT -Wno-deprecated
LDFLAGS = $(shell root-config --glibs --cflags)
LD = g++

LIB = ${INSTALLREPOSITORY}/libs/libsPM/lib_5.34

LIBS 		= $(ROOTGLIBS) \
	-lm $(LIB)/INGRIDEVENTSUMMARY.so $(LIB)/IngridHitSummary.so \
	$(LIB)/IngridSimHitSummary.so $(LIB)/IngridSimVertexSummary.so \
        $(LIB)/IngridSimParticleSummary.so $(LIB)/BeamInfoSummary.so \
        $(LIB)/IngridBasicReconSummary.so $(LIB)/Ingrid1stReducSummary.so \
        $(LIB)/IngridTrackSummary.so $(LIB)/NeutInfoSummary.so \
        $(LIB)/PMReconSummary.so $(LIB)/PMAnaSummary.so \



###############################################################################

SRCS    := event_viewer.$(SrcSuf) event_viewerDict.$(SrcSuf)
OBJS    := $(SRCS:.$(SrcSuf)=.$(ObjSuf)) event_viewerDict.$(ObjSuf)
EXE      = GUI$(ExeSuf)

PROGRAMS      = $(EXE)

###############################################################################

.SUFFIXES: .$(SrcSuf) .$(ObjSuf)

all:        $(PROGRAMS)

$(EXE):  $(OBJS)
	$(LD) $(LDFLAGS) $^ $(ROOTGLIBS) $(LIBS) -Wl  $(OutPutOpt)$@
	@echo "$@ done"

event_viewerDict.cxx: event_viewer.h LinkDef.h
#		@echo $(LIBS) 
		@echo "Generating dictionary $@..."
		@rootcint -f $@ -c -I$(LIB) $^ 

clean:
		rm -f *.o *.d *Dict.* $(EXE)

###############################################################################

.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) $(EXTRAFLAGS) -I$(LIB) -c $<

###############################################################################

