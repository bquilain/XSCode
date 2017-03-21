CPP 		= g++
CXXFLAGS	= -g -O3 -Wall -fPIC -D_REENTRANT -Wno-deprecated

ROOTCFLAGS	:= $(shell root-config --cflags)
ROOTLIBS     	:= $(shell root-config --libs)
ROOTGLIBS    	:= $(shell root-config --glibs)
CXXFLAGS	+= $(ROOTCFLAGS)

LIB = ${INSTALLREPOSITORY}/libs/libsWM
LIBSO = $(LIB)/lib_5.34

LIBS 		= $(ROOTGLIBS) \
          -lm $(LIBSO)/INGRIDEVENTSUMMARY.so $(LIBSO)/IngridHitSummary.so \
          $(LIBSO)/IngridSimHitSummary.so $(LIBSO)/IngridSimVertexSummary.so \
          $(LIBSO)/IngridSimParticleSummary.so $(LIBSO)/BeamInfoSummary.so \
          $(LIBSO)/IngridBasicReconSummary.so $(LIBSO)/Ingrid1stReducSummary.so \
          $(LIBSO)/IngridTrackSummary.so $(LIBSO)/NeutInfoSummary.so \
          $(LIBSO)/PMReconSummary.so $(LIBSO)/PMAnaSummary.so \



SRC = ${INSTALLREPOSITORY}/Reconstruction/src
INC = ${INSTALLREPOSITORY}/Reconstruction/inc

CXXFLAGS += -I$(LIB)/src -I$(SRC) -I$(INC)

TARGET = IngTClster IngBasicRecon IngReconMuBasic IngGetOnlyBasic \
         IngGetOnlyBasicMC get_daily_evt get_daily_evt_corr \
         IngAnaVertex PMGetOnlyBasic PMGetOnlyBasicMC \
         IngMergePM IngAddNoisePMMC_new IngAddCutMod3 \
	 Lolirecon LoliAna Loli_addcrosstalk IngAddEasirochit \
	 Loli_addcrosstalk_slit

all: $(TARGET)

IngTClster: IngTClster.o
IngBasicRecon: IngBasicRecon.o
IngReconMuBasic: IngReconMuBasic.o
IngGetOnlyBasic: IngGetOnlyBasic.o
IngGetOnlyBasicMC: IngGetOnlyBasicMC.o
get_daily_evt: get_daily_evt.o
get_daily_evt_corr: get_daily_evt_corr.o
IngMergePM: IngMergePM.o
IngAnaVertex: IngAnaVertex.o
PMGetOnlyBasic: PMGetOnlyBasic.o
PMGetOnlyBasicMC: PMGetOnlyBasicMC.o
IngAddNoisePMMC_new: IngAddNoisePMMC_new.o
IngAddCutMod3: IngAddCutMod3.o
Lolirecon: Lolirecon.o
LoliAna: LoliAna.o
Loli_addcrosstalk: Loli_addcrosstalk.o
Loli_addcrosstalk_slit: Loli_addcrosstalk_slit.o
#Loli_addnoise: Loli_addnoise.o
IngAddEasirochit: IngAddEasirochit.o


%: %.o
	@echo "Now make $@"
	@$(CPP) -o $@ $< $(CXXFLAGS) $(LIBS)
	@echo "..Compile done! "

%.o: %.cc
	@echo "Start Compiling $<"
	@$(CPP) $(CXXFLAGS) -c $<
	@echo ".. Compiling Object Files $<   --> done"
	@echo ""


clean: 
	@echo "Now Clean Up"
	rm -f $(TARGET) *~ *.o *.o~ core
