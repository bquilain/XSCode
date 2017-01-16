CPP 		= g++
CXXFLAGS	= -g -O3 -Wall -fPIC -D_REENTRANT -Wno-deprecated

ROOTCFLAGS	:= $(shell root-config --cflags)
ROOTLIBS     	:= $(shell root-config --libs)
ROOTGLIBS    	:= $(shell root-config --glibs)
CXXFLAGS	+= $(ROOTCFLAGS)

LIB = ${INSTALLREPOSITORY}/libs/libsPM/lib_5.34
#LIB = /home/bquilain/CC0pi_XS/MC/ingmc_IngridRevReWeightFinal/lib_New

LIBS 		= $(ROOTGLIBS) \
          -lm $(LIB)/INGRIDEVENTSUMMARY.so $(LIB)/IngridHitSummary.so \
          $(LIB)/IngridSimHitSummary.so $(LIB)/IngridSimVertexSummary.so \
          $(LIB)/IngridSimParticleSummary.so $(LIB)/BeamInfoSummary.so \
          $(LIB)/IngridBasicReconSummary.so $(LIB)/Ingrid1stReducSummary.so \
          $(LIB)/IngridTrackSummary.so $(LIB)/NeutInfoSummary.so \
          $(LIB)/PMReconSummary.so $(LIB)/PMAnaSummary.so \


SRC = ${INSTALLREPOSITORY}/Reconstruction/src
INC = ${INSTALLREPOSITORY}/Reconstruction/inc

CXXFLAGS += -I$(LIB) -I$(SRC) -I$(INC)

TARGET = IngTClster IngBasicRecon IngReconMuBasic IngGetOnlyBasic \
         IngGetOnlyBasicMC get_daily_evt get_daily_evt_corr \
         PMrecon IngAnaVertex PMAna PMGetOnlyBasic PMGetOnlyBasicMC \
         IngMergePM PMrecon_HitInfo IngAnaVertex_HitInfo PMrecon_HitInfo_Sand PMreconRev PMAna_HitInfo PMAna_HitInfo_Xsec PMAna_woIngrid PMAna_ForXsections PMreconRev2 PEdistBen PMPEinPMsci TrackMatching3D VetoAndFV PMreconRev_woXTalk PMreconRev_KS PMAnaRev_KS_woXTalk PMreconRev_KS_woXTalk SandEventReductionPM PMAnaRevRev PMAnaRevOfficial PMreconRevOfficial PMreconRevOfficial_shortingridtracks

all: $(TARGET)

IngTClster: IngTClster.o

IngBasicRecon: IngBasicRecon.o

IngReconMuBasic: IngReconMuBasic.o

IngGetOnlyBasic: IngGetOnlyBasic.o

IngGetOnlyBasicMC: IngGetOnlyBasicMC.o

get_daily_evt: get_daily_evt.o

get_daily_evt_corr: get_daily_evt_corr.o

IngMergePM: IngMergePM.o

PMrecon: PMrecon.o

IngAnaVertex: IngAnaVertex.o

PMAna: PMAna.o

PMGetOnlyBasic: PMGetOnlyBasic.o

PMGetOnlyBasicMC: PMGetOnlyBasicMC.o

PMrecon_HitInfo: PMrecon_HitInfo.o

IngAnaVertex_HitInfo: IngAnaVertex_HitInfo.o

PMrecon_HitInfo_Sand: PMrecon_HitInfo_Sand.o

PMreconRev: PMreconRev.o

PMAna_HitInfo: PMAna_HitInfo.o

PMAna_HitInfo_Xsec: PMAna_HitInfo_Xsec.o

PMAna_HitInfo_Xsec: PMAna_woIngrid.o

PMAna_ForXsections: PMAna_ForXsections.o

PMreconRev2: PMreconRev2.o

PEdistBen: PEdistBen.o

PMPEinPMsci: PMPEinPMsci.o

TrackMatching3D: TrackMatching3D.o

VetoAndFV: VetoAndFV.o

PMreconRev_woXTalk: PMreconRev_woXTalk.o

PMreconRev_KS: PMreconRev_KS.o

PMAnaRev_KS_woXTalk: PMAnaRev_KS_woXTalk.o

PMreconRev_KS_woXTalk: PMreconRev_KS_woXTalk.o

SandEventReductionPM: SandEventReductionPM.o

PMAnaRevRev: PMAnaRevRev.o

PMAnaRevOfficial: PMAnaRevOfficial.o

PMreconRevOfficial: PMreconRevOfficial.o

PMreconRevOfficial_shortingridtracks: PMreconRevOfficial_shortingridtracks.o

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
