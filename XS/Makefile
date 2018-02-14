CPP 		= g++
CXXFLAGS	= -g -O3 -Wall -fPIC -D_REENTRANT -Wno-deprecated

ROOTCFLAGS	:= $(shell root-config --cflags)
ROOTLIBS     	:= $(shell root-config --libs)
ROOTGLIBS    	:= $(shell root-config --glibs)
CXXFLAGS	+= $(ROOTCFLAGS)

TMVACFLAGS       = -I$(TMVASYS)/
CXXFLAGS        += $(TMVACFLAGS)
TMVALIBS         = -L$(TMVASYS)/lib $(TMVASYS)/lib/libTMVA.1.so $(shell root-config --libs) -lMLP -lTreePlayer -lMinuit

#LIB = $(INSTALLREPOSITORY)/CC0pi_XS/MC/ingmc_IngridRevReWeightFinal/lib_New
LIB = $(INSTALLREPOSITORY)/libs/libsPM/lib_5.34
#lib

LIBS 		= $(ROOTGLIBS) $(TMVALIBS)\
          -lm $(LIB)/INGRIDEVENTSUMMARY.so $(LIB)/IngridHitSummary.so \
          $(LIB)/IngridSimHitSummary.so $(LIB)/IngridSimVertexSummary.so \
          $(LIB)/IngridSimParticleSummary.so $(LIB)/BeamInfoSummary.so \
          $(LIB)/IngridBasicReconSummary.so $(LIB)/Ingrid1stReducSummary.so \
          $(LIB)/IngridTrackSummary.so $(LIB)/NeutInfoSummary.so \
          $(LIB)/PMReconSummary.so $(LIB)/PMAnaSummary.so \

INC = inc
SRC= src

CXXFLAGS += -I$(LIB) -I$(SRC) -I$(INC)

TARGET=XSFileGenerator FinalMCOutputMaker CC0piPlots CC0piSelection Draw_CompareCalibrationsPM DrawXSPlots DrawXSPlots_Toys GeneratePEAngleDistributions XS_CC0pi_Plan DrawReconstructedEventsPlots_Dvt UnfoldingOptimisation_Dvt SequentialSand SandSelection_Plan GeneratePDFMuCL_Likelihood ShrinkXSFormat GenieReweight

all: $(TARGET)


GeneratePEAngleDistributions: GeneratePEAngleDistributions.o
Draw_CompareCalibrationsPM: Draw_CompareCalibrationsPM.o
SandSelection_Plan: SandSelection_Plan.o
CC0piSelection: CC0piSelection.o
CC0piSelection_Likelihood: CC0piSelection_Likelihood.o
DrawXSPlots: DrawXSPlots.o
DrawXSPlots_Toys: DrawXSPlots_Toys.o
CC0piPlots: CC0piPlots.o
XSFileGenerator: XSFileGenerator.o
FinalMCOutputMaker: FinalMCOutputMaker.o
SequentialSand: SequentialSand.o
XS_CC0pi_Plan: XS_CC0pi_Plan.o
DrawReconstructedEventsPlots_Dvt: DrawReconstructedEventsPlots_Dvt.o
UnfoldingOptimisation_Dvt: UnfoldingOptimisation_Dvt.o
GeneratePDFMuCL_Likelihood: GeneratePDFMuCL_Likelihood.o
ShrinkXSFormat: ShrinkXSFormat.o
GenieReweight: GenieReweight.o	

%: %.o
	@echo "Now make $@"
	@$(CPP) -o $@ $< $(CXXFLAGS) $(LIBS) 
	@echo "..Compile done! "

%.o: %.c
	@echo "$<"
	@echo "Start Compiling $<"
	@$(CPP) $(CXXFLAGS) -c $< -o $@
	@echo ".. Compiling Object Files $<   --> done"
	@echo ""

clean: 
	@echo "Now Clean Up"
	rm -f $(TARGET) *~ *.o *.o~ core
