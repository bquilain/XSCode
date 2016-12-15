CPP 		= g++
CXXFLAGS	= -g -O3 -Wall -fPIC -D_REENTRANT -Wno-deprecated

ROOTCFLAGS	:= $(shell root-config --cflags)
ROOTLIBS     	:= $(shell root-config --libs)
ROOTGLIBS    	:= $(shell root-config --glibs)
CXXFLAGS	+= $(ROOTCFLAGS)

#LIB = $(INSTALLREPOSITORY)/CC0pi_XS/MC/ingmc_IngridRevReWeightFinal/lib_New
LIB = $(INSTALLREPOSITORY)/libs/lib_5.34
#lib

LIBS 		= $(ROOTGLIBS) \
          -lm $(LIB)/INGRIDEVENTSUMMARY.so $(LIB)/IngridHitSummary.so \
          $(LIB)/IngridSimHitSummary.so $(LIB)/IngridSimVertexSummary.so \
          $(LIB)/IngridSimParticleSummary.so $(LIB)/BeamInfoSummary.so \
          $(LIB)/IngridBasicReconSummary.so $(LIB)/Ingrid1stReducSummary.so \
          $(LIB)/IngridTrackSummary.so $(LIB)/NeutInfoSummary.so \
          $(LIB)/PMReconSummary.so $(LIB)/PMAnaSummary.so \

INC = inc
SRC= src

CXXFLAGS += -I$(LIB) -I$(SRC) -I$(INC) 

#TARGET= IngAddNoisePMMC CompareCalibrationsPM GeneratePEAngleDistributions HitEfficiency GenerateHitEfficiencyMask Draw_CompareCalibrationsPM GeneratePDFMuCL_Plan XS_CC0pi_Plan SandSelection_Plan CC0piSelection UnfoldingOptimisation DrawXSPlots XSFileGenerator
#TARGET= GenerateHitEfficiencyMask
TARGET= CC0piSelection XS_CC0pi_Plan XSFileGenerator UnfoldingOptimisation UnfoldingOptimisation_Dvt DrawXSPlots_Toys GeneratePDFMuCL_Likelihood GeneratePDFMuCL_Plan IngAddNoisePMMC HitEfficiency CompareCalibrationsPM GeneratePEAngleDistributions GenerateHitEfficiencyMask DrawReconstructedEventsPlots_Dvt FinalMCOutputMaker
#DrawTrueEventsPlots_Dvt

all: $(TARGET)
IngAddNoisePMMC: IngAddNoisePMMC.o
XS_CC0pi_Plan: XS_CC0pi_Plan.o
UnfoldingOptimisation: UnfoldingOptimisation.o
GeneratePDFMuCL_Plan: GeneratePDFMuCL_Plan.o
GeneratePDFMuCL_Likelihood: GeneratePDFMuCL_Likelihood.o
CompareCalibrationsPM: CompareCalibrationsPM.o
GeneratePEAngleDistributions: GeneratePEAngleDistributions.o
HitEfficiency: HitEfficiency.o
GenerateHitEfficiencyMask: GenerateHitEfficiencyMask.o
Draw_CompareCalibrationsPM: Draw_CompareCalibrationsPM.o
SandSelection_Plan: SandSelection_Plan.o
CC0piSelection: CC0piSelection.o
DrawXSPlots: DrawXSPlots.o
DrawXSPlots_Toys: DrawXSPlots_Toys.o
XSFileGenerator: XSFileGenerator.o
UnfoldingOptimisation_Dvt: UnfoldingOptimisation_Dvt.o
DrawReconstructedEventsPlots_Dvt: DrawReconstructedEventsPlots_Dvt.o
DrawTrueEventsPlots_Dvt: DrawTrueEventsPlots_Dvt.o
FinalMCOutputMaker: FinalMCOutputMaker.o

%: %.o
	@echo "Now make $@"
	@$(CPP) -o $@ $< $(CXXFLAGS) $(LIBS) 
	@echo "..Compile done! "

%.o: %.c
	@echo "$<"
	@echo "Start Compiling $<"
	@$(CPP) $(CXXFLAGS) -c $<
	@echo ".. Compiling Object Files $<   --> done"
	@echo ""

clean: 
	@echo "Now Clean Up"
	rm -f $(TARGET) *~ *.o *.o~ core
