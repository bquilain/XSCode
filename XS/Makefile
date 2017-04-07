CPP 		= g++
CXXFLAGS	= -g -O3 -Wall -fPIC -D_REENTRANT -Wno-deprecated

ROOTCFLAGS	:= $(shell root-config --cflags)
ROOTLIBS     	:= $(shell root-config --libs)
ROOTGLIBS    	:= $(shell root-config --glibs)
CXXFLAGS	+= $(ROOTCFLAGS)

<<<<<<< HEAD
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
=======
#LIB = $(INSTALLREPOSITORY)/libs/libsPM/lib_5.34

LIBS 		= $(ROOTGLIBS) \
>>>>>>> d3a82e5c713324d10e99b133c2c5ed1e3f825d05

INC = inc
SRC= src

CXXFLAGS += -I$(SRC) -I$(INC) 

<<<<<<< HEAD
#TARGET= IngAddNoisePMMC CompareCalibrationsPM GeneratePEAngleDistributions HitEfficiency GenerateHitEfficiencyMask Draw_CompareCalibrationsPM GeneratePDFMuCL_Plan XS_CC0pi_Plan SandSelection_Plan CC0piSelection UnfoldingOptimisation DrawXSPlots XSFileGenerator
#TARGET= GenerateHitEfficiencyMask
TARGET= CC0piSelection_Dvt XS_CC0pi_Plan_Dvt XSFileGenerator UnfoldingOptimisation UnfoldingOptimisation_Dvt DrawXSPlots_Toys GeneratePDFMuCL_Likelihood_Dvt GeneratePDFMuCL_Plan IngAddNoisePMMC HitEfficiency CompareCalibrationsPM GeneratePEAngleDistributions GenerateHitEfficiencyMask DrawReconstructedEventsPlots_Dvt FinalMCOutputMaker
#DrawTrueEventsPlots_Dvt

all: $(TARGET)
IngAddNoisePMMC: IngAddNoisePMMC.o
XS_CC0pi_Plan_Dvt: XS_CC0pi_Plan_Dvt.o
UnfoldingOptimisation: UnfoldingOptimisation.o
GeneratePDFMuCL_Plan: GeneratePDFMuCL_Plan.o
GeneratePDFMuCL_Likelihood_Dvt: GeneratePDFMuCL_Likelihood_Dvt.o
CompareCalibrationsPM: CompareCalibrationsPM.o
=======
TARGET=XSFileGenerator FinalMCOutputMaker CC0piPlots CC0piSelection CC0piSelection_Likelihood Draw_CompareCalibrationsPM DrawXSPlots DrawXSPlots_Toys GeneratePEAngleDistributions SandSelection_Plan SequentialSand

all: $(TARGET)

>>>>>>> d3a82e5c713324d10e99b133c2c5ed1e3f825d05
GeneratePEAngleDistributions: GeneratePEAngleDistributions.o
Draw_CompareCalibrationsPM: Draw_CompareCalibrationsPM.o
SandSelection_Plan: SandSelection_Plan.o
<<<<<<< HEAD
CC0piSelection_Dvt: CC0piSelection_Dvt.o
=======
CC0piSelection: CC0piSelection.o
CC0piSelection_Likelihood: CC0piSelection_Likelihood.o
>>>>>>> d3a82e5c713324d10e99b133c2c5ed1e3f825d05
DrawXSPlots: DrawXSPlots.o
DrawXSPlots_Toys: DrawXSPlots_Toys.o
CC0piPlots: CC0piPlots.o
XSFileGenerator: XSFileGenerator.o
FinalMCOutputMaker: FinalMCOutputMaker.o
SequentialSand: SequentialSand.o

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
