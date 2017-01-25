CPP 		= g++
CXXFLAGS	= -g -O3 -Wall -fPIC -D_REENTRANT -Wno-deprecated

ROOTCFLAGS	:= $(shell root-config --cflags)
ROOTLIBS     	:= $(shell root-config --libs)
ROOTGLIBS    	:= $(shell root-config --glibs)
CXXFLAGS	+= $(ROOTCFLAGS)

#LIB = $(INSTALLREPOSITORY)/libs/libsPM/lib_5.34

LIBS 		= $(ROOTGLIBS) \

INC = inc
SRC= src

CXXFLAGS += -I$(SRC) -I$(INC) 

TARGET=XSFileGenerator FinalMCOutputMaker CC0piPlots CC0piSelection CC0piSelection_Likelihood Draw_CompareCalibrationsPM DrawXSPlots DrawXSPlots_Toys GeneratePEAngleDistributions SandSelection_Plan SequentialSand

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
