ROOTINC = -I`root-config --incdir`
CLASSDIR = .
CLASSINC= -I${CLASSDIR}
INCDIRS = ${ROOTINC}${CLASSINC}

CXXFLAGS  = `root-config --cflags`
COPTFLAGS = -g

SOBJS = NeutInfoSummary.so INGRIDEVENTSUMMARY.so IngridHitSummary.so IngridSimHitSummary.so IngridSimVertexSummary.so IngridSimParticleSummary.so BeamInfoSummary.so IngridBasicReconSummary.so PMReconSummary.so PMAnaSummary.so Ingrid1stReducSummary.so IngridTrackSummary.so

OBJS  =  ${SOBJS} 

.SUFFIXES:      .so

GENROOTSO = env COPTFLAGS="${COPTFLAGS}" INCDIRS="${INCDIRS}" \
                        ./gen_root_so.sh


.cc.o:
	$(CXX) -c $(COPTFLAGS) $(CXXFLAGS) $(INCDIRS) -o $@ $<

.cc.so:
	$(GENROOTSO) $(basename $<)

.F.o:
	$(FC) -c $(FOPTFLAGS) $(FINCDIRS) -o $@ $<

all: $(OBJS)

clean:
	$(RM) -f *.o *~ *LinkDef.h *Dict.* ${OBJS} chkfort
