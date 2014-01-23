UNAME := $(shell uname)

SRC=src

ifeq ($(UNAME), Darwin)
	FLAG=-O2 -Wall -I `root-config --incdir` -I /opt/local/include/ -I ${SRC}
	LIBS = -L ${ROOTSYS}/lib -L /opt/local/lib/
   BOOSTLIBS = -lboost_serialization-mt -lboost_program_options-mt 
endif
ifeq ($(UNAME), Linux)
	FLAG=-O3 -Wall -Werror=type-limits -I `root-config --incdir` -I ${SRC} -I `scram tool tag boost INCLUDE`
	LIBS = -L ${ROOTSYS}/lib -L `scram tool tag boost LIBDIR`
   BOOSTLIBS = -lboost_serialization -lboost_program_options 
endif

ifeq ($(CUDA_ENABLED),true)
        OBJECTS=AMSimulation.o SuperStrip.o Hit.o Pattern.o PatternLayer.o GradedPattern.o PatternTrunk.o PatternTree.o PatternGenerator.o Sector.o SectorTree.o CMSPatternLayer.o Segment.o Module.o Ladder.o Layer.o Detector.o PatternFinder.o Track.o TrackFitter.o FitParams.o PrincipalTrackFitter.o PrincipalFitGenerator.o MultiDimFitData.o KarimakiTrackFitter.o HoughFitter.o ComputerHough.o libhoughCPU.o gpu.o
else
        OBJECTS=AMSimulation.o SuperStrip.o Hit.o Pattern.o PatternLayer.o GradedPattern.o PatternTrunk.o PatternTree.o PatternGenerator.o Sector.o SectorTree.o CMSPatternLayer.o Segment.o Module.o Ladder.o Layer.o Detector.o PatternFinder.o Track.o TrackFitter.o FitParams.o PrincipalTrackFitter.o PrincipalFitGenerator.o MultiDimFitData.o KarimakiTrackFitter.o HoughFitter.o ComputerHough.o libhoughCPU.o
endif

AMSimulation:$(OBJECTS)
	g++ -o AMSimulation $(OBJECTS) ${LIBS} ${BOOSTLIBS} -lCore -lCint -lRIO -lHist -lTree -lMatrix

clean:	
	rm -rf *.o;rm -f ${SRC}/*~;rm -f ${SRC}/*#

Hit.o:${SRC}/Hit.h ${SRC}/Hit.cc
	g++ -c ${FLAG} ${SRC}/Hit.cc

Segment.o:${SRC}/SuperStrip.h ${SRC}/Segment.h ${SRC}/Segment.cc
	g++ -c ${FLAG} ${SRC}/Segment.cc

Module.o:${SRC}/Segment.h ${SRC}/Module.h ${SRC}/Module.cc
	g++ -c ${FLAG} ${SRC}/Module.cc

Ladder.o:${SRC}/Module.h ${SRC}/Ladder.h ${SRC}/Ladder.cc
	g++ -c ${FLAG} ${SRC}/Ladder.cc

Layer.o:${SRC}/Ladder.h ${SRC}/Layer.h ${SRC}/Layer.cc
	g++ -c ${FLAG} ${SRC}/Layer.cc

Detector.o:${SRC}/Layer.h ${SRC}/Detector.h ${SRC}/Detector.cc
	g++ -c ${FLAG} ${SRC}/Detector.cc

SuperStrip.o:${SRC}/SuperStrip.h ${SRC}/SuperStrip.cc
	g++ -c ${FLAG} ${SRC}/SuperStrip.cc

PatternLayer.o:${SRC}/PatternLayer.h ${SRC}/PatternLayer.cc
	g++ -c ${FLAG} ${SRC}/PatternLayer.cc

CMSPatternLayer.o:${SRC}/PatternLayer.h ${SRC}/CMSPatternLayer.h ${SRC}/CMSPatternLayer.cc
	g++ -c ${FLAG} ${SRC}/CMSPatternLayer.cc

Sector.o:${SRC}/Sector.h ${SRC}/Sector.cc
	g++ -c ${FLAG} ${SRC}/Sector.cc

SectorTree.o:${SRC}/SectorTree.h ${SRC}/SectorTree.cc
	g++ -c ${FLAG} ${SRC}/SectorTree.cc

Pattern.o:${SRC}/Pattern.h ${SRC}/Pattern.cc
	g++ -c ${FLAG} ${SRC}/Pattern.cc

GradedPattern.o:${SRC}/GradedPattern.h ${SRC}/GradedPattern.cc
	g++ -c ${FLAG} ${SRC}/GradedPattern.cc

PatternTrunk.o:${SRC}/PatternTrunk.h ${SRC}/PatternTrunk.cc
	g++ -c ${FLAG} ${SRC}/PatternTrunk.cc

PatternTree.o:${SRC}/PatternTree.h ${SRC}/PatternTree.cc
	g++ -c ${FLAG} ${SRC}/PatternTree.cc

PatternGenerator.o:${SRC}/PatternGenerator.h ${SRC}/PatternGenerator.cc
	g++ -c ${FLAG} ${SRC}/PatternGenerator.cc

PatternFinder.o:${SRC}/PatternFinder.h ${SRC}/PatternFinder.cc
	g++ -c ${FLAG} ${SRC}/PatternFinder.cc

Track.o:${SRC}/Track.h ${SRC}/Track.cc
	g++ -c ${FLAG} ${SRC}/Track.cc

TrackFitter.o:${SRC}/TrackFitter.h ${SRC}/TrackFitter.cc
	g++ -c ${FLAG} ${SRC}/TrackFitter.cc

FitParams.o:${SRC}/FitParams.h ${SRC}/FitParams.cc
	g++ -c ${FLAG} ${SRC}/FitParams.cc

PrincipalTrackFitter.o:${SRC}/PrincipalTrackFitter.h ${SRC}/PrincipalTrackFitter.cc
	g++ -c ${FLAG} ${SRC}/PrincipalTrackFitter.cc

PrincipalFitGenerator.o:${SRC}/PrincipalFitGenerator.h ${SRC}/PrincipalFitGenerator.cc
	g++ -c ${FLAG} ${SRC}/PrincipalFitGenerator.cc

MultiDimFitData.o:${SRC}/MultiDimFitData.h ${SRC}/MultiDimFitData.cc
	g++ -c ${FLAG} ${SRC}/MultiDimFitData.cc

KarimakiTrackFitter.o:${SRC}/KarimakiTrackFitter.h ${SRC}/KarimakiTrackFitter.cc
	g++ -c ${FLAG} ${SRC}/KarimakiTrackFitter.cc

HoughFitter.o:${SRC}/HoughFitter.h ${SRC}/HoughFitter.cc
	g++ -c ${FLAG} ${SRC}/HoughFitter.cc

ComputerHough.o:${SRC}/ComputerHough.h ${SRC}/ComputerHough.cc $(SRC)/libhoughStruct.h $(SRC)/HoughStruct.h libhoughCPU.o
	g++ -c ${FLAG} ${SRC}/ComputerHough.cc

libhoughCPU.o:${SRC}/libhoughCPU.h ${SRC}/libhoughCPU.c
	g++ -c ${FLAG} ${SRC}/libhoughCPU.c

AMSimulation.o:${SRC}/AMSimulation.cc
	g++ -c ${FLAG} ${SRC}/AMSimulation.cc

doc:doxygen.cfg
	doxygen doxygen.cfg

AnalysePatternsPU:Analysis/AnalysePatternsPU.cc
	g++ -c ${FLAG} Analysis/AnalysePatternsPU.cc;g++ -o AnalysePatternsPU AnalysePatternsPU.o ${LIBS} `root-config --libs`

ExtractTTree:Analysis/ExtractTTree.cc
	g++ -c ${FLAG} Analysis/ExtractTTree.cc;g++ -o ExtractTTree ExtractTTree.o ${LIBS} `root-config --libs`
