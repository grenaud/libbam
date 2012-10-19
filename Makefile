
CXX      = g++
BAMTOOLS= /mnt/solexa/bin/bamtools-1.0.2/
LIBGAB   = /home/gabriel_renaud/lib/

CXXFLAGS  = -lm -O3 -Wall -I${LIBGAB}   -I/include/ -c
LDFLAGS  += -lbamtools
LDLIBS   += ${BAMTOOLS}/lib/libbamtools.a -lm -lz


all: allFailqc allPassqc cutDeaminated decrQualDeaminated decrQualDeaminatedDoubleStranded failQualPair filterDeaminated removeRG retrieveRG subSampleBAM transBAM transBAMperRead ReconsReferenceBAM.o




allFailqc.o: allFailqc.cpp
	${CXX} ${CXXFLAGS} allFailqc.cpp

allFailqc: allFailqc.o  ${LIBGAB}utils.o
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)

allPassqc.o: allPassqc.cpp
	${CXX} ${CXXFLAGS} allPassqc.cpp

allPassqc: allPassqc.o  ${LIBGAB}utils.o
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)

cutDeaminated.o: cutDeaminated.cpp
	${CXX} ${CXXFLAGS} cutDeaminated.cpp

cutDeaminated: cutDeaminated.o  ${LIBGAB}utils.o
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)

decrQualDeaminated.o: decrQualDeaminated.cpp
	${CXX} ${CXXFLAGS} decrQualDeaminated.cpp

decrQualDeaminated: decrQualDeaminated.o  ${LIBGAB}utils.o
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)

decrQualDeaminatedDoubleStranded.o: decrQualDeaminatedDoubleStranded.cpp
	${CXX} ${CXXFLAGS} decrQualDeaminatedDoubleStranded.cpp

decrQualDeaminatedDoubleStranded: decrQualDeaminatedDoubleStranded.o  ${LIBGAB}utils.o
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)

failQualPair.o: failQualPair.cpp
	${CXX} ${CXXFLAGS} failQualPair.cpp

failQualPair: failQualPair.o  ${LIBGAB}utils.o
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)

filterDeaminated.o: filterDeaminated.cpp
	${CXX} ${CXXFLAGS} filterDeaminated.cpp

filterDeaminated: filterDeaminated.o  ${LIBGAB}utils.o ReconsReferenceBAM.o
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)

ReconsReferenceBAM.o: ReconsReferenceBAM.cpp
	${CXX} ${CXXFLAGS} ReconsReferenceBAM.cpp

removeRG.o: removeRG.cpp
	${CXX} ${CXXFLAGS} removeRG.cpp

removeRG: removeRG.o  ${LIBGAB}utils.o
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)

retrieveRG.o: retrieveRG.cpp
	${CXX} ${CXXFLAGS} retrieveRG.cpp

retrieveRG: retrieveRG.o  ${LIBGAB}utils.o
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)

subSampleBAM.o: subSampleBAM.cpp
	${CXX} ${CXXFLAGS} subSampleBAM.cpp

subSampleBAM: subSampleBAM.o  ${LIBGAB}utils.o
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)

transBAM.o: transBAM.cpp
	${CXX} ${CXXFLAGS} transBAM.cpp

transBAM: transBAM.o  ${LIBGAB}utils.o ReconsReferenceBAM.o
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)

transBAMperRead.o: transBAMperRead.cpp
	${CXX} ${CXXFLAGS} transBAMperRead.cpp

transBAMperRead: transBAMperRead.o  ${LIBGAB}utils.o ReconsReferenceBAM.o
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)




clean :
	rm -f  allFailqc allPassqc cutDeaminated decrQualDeaminated decrQualDeaminatedDoubleStranded failQualPair filterDeaminated getCtrlReadsBAM removeRG retrieveRG subSampleBAM transBAM transBAMperRead ReconsReferenceBAM.o

