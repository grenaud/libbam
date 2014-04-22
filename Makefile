
CXX      = g++ 
BAMTOOLS= /mnt/solexa/bin/bamtools-2.2.2/
LIBGAB   = /home/gabriel_renaud/lib/

CXXFLAGS  = -lm -O3 -Wall -I${LIBGAB} -I${LIBGAB}/VCFparser/ -I${LIBGAB}/VCFparser/gzstream/  -I/home/gabriel_renaud/Software/tabix-0.2.6/ -I${BAMTOOLS}/include/ -c
#LDFLAGS  += -lbamtools
LDLIBS   += ${BAMTOOLS}/lib/libbamtools.a ${LIBGAB}/PutProgramInHeader.o  -lm -lz 


all: allFailqc allPassqc cutDeaminated decrQualDeaminated decrQualDeaminatedDoubleStranded failQualPair filterDeaminated filterDeaminatedVCF removeRG retrieveRG subSampleBAM transBAM transBAMperRead  filterHighEditDistance.o filterHighEditDistance  editDist removeUnalignedANDWrongCigar dumpLoneMates retrieveMapped_single_and_ProperlyPair setAsUnpaired compareRG  filterDeaminatedFasta filterDeaminatedVCFpreload filterDeaminatedVCFpreload1000g filterDeaminatedpreload1000g subsamplebamFixedNumber subsamplebamFixedNumberPair addRGinHeaderHack splitByChr splitByRG filterEditDist bamCat tallyByRG addRG removeTagsMapping addRG_CTEAM baseQualScorePerCycle insertSize singleAndFirstMate cutReadsDistribution cutStart


%.o: %.cpp
	${CXX} ${CXXFLAGS} $^ -o $@


compareRG: compareRG.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

addRG: addRG.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

cutStart: cutStart.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

addRG_CTEAM: addRG_CTEAM.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

removeTagsMapping: removeTagsMapping.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

cutReadsDistribution: cutReadsDistribution.o  ${LIBGAB}utils.o ${LIBGAB}VCFparser/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) 

filterEditDist: filterEditDist.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 


bamCat: bamCat.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 


dumpLoneMates: dumpLoneMates.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 


removeUnalignedANDWrongCigar: removeUnalignedANDWrongCigar.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 


retrieveMapped_single_and_ProperlyPair: retrieveMapped_single_and_ProperlyPair.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 


editDist: editDist.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

insertSize: insertSize.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 


allFailqc: allFailqc.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

singleAndFirstMate: singleAndFirstMate.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 


setAsUnpaired: setAsUnpaired.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 



allPassqc: allPassqc.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 


splitByChr: splitByChr.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

splitByRG: splitByRG.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

tallyByRG: tallyByRG.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

cutDeaminated: cutDeaminated.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 


decrQualDeaminated: decrQualDeaminated.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 


decrQualDeaminatedDoubleStranded: decrQualDeaminatedDoubleStranded.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

failQualPair.o: failQualPair.cpp
	${CXX} ${CXXFLAGS} failQualPair.cpp

failQualPair: failQualPair.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 


filterHighEditDistance: filterHighEditDistance.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 


baseQualScorePerCycle: baseQualScorePerCycle.o  ${LIBGAB}utils.o ${LIBGAB}ReconsReferenceBAM.o
	${CXX} -o $@ $^ $(LDLIBS) 


filterDeaminatedVCF: filterDeaminatedVCF.o  ${LIBGAB}utils.o ${LIBGAB}/VCFparser/VCFreader.o ${LIBGAB}/VCFparser/gzstream/libgzstream.a  ${LIBGAB}/VCFparser/SimpleVCF.o ${LIBGAB}/VCFparser/CoreVCF.o ${LIBGAB}/VCFparser/ReadTabix.o ${LIBGAB}ReconsReferenceBAM.o /home/gabriel_renaud/Software/tabix-0.2.6/libtabix.a
	${CXX} -o $@ $^ $(LDLIBS) 

filterDeaminatedVCFpreload: filterDeaminatedVCFpreload.o  ${LIBGAB}utils.o ${LIBGAB}/VCFparser/VCFreader.o ${LIBGAB}/VCFparser/gzstream/libgzstream.a  ${LIBGAB}/VCFparser/SimpleVCF.o ${LIBGAB}/VCFparser/CoreVCF.o ${LIBGAB}/VCFparser/ReadTabix.o ${LIBGAB}ReconsReferenceBAM.o /home/gabriel_renaud/Software/tabix-0.2.6/libtabix.a 
	${CXX} -o $@ $^ $(LDLIBS) 

filterDeaminatedVCFpreload1000g: filterDeaminatedVCFpreload1000g.o  ${LIBGAB}utils.o ${LIBGAB}/VCFparser/VCFreader.o ${LIBGAB}/VCFparser/gzstream/libgzstream.a  ${LIBGAB}/VCFparser/SimpleVCF.o ${LIBGAB}/VCFparser/CoreVCF.o ${LIBGAB}/VCFparser/ReadTabix.o ${LIBGAB}ReconsReferenceBAM.o /home/gabriel_renaud/Software/tabix-0.2.6/libtabix.a
	${CXX} -o $@ $^ $(LDLIBS) 

filterDeaminatedpreload1000g: filterDeaminatedpreload1000g.o  ${LIBGAB}utils.o ${LIBGAB}/VCFparser/VCFreader.o ${LIBGAB}/VCFparser/gzstream/libgzstream.a  ${LIBGAB}/VCFparser/SimpleVCF.o ${LIBGAB}/VCFparser/CoreVCF.o ${LIBGAB}/VCFparser/ReadTabix.o ${LIBGAB}ReconsReferenceBAM.o /home/gabriel_renaud/Software/tabix-0.2.6/libtabix.a
	${CXX} -o $@ $^ $(LDLIBS) 

filterDeaminatedFasta: filterDeaminatedFasta.o  ${LIBGAB}utils.o ${LIBGAB}ReconsReferenceBAM.o ${LIBGAB}VCFparser/FastQParser.o ${LIBGAB}VCFparser/FastQObj.o ${LIBGAB}/VCFparser/gzstream/libgzstream.a 
	${CXX} -o $@ $^ $(LDLIBS) 

removeRG: removeRG.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

addRGinHeaderHack: addRGinHeaderHack.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 


retrieveRG: retrieveRG.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 


subSampleBAM: subSampleBAM.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

subsamplebamFixedNumber: subsamplebamFixedNumber.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

subsamplebamFixedNumberPair: subsamplebamFixedNumberPair.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

transBAM: transBAM.o  ${LIBGAB}utils.o ${LIBGAB}ReconsReferenceBAM.o
	${CXX} -o $@ $^ $(LDLIBS) 


transBAMperRead: transBAMperRead.o  ${LIBGAB}utils.o ${LIBGAB}ReconsReferenceBAM.o
	${CXX} -o $@ $^ $(LDLIBS) 




clean :
	rm -f  allFailqc allPassqc cutDeaminated decrQualDeaminated decrQualDeaminatedDoubleStranded failQualPair filterDeaminated filterDeaminatedVCF getCtrlReadsBAM removeRG retrieveRG subSampleBAM transBAM transBAMperRead ${LIBGAB}ReconsReferenceBAM.o filterHighEditDistance.o filterHighEditDistance  removeUnalignedANDWrongCigar.o removeUnalignedANDWrongCigar retrieveMapped_single_and_ProperlyPair.o retrieveMapped_single_and_ProperlyPair setAsUnpaired filterDeaminatedFasta compareRG filterDeaminatedFasta  compareRG filterDeaminatedVCFpreload filterDeaminatedVCFpreload.o subsamplebamFixedNumber subsamplebamFixedNumberPair addRGinHeaderHack splitByChr filterDeaminatedVCFpreload1000g filterDeaminatedpreload1000g bamCat splitByRG tallyByRG addRG removeTagsMapping addRG_CTEAM insertSize singleAndFirstMate cutReadsDistribution  cutStart *.o 
