
CXX      = g++ 
BAMTOOLS=  bamtools/
LIBGAB   = libgab/
#MISTARTOOLS   = mistartools/

CXXFLAGS  = -lm -O3 -Wall -I${LIBGAB} -I${LIBGAB}/gzstream/   -I${BAMTOOLS}/src// -c
#LDFLAGS  += -lbamtools
LDLIBS   += ${LIBGAB}/PutProgramInHeader.o ${BAMTOOLS}/build/src/api/libbamtools.a   -lm -lz 

#${MISTARTOOLS}/VCFreader.o filterDeaminatedVCF filterDeaminatedVCFpreload filterDeaminatedVCFpreload1000g  filterDeaminatedFasta filterDeaminatedpreload1000g
all: ${BAMTOOLS}/build/src/api/libbamtools.a  ${LIBGAB}/utils.o  allFailqc allPassqc cutDeaminated decrQualDeaminated decrQualDeaminatedDoubleStranded failQualPair filterDeaminated  removeRG retrieveRG subSampleBAM transBAM transBAMperRead  filterHighEditDistance.o filterHighEditDistance  editDist removeUnalignedANDWrongCigar dumpLoneMates  retrieveMapped_single_and_ProperlyPair retrieveMapped_single_and_ProperlyPair_NoPCR setAsUnpaired compareRG  subsamplebamFixedNumber subsamplebamFixedNumberPair addFFtomakeUdosDamagePatternHappy addRGinHeaderHack splitByChr splitByRG filterEditDist bamCat tallyByRG addRG removeTagsMapping addRG_CTEAM baseQualScorePerCycle insertSize singleAndFirstMate cutReadsDistribution cutStart cutEndKeepBeginning addRGInReadAndHeader removeIndices retrieveReadsWithName readWhereEitherIsMapped sumBases filterOnlyThoseWithRG errorRate errorQCScores filterDeaminatedDouble crossSampleStat

${LIBGAB}/utils.h:
	rm -rf ${LIBGAB}
	git clone --recursive --depth 1 https://github.com/grenaud/libgab.git

${LIBGAB}/utils.o:   ${LIBGAB}/utils.h
	make -C ${LIBGAB}

${BAMTOOLS}/src/api/BamAlignment.h:
	rm -rf ${BAMTOOLS}/
	git clone --recursive --depth 1 https://github.com/pezmaster31/bamtools.git

${BAMTOOLS}/build/src/api/libbamtools.a: ${BAMTOOLS}/src/api/BamAlignment.h
	cd ${BAMTOOLS}/ && mkdir -p build/  && cd build/ && cmake .. && make && cd ../..

# ${MISTARTOOLS}/VCFreader.o: ${MISTARTOOLS}/VCFreader.cpp
# 	cd ${MISTARTOOLS}/ && make

# ${MISTARTOOLS}/VCFreader.cpp:
# 	rm -rf ${MISTARTOOLS}/
# 	git clone --recursive --depth 1 https://github.com/grenaud/mistartools.git



%.o: %.cpp
	${CXX} ${CXXFLAGS} $^ -o $@


crossSampleStat: crossSampleStat.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 


sumBases: sumBases.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

readWhereEitherIsMapped: readWhereEitherIsMapped.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

compareRG: compareRG.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

addRG: addRG.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

filterOnlyThoseWithRG: filterOnlyThoseWithRG.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

cutStart: cutStart.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

cutEndKeepBeginning: cutEndKeepBeginning.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

addRG_CTEAM: addRG_CTEAM.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

removeTagsMapping: removeTagsMapping.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

retrieveReadsWithName: retrieveReadsWithName.o  ${LIBGAB}utils.o ${LIBGAB}/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) 

cutReadsDistribution: cutReadsDistribution.o  ${LIBGAB}utils.o ${LIBGAB}/gzstream/libgzstream.a
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


retrieveMapped_single_and_ProperlyPair_NoPCR: retrieveMapped_single_and_ProperlyPair_NoPCR.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 


editDist: editDist.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

errorRate: errorRate.o  ${LIBGAB}utils.o
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


#filterDeaminatedVCF: filterDeaminatedVCF.o  ${LIBGAB}utils.o ${MISTARTOOLS}/VCFreader.o ${LIBGAB}/gzstream/libgzstream.a  ${MISTARTOOLS}/SimpleVCF.o ${MISTARTOOLS}/CoreVCF.o ${MISTARTOOLS}/ReadTabix.o ${LIBGAB}ReconsReferenceBAM.o ${MISTARTOOLS}/tabix/libtabix.a
#	${CXX} -o $@ $^ $(LDLIBS) 

#filterDeaminatedVCFpreload: filterDeaminatedVCFpreload.o  ${LIBGAB}utils.o ${MISTARTOOLS}/VCFreader.o ${LIBGAB}/gzstream/libgzstream.a  ${MISTARTOOLS}/SimpleVCF.o ${MISTARTOOLS}/CoreVCF.o ${MISTARTOOLS}/ReadTabix.o ${LIBGAB}ReconsReferenceBAM.o ${MISTARTOOLS}/tabix/libtabix.a 
	${CXX} -o $@ $^ $(LDLIBS) 
#
#filterDeaminatedVCFpreload1000g: filterDeaminatedVCFpreload1000g.o  ${LIBGAB}utils.o ${MISTARTOOLS}/VCFreader.o ${LIBGAB}/gzstream/libgzstream.a  ${MISTARTOOLS}/SimpleVCF.o ${MISTARTOOLS}/CoreVCF.o ${MISTARTOOLS}/ReadTabix.o ${LIBGAB}ReconsReferenceBAM.o ${MISTARTOOLS}/tabix/libtabix.a
	${CXX} -o $@ $^ $(LDLIBS) 

#filterDeaminatedpreload1000g: filterDeaminatedpreload1000g.o  ${LIBGAB}utils.o ${MISTARTOOLS}/VCFreader.o ${LIBGAB}/gzstream/libgzstream.a  ${MISTARTOOLS}/SimpleVCF.o ${MISTARTOOLS}/CoreVCF.o ${MISTARTOOLS}/ReadTabix.o ${LIBGAB}ReconsReferenceBAM.o ${MISTARTOOLS}/tabix/libtabix.a
	${CXX} -o $@ $^ $(LDLIBS) 

#filterDeaminatedFasta: filterDeaminatedFasta.o  ${LIBGAB}utils.o ${LIBGAB}ReconsReferenceBAM.o ${LIBGAB}FastQParser.o ${LIBGAB}FastQObj.o ${LIBGAB}/gzstream/libgzstream.a 
#	${CXX} -o $@ $^ $(LDLIBS) 

filterDeaminated: filterDeaminated.o  ${LIBGAB}utils.o ${LIBGAB}ReconsReferenceBAM.o ${LIBGAB}FastQParser.o ${LIBGAB}FastQObj.o ${LIBGAB}/gzstream/libgzstream.a 
	${CXX} -o $@ $^ $(LDLIBS) 

filterDeaminatedDouble: filterDeaminatedDouble.o  ${LIBGAB}utils.o ${LIBGAB}ReconsReferenceBAM.o ${LIBGAB}FastQParser.o ${LIBGAB}FastQObj.o ${LIBGAB}/gzstream/libgzstream.a 
	${CXX} -o $@ $^ $(LDLIBS) 

removeIndices: removeIndices.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

removeRG: removeRG.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

addRGinHeaderHack: addRGinHeaderHack.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

addFFtomakeUdosDamagePatternHappy: addFFtomakeUdosDamagePatternHappy.o ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

addRGInReadAndHeader: addRGInReadAndHeader.o  ${LIBGAB}utils.o
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

errorQCScores: errorQCScores.o  ${LIBGAB}utils.o ${LIBGAB}ReconsReferenceBAM.o
	${CXX} -o $@ $^ $(LDLIBS) 





clean :
	rm -f  allFailqc allPassqc cutDeaminated decrQualDeaminated decrQualDeaminatedDoubleStranded failQualPair filterDeaminated filterDeaminatedVCF getCtrlReadsBAM removeRG removeIndices retrieveRG subSampleBAM transBAM transBAMperRead ${LIBGAB}ReconsReferenceBAM.o filterHighEditDistance.o filterHighEditDistance  removeUnalignedANDWrongCigar.o removeUnalignedANDWrongCigar retrieveMapped_single_and_ProperlyPair.o retrieveMapped_single_and_ProperlyPair retrieveMapped_single_and_ProperlyPair_NoPCR setAsUnpaired filterDeaminatedFasta compareRG filterDeaminatedFasta  compareRG filterDeaminatedVCFpreload filterDeaminatedVCFpreload.o subsamplebamFixedNumber subsamplebamFixedNumberPair addRGinHeaderHack splitByChr filterDeaminatedVCFpreload1000g filterDeaminatedpreload1000g bamCat splitByRG tallyByRG addRG removeTagsMapping addRG_CTEAM insertSize singleAndFirstMate cutReadsDistribution  cutStart cutEndKeepBeginning addRGInReadAndHeader removeTagsMapping readWhereEitherIsMapped sumBases addFFtomakeUdosDamagePatternHappy filterOnlyThoseWithRG errorRate filterDeaminatedDouble crossSampleStat *.o 
