
CXX      = g++ 
BAMTOOLS=  bamtools/
LIBGAB   = libgab/
#MISTARTOOLS   = mistartools/

CXXFLAGS  = -lm -O3 -Wall -I${LIBGAB} -I${LIBGAB}/gzstream/   -I${BAMTOOLS}/src//  -I${BAMTOOLS}/build/src/ -c
#LDFLAGS  += -lbamtools
LDLIBS   += ${LIBGAB}/PutProgramInHeader.o ${BAMTOOLS}/build/src/libbamtools.a   -lm -lz 

#${MISTARTOOLS}/VCFreader.o filterDeaminatedVCF filterDeaminatedVCFpreload filterDeaminatedVCFpreload1000g  filterDeaminatedFasta filterDeaminatedpreload1000g

all: ${BAMTOOLS}/build/src/libbamtools.a  ${LIBGAB}/libgab.a  allFailqc allPassqc allNoDUPRM cutDeaminated decrQualNs decrQualDeaminated decrQualDeaminatedDoubleStranded failQualPair filterDeaminated  removeRG replaceRG cutUMI removeThoseWithoutMD retrieveRG subSampleBAM transBAM transBAMperRead  filterHighEditDistance.o filterHighEditDistance  editDist removeUnalignedANDWrongCigar dumpLoneMates  retrieveMapped_single_and_ProperlyPair retrieveMapped_single_and_ProperlyPair_NoPCR setAsUnpaired compareRG  subsamplebamFixedNumber subsamplebamFixedNumberPair addFFtomakeUdosDamagePatternHappy addRGinHeaderHack splitByChr splitByRG filterEditDist bamCat tallyByRG addRG removeTagsMapping addRG_CTEAM baseQualScorePerCycle insertSize singleAndFirstMate cutReadsDistribution cutStart cutEndKeepBeginning addRGInReadAndHeader removeIndices retrieveReadsWithName readWhereEitherIsMapped sumBases filterOnlyThoseWithRG errorRate errorQCScores filterDeaminatedDouble crossSampleStat retrieveMappedCertainIsize subsamplebamFixedNumberCollate detectUMIvarLength filterNoNonDNAbase

${LIBGAB}/libgab.h:
	rm -rf ${LIBGAB}
	git clone --recursive --depth 1 https://github.com/grenaud/libgab.git

${LIBGAB}/libgab.a:   ${LIBGAB}/libgab.h
	make -C ${LIBGAB}

${BAMTOOLS}/src/api/BamAlignment.h:
	rm -rf ${BAMTOOLS}/
	git clone --recursive --depth 1 https://github.com/pezmaster31/bamtools.git

${BAMTOOLS}/build/src/libbamtools.a: ${BAMTOOLS}/src/api/BamAlignment.h
	cd ${BAMTOOLS}/ && mkdir -p build/  && cd build/ && cmake .. && make && cd ../..

# ${MISTARTOOLS}/VCFreader.o: ${MISTARTOOLS}/VCFreader.cpp
# 	cd ${MISTARTOOLS}/ && make

# ${MISTARTOOLS}/VCFreader.cpp:
# 	rm -rf ${MISTARTOOLS}/
# 	git clone --recursive --depth 1 https://github.com/grenaud/mistartools.git



%.o: %.cpp
	${CXX} ${CXXFLAGS} $^ -o $@


crossSampleStat: crossSampleStat.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 


sumBases: sumBases.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 

readWhereEitherIsMapped: readWhereEitherIsMapped.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 

compareRG: compareRG.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 

addRG: addRG.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 

filterOnlyThoseWithRG: filterOnlyThoseWithRG.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 

filterNoNonDNAbase: filterNoNonDNAbase.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 

cutStart: cutStart.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 

cutEndKeepBeginning: cutEndKeepBeginning.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 

addRG_CTEAM: addRG_CTEAM.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 

removeTagsMapping: removeTagsMapping.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 

pairsAsSingle: pairsAsSingle.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 

retrieveReadsWithName: retrieveReadsWithName.o  ${LIBGAB}libgab.a ${LIBGAB}/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) 

cutReadsDistribution: cutReadsDistribution.o  ${LIBGAB}libgab.a ${LIBGAB}/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) 

filterEditDist: filterEditDist.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 


bamCat: bamCat.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 


dumpLoneMates: dumpLoneMates.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 


removeUnalignedANDWrongCigar: removeUnalignedANDWrongCigar.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 


retrieveMappedCertainIsize: retrieveMappedCertainIsize.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 

retrieveMapped_single_and_ProperlyPair: retrieveMapped_single_and_ProperlyPair.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 


retrieveMapped_single_and_ProperlyPair_NoPCR: retrieveMapped_single_and_ProperlyPair_NoPCR.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 


editDist: editDist.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 

errorRate: errorRate.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 

insertSize: insertSize.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 


allFailqc: allFailqc.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 

detectUMIvarLength: detectUMIvarLength.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 

singleAndFirstMate: singleAndFirstMate.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 


setAsUnpaired: setAsUnpaired.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 


allPassqc: allPassqc.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 

allNoDUPRM: allNoDUPRM.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 

splitByChr: splitByChr.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 

splitByRG: splitByRG.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 

tallyByRG: tallyByRG.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 

cutDeaminated: cutDeaminated.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 

decrQualNs: decrQualNs.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 

decrQualDeaminated: decrQualDeaminated.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 

decrQualDeaminatedDoubleStranded: decrQualDeaminatedDoubleStranded.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 

failQualPair.o: failQualPair.cpp
	${CXX} ${CXXFLAGS} failQualPair.cpp

failQualPair: failQualPair.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 


filterHighEditDistance: filterHighEditDistance.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 


baseQualScorePerCycle: baseQualScorePerCycle.o  ${LIBGAB}libgab.a ${LIBGAB}ReconsReferenceBAM.o
	${CXX} -o $@ $^ $(LDLIBS) 


#filterDeaminatedVCF: filterDeaminatedVCF.o  ${LIBGAB}libgab.a ${MISTARTOOLS}/VCFreader.o ${LIBGAB}/gzstream/libgzstream.a  ${MISTARTOOLS}/SimpleVCF.o ${MISTARTOOLS}/CoreVCF.o ${MISTARTOOLS}/ReadTabix.o ${LIBGAB}ReconsReferenceBAM.o ${MISTARTOOLS}/tabix/libtabix.a
#	${CXX} -o $@ $^ $(LDLIBS) 

#filterDeaminatedVCFpreload: filterDeaminatedVCFpreload.o  ${LIBGAB}libgab.a ${MISTARTOOLS}/VCFreader.o ${LIBGAB}/gzstream/libgzstream.a  ${MISTARTOOLS}/SimpleVCF.o ${MISTARTOOLS}/CoreVCF.o ${MISTARTOOLS}/ReadTabix.o ${LIBGAB}ReconsReferenceBAM.o ${MISTARTOOLS}/tabix/libtabix.a 
	${CXX} -o $@ $^ $(LDLIBS) 
#
#filterDeaminatedVCFpreload1000g: filterDeaminatedVCFpreload1000g.o  ${LIBGAB}libgab.a ${MISTARTOOLS}/VCFreader.o ${LIBGAB}/gzstream/libgzstream.a  ${MISTARTOOLS}/SimpleVCF.o ${MISTARTOOLS}/CoreVCF.o ${MISTARTOOLS}/ReadTabix.o ${LIBGAB}ReconsReferenceBAM.o ${MISTARTOOLS}/tabix/libtabix.a
	${CXX} -o $@ $^ $(LDLIBS) 

#filterDeaminatedpreload1000g: filterDeaminatedpreload1000g.o  ${LIBGAB}libgab.a ${MISTARTOOLS}/VCFreader.o ${LIBGAB}/gzstream/libgzstream.a  ${MISTARTOOLS}/SimpleVCF.o ${MISTARTOOLS}/CoreVCF.o ${MISTARTOOLS}/ReadTabix.o ${LIBGAB}ReconsReferenceBAM.o ${MISTARTOOLS}/tabix/libtabix.a
	${CXX} -o $@ $^ $(LDLIBS) 

#filterDeaminatedFasta: filterDeaminatedFasta.o  ${LIBGAB}libgab.a ${LIBGAB}ReconsReferenceBAM.o ${LIBGAB}FastQParser.o ${LIBGAB}FastQObj.o ${LIBGAB}/gzstream/libgzstream.a 
#	${CXX} -o $@ $^ $(LDLIBS) 

filterDeaminated: filterDeaminated.o  ${LIBGAB}libgab.a ${LIBGAB}ReconsReferenceBAM.o ${LIBGAB}FastQParser.o ${LIBGAB}FastQObj.o ${LIBGAB}/gzstream/libgzstream.a 
	${CXX} -o $@ $^ $(LDLIBS) 

filterDeaminatedDouble: filterDeaminatedDouble.o  ${LIBGAB}libgab.a ${LIBGAB}ReconsReferenceBAM.o ${LIBGAB}FastQParser.o ${LIBGAB}FastQObj.o ${LIBGAB}/gzstream/libgzstream.a 
	${CXX} -o $@ $^ $(LDLIBS) 

removeIndices: removeIndices.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 

removeRG: removeRG.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 

replaceRG: replaceRG.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 

cutUMI: cutUMI.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 

removeThoseWithoutMD: removeThoseWithoutMD.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 

addRGinHeaderHack: addRGinHeaderHack.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 

addFFtomakeUdosDamagePatternHappy: addFFtomakeUdosDamagePatternHappy.o ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 

addRGInReadAndHeader: addRGInReadAndHeader.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 

retrieveRG: retrieveRG.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 


subSampleBAM: subSampleBAM.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 

subsamplebamFixedNumber: subsamplebamFixedNumber.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 

subsamplebamFixedNumberPair: subsamplebamFixedNumberPair.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 

subsamplebamFixedNumberCollate: subsamplebamFixedNumberCollate.o  ${LIBGAB}libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 

transBAM: transBAM.o  ${LIBGAB}libgab.a ${LIBGAB}ReconsReferenceBAM.o
	${CXX} -o $@ $^ $(LDLIBS) 



transBAMperRead: transBAMperRead.o  ${LIBGAB}libgab.a ${LIBGAB}ReconsReferenceBAM.o
	${CXX} -o $@ $^ $(LDLIBS) 

errorQCScores: errorQCScores.o  ${LIBGAB}libgab.a ${LIBGAB}ReconsReferenceBAM.o
	${CXX} -o $@ $^ $(LDLIBS) 





clean :
	rm -f  allFailqc detectUMIvarLength allNoDUPRM allPassqc cutDeaminated decrQualNs decrQualDeaminated decrQualDeaminatedDoubleStranded failQualPair filterDeaminated filterDeaminatedVCF getCtrlReadsBAM removeRG replaceRG cutUMI removeThoseWithoutMD removeIndices retrieveRG subSampleBAM transBAM transBAMperRead ${LIBGAB}ReconsReferenceBAM.o filterHighEditDistance.o filterHighEditDistance  removeUnalignedANDWrongCigar.o removeUnalignedANDWrongCigar retrieveMapped_single_and_ProperlyPair.o retrieveMapped_single_and_ProperlyPair retrieveMapped_single_and_ProperlyPair_NoPCR setAsUnpaired filterDeaminatedFasta compareRG filterDeaminatedFasta  compareRG filterDeaminatedVCFpreload filterDeaminatedVCFpreload.o subsamplebamFixedNumber subsamplebamFixedNumberPair addRGinHeaderHack splitByChr filterDeaminatedVCFpreload1000g filterDeaminatedpreload1000g bamCat splitByRG tallyByRG addRG removeTagsMapping addRG_CTEAM insertSize singleAndFirstMate cutReadsDistribution  cutStart cutEndKeepBeginning addRGInReadAndHeader removeTagsMapping readWhereEitherIsMapped sumBases addFFtomakeUdosDamagePatternHappy filterOnlyThoseWithRG errorRate filterDeaminatedDouble crossSampleStat retrieveMappedCertainIsize subsamplebamFixedNumberPair *.o 

