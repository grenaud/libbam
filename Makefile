
CXX      = g++ 
BAMTOOLS=  bamtools/
LIBGAB   = libgab/
MISTARTOOLS   = mistartools/

CXXFLAGS  = -lm -O3 -Wall -I${LIBGAB} -I${MISTARTOOLS} -I${MISTARTOOLS}/tabix/ -I${LIBGAB}/gzstream/   -I${BAMTOOLS}/include/ -c
#LDFLAGS  += -lbamtools
LDLIBS   += ${LIBGAB}/PutProgramInHeader.o ${BAMTOOLS}/lib/libbamtools.so   -lm -lz 


all: ${LIBGAB}utils.o ${MISTARTOOLS}/VCFreader.o allFailqc allPassqc cutDeaminated decrQualDeaminated decrQualDeaminatedDoubleStranded failQualPair filterDeaminated filterDeaminatedVCF removeRG retrieveRG subSampleBAM transBAM transBAMperRead  filterHighEditDistance.o filterHighEditDistance  editDist removeUnalignedANDWrongCigar dumpLoneMates retrieveMapped_single_and_ProperlyPair setAsUnpaired compareRG  filterDeaminatedFasta filterDeaminatedVCFpreload filterDeaminatedVCFpreload1000g filterDeaminatedpreload1000g subsamplebamFixedNumber subsamplebamFixedNumberPair addFFtomakeUdosDamagePatternHappy addRGinHeaderHack splitByChr splitByRG filterEditDist bamCat tallyByRG addRG removeTagsMapping addRG_CTEAM baseQualScorePerCycle insertSize singleAndFirstMate cutReadsDistribution cutStart cutEndKeepBeginning addRGInReadAndHeader removeIndices retrieveReadsWithName readWhereEitherIsMapped sumBases filterOnlyThoseWithRG errorRate errorQCScores filterDeaminatedDouble crossSampleStat

${LIBGAB}/utils.h:
	rm -rf ${LIBGAB}
	git clone --recursive --depth 1 https://github.com/grenaud/libgab.git


${LIBGAB}utils.o: ${BAMTOOLS}/lib/libbamtools.so  ${LIBGAB}/utils.h
	make -C ${LIBGAB}

${BAMTOOLS}/src/api/BamAlignment.h:
	rm -rf ${BAMTOOLS}/
	git clone --recursive --depth 1 https://github.com/pezmaster31/bamtools.git

${BAMTOOLS}/lib/libbamtools.so: ${BAMTOOLS}/src/api/BamAlignment.h
	cd ${BAMTOOLS}/ && mkdir -p build/  && cd build/ && cmake .. && make && cd ../..

${MISTARTOOLS}/VCFreader.o: ${MISTARTOOLS}/VCFreader.cpp
	cd ${MISTARTOOLS}/ && make

${MISTARTOOLS}/VCFreader.cpp:
	rm -rf ${MISTARTOOLS}/
	git clone --recursive --depth 1 https://github.com/grenaud/mistartools.git



%.o: %.cpp
	${CXX} ${CXXFLAGS} $^ -o $@


crossSampleStat: ${BAMTOOLS}/lib/libbamtools.so  crossSampleStat.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 


sumBases: ${BAMTOOLS}/lib/libbamtools.so  sumBases.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

readWhereEitherIsMapped: ${BAMTOOLS}/lib/libbamtools.so  readWhereEitherIsMapped.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

compareRG: ${BAMTOOLS}/lib/libbamtools.so  compareRG.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

addRG: ${BAMTOOLS}/lib/libbamtools.so  addRG.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

filterOnlyThoseWithRG: ${BAMTOOLS}/lib/libbamtools.so  filterOnlyThoseWithRG.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

cutStart: ${BAMTOOLS}/lib/libbamtools.so  cutStart.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

cutEndKeepBeginning: ${BAMTOOLS}/lib/libbamtools.so  cutEndKeepBeginning.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

addRG_CTEAM: ${BAMTOOLS}/lib/libbamtools.so  addRG_CTEAM.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

removeTagsMapping: ${BAMTOOLS}/lib/libbamtools.so  removeTagsMapping.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

retrieveReadsWithName: ${BAMTOOLS}/lib/libbamtools.so  retrieveReadsWithName.o  ${LIBGAB}utils.o ${LIBGAB}/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) 

cutReadsDistribution: ${BAMTOOLS}/lib/libbamtools.so  cutReadsDistribution.o  ${LIBGAB}utils.o ${LIBGAB}/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) 

filterEditDist: ${BAMTOOLS}/lib/libbamtools.so  filterEditDist.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 


bamCat: ${BAMTOOLS}/lib/libbamtools.so  bamCat.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 


dumpLoneMates: ${BAMTOOLS}/lib/libbamtools.so  dumpLoneMates.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 


removeUnalignedANDWrongCigar: ${BAMTOOLS}/lib/libbamtools.so  removeUnalignedANDWrongCigar.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 


retrieveMapped_single_and_ProperlyPair: ${BAMTOOLS}/lib/libbamtools.so  retrieveMapped_single_and_ProperlyPair.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 


editDist: ${BAMTOOLS}/lib/libbamtools.so  editDist.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

errorRate: ${BAMTOOLS}/lib/libbamtools.so  errorRate.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

insertSize: ${BAMTOOLS}/lib/libbamtools.so  insertSize.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 


allFailqc: ${BAMTOOLS}/lib/libbamtools.so  allFailqc.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

singleAndFirstMate: ${BAMTOOLS}/lib/libbamtools.so  singleAndFirstMate.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 


setAsUnpaired: ${BAMTOOLS}/lib/libbamtools.so  setAsUnpaired.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 



allPassqc: ${BAMTOOLS}/lib/libbamtools.so  allPassqc.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 


splitByChr: ${BAMTOOLS}/lib/libbamtools.so  splitByChr.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

splitByRG: ${BAMTOOLS}/lib/libbamtools.so  splitByRG.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

tallyByRG: ${BAMTOOLS}/lib/libbamtools.so  tallyByRG.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

cutDeaminated: ${BAMTOOLS}/lib/libbamtools.so  cutDeaminated.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 


decrQualDeaminated: ${BAMTOOLS}/lib/libbamtools.so  decrQualDeaminated.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 


decrQualDeaminatedDoubleStranded: ${BAMTOOLS}/lib/libbamtools.so  decrQualDeaminatedDoubleStranded.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

failQualPair.o: ${BAMTOOLS}/lib/libbamtools.so  failQualPair.cpp
	${CXX} ${CXXFLAGS} failQualPair.cpp

failQualPair: ${BAMTOOLS}/lib/libbamtools.so  failQualPair.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 


filterHighEditDistance: ${BAMTOOLS}/lib/libbamtools.so  filterHighEditDistance.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 


baseQualScorePerCycle: ${BAMTOOLS}/lib/libbamtools.so  baseQualScorePerCycle.o  ${LIBGAB}utils.o ${LIBGAB}ReconsReferenceBAM.o
	${CXX} -o $@ $^ $(LDLIBS) 


filterDeaminatedVCF: ${BAMTOOLS}/lib/libbamtools.so  filterDeaminatedVCF.o  ${LIBGAB}utils.o ${MISTARTOOLS}/VCFreader.o ${LIBGAB}/gzstream/libgzstream.a  ${MISTARTOOLS}/SimpleVCF.o ${MISTARTOOLS}/CoreVCF.o ${MISTARTOOLS}/ReadTabix.o ${LIBGAB}ReconsReferenceBAM.o ${MISTARTOOLS}/tabix/libtabix.a
	${CXX} -o $@ $^ $(LDLIBS) 

filterDeaminatedVCFpreload: ${BAMTOOLS}/lib/libbamtools.so  filterDeaminatedVCFpreload.o  ${LIBGAB}utils.o ${MISTARTOOLS}/VCFreader.o ${LIBGAB}/gzstream/libgzstream.a  ${MISTARTOOLS}/SimpleVCF.o ${MISTARTOOLS}/CoreVCF.o ${MISTARTOOLS}/ReadTabix.o ${LIBGAB}ReconsReferenceBAM.o ${MISTARTOOLS}/tabix/libtabix.a 
	${CXX} -o $@ $^ $(LDLIBS) 

filterDeaminatedVCFpreload1000g: ${BAMTOOLS}/lib/libbamtools.so  filterDeaminatedVCFpreload1000g.o  ${LIBGAB}utils.o ${MISTARTOOLS}/VCFreader.o ${LIBGAB}/gzstream/libgzstream.a  ${MISTARTOOLS}/SimpleVCF.o ${MISTARTOOLS}/CoreVCF.o ${MISTARTOOLS}/ReadTabix.o ${LIBGAB}ReconsReferenceBAM.o ${MISTARTOOLS}/tabix/libtabix.a
	${CXX} -o $@ $^ $(LDLIBS) 

filterDeaminatedpreload1000g: ${BAMTOOLS}/lib/libbamtools.so  filterDeaminatedpreload1000g.o  ${LIBGAB}utils.o ${MISTARTOOLS}/VCFreader.o ${LIBGAB}/gzstream/libgzstream.a  ${MISTARTOOLS}/SimpleVCF.o ${MISTARTOOLS}/CoreVCF.o ${MISTARTOOLS}/ReadTabix.o ${LIBGAB}ReconsReferenceBAM.o ${MISTARTOOLS}/tabix/libtabix.a
	${CXX} -o $@ $^ $(LDLIBS) 

filterDeaminatedFasta: ${BAMTOOLS}/lib/libbamtools.so  filterDeaminatedFasta.o  ${LIBGAB}utils.o ${LIBGAB}ReconsReferenceBAM.o ${LIBGAB}FastQParser.o ${LIBGAB}FastQObj.o ${LIBGAB}/gzstream/libgzstream.a 
	${CXX} -o $@ $^ $(LDLIBS) 

filterDeaminated: ${BAMTOOLS}/lib/libbamtools.so  filterDeaminated.o  ${LIBGAB}utils.o ${LIBGAB}ReconsReferenceBAM.o ${LIBGAB}FastQParser.o ${LIBGAB}FastQObj.o ${LIBGAB}/gzstream/libgzstream.a 
	${CXX} -o $@ $^ $(LDLIBS) 

filterDeaminatedDouble: ${BAMTOOLS}/lib/libbamtools.so  filterDeaminatedDouble.o  ${LIBGAB}utils.o ${LIBGAB}ReconsReferenceBAM.o ${LIBGAB}FastQParser.o ${LIBGAB}FastQObj.o ${LIBGAB}/gzstream/libgzstream.a 
	${CXX} -o $@ $^ $(LDLIBS) 

removeIndices: ${BAMTOOLS}/lib/libbamtools.so  removeIndices.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

removeRG: ${BAMTOOLS}/lib/libbamtools.so  removeRG.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

addRGinHeaderHack: ${BAMTOOLS}/lib/libbamtools.so  addRGinHeaderHack.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

addFFtomakeUdosDamagePatternHappy: ${BAMTOOLS}/lib/libbamtools.so  addFFtomakeUdosDamagePatternHappy.o ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

addRGInReadAndHeader: ${BAMTOOLS}/lib/libbamtools.so  addRGInReadAndHeader.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

retrieveRG: ${BAMTOOLS}/lib/libbamtools.so  retrieveRG.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 


subSampleBAM: ${BAMTOOLS}/lib/libbamtools.so  subSampleBAM.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

subsamplebamFixedNumber: ${BAMTOOLS}/lib/libbamtools.so  subsamplebamFixedNumber.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

subsamplebamFixedNumberPair: ${BAMTOOLS}/lib/libbamtools.so  subsamplebamFixedNumberPair.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

transBAM: ${BAMTOOLS}/lib/libbamtools.so  transBAM.o  ${LIBGAB}utils.o ${LIBGAB}ReconsReferenceBAM.o
	${CXX} -o $@ $^ $(LDLIBS) 



transBAMperRead: ${BAMTOOLS}/lib/libbamtools.so  transBAMperRead.o  ${LIBGAB}utils.o ${LIBGAB}ReconsReferenceBAM.o
	${CXX} -o $@ $^ $(LDLIBS) 

errorQCScores: ${BAMTOOLS}/lib/libbamtools.so  errorQCScores.o  ${LIBGAB}utils.o ${LIBGAB}ReconsReferenceBAM.o
	${CXX} -o $@ $^ $(LDLIBS) 





clean :
	rm -f  allFailqc allPassqc cutDeaminated decrQualDeaminated decrQualDeaminatedDoubleStranded failQualPair filterDeaminated filterDeaminatedVCF getCtrlReadsBAM removeRG removeIndices retrieveRG subSampleBAM transBAM transBAMperRead ${LIBGAB}ReconsReferenceBAM.o filterHighEditDistance.o filterHighEditDistance  removeUnalignedANDWrongCigar.o removeUnalignedANDWrongCigar retrieveMapped_single_and_ProperlyPair.o retrieveMapped_single_and_ProperlyPair setAsUnpaired filterDeaminatedFasta compareRG filterDeaminatedFasta  compareRG filterDeaminatedVCFpreload filterDeaminatedVCFpreload.o subsamplebamFixedNumber subsamplebamFixedNumberPair addRGinHeaderHack splitByChr filterDeaminatedVCFpreload1000g filterDeaminatedpreload1000g bamCat splitByRG tallyByRG addRG removeTagsMapping addRG_CTEAM insertSize singleAndFirstMate cutReadsDistribution  cutStart cutEndKeepBeginning addRGInReadAndHeader removeTagsMapping readWhereEitherIsMapped sumBases addFFtomakeUdosDamagePatternHappy filterOnlyThoseWithRG errorRate filterDeaminatedDouble crossSampleStat *.o 
