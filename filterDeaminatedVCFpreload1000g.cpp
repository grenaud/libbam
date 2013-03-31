#include <iostream>
#include <vector>
#include <set>
#include <ctype.h>
#include <stdlib.h>

#include "api/BamMultiReader.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"
#include "api/SamSequenceDictionary.h"

#include "utils.h"
#include "ReconsReferenceBAM.h"
#include "VCFreader.h"

using namespace std;
using namespace BamTools;

const int offset=33;
string vcf1000g="/mnt/sequencedb/1000Genomes/ftp/phase1/analysis_results/integrated_call_sets/ALL.wgs.integrated_phase1_v3.20101123.snps_indels_sv.sites.vcf.gz";

inline bool isTransition(char reference,char read){

    if(  (reference == 'A' && read == 'G') ||
	 (reference == 'G' && read == 'A') ||
	 (reference == 'C' && read == 'T') ||
	 (reference == 'T' && read == 'C') ){
	return true;
    }
	 
    return false;
}



inline bool isTransversions(char reference,char read){
    if(  (reference == 'A' && read == 'C') ||
	 (reference == 'A' && read == 'T') ||
	 
	 (reference == 'G' && read == 'C') ||
	 (reference == 'G' && read == 'T') ||
	 
	 (reference == 'C' && read == 'A') ||
	 (reference == 'C' && read == 'G') ||

	 (reference == 'T' && read == 'A') ||
	 (reference == 'T' && read == 'G') ){
	return true;
    }

    return false;	
}


inline bool hasIinfirstOrLastTwoBases(const string & reconstructedReference){
    if(reconstructedReference.length() <= 4){
	cerr<<"ERROR read has length less than 4 bp"<<endl;
	exit(1);
    }

    for(unsigned int j=0;j<2;j++){
	if(reconstructedReference[j] == 'I')
	    return true;
    }


    for(unsigned int j=(reconstructedReference.length()-2);
	j<(reconstructedReference.length());
	j++){
	if(reconstructedReference[j] == 'I')
	    return true;
    }

    return false;
}

inline bool deletionsNextToTwo(const BamAlignment  * al){
    vector<int> lengthOfNonDels;
    vector<CigarOp> cigarData=al->CigarData;
    bool foundDel=false;
    for(unsigned int i=0;i<cigarData.size();i++){
        if(cigarData[i].Type == 'D'){
	    foundDel=true;
	}else{
	    lengthOfNonDels.push_back(cigarData[i].Length);
	}

	//toReturn+=cigarData[i].Length;      
	//reconstructedTemp+=string(cigarData[i].Length,cigarData[i].Type);
    }

    if(foundDel){
	if(lengthOfNonDels[0]<=2)
	    return true;
	if(lengthOfNonDels[ lengthOfNonDels.size() -1 ]<=2)
	    return true;
    }
    
    return false;
}

 

inline int countMatchesRecons(const string & reconstructedReference,int minusIndex){
    int lengthMatches=0;
	

    for(unsigned int j=0;j<(reconstructedReference.length()-minusIndex);j++){
	if(reconstructedReference[j] == 'M' ||
	   reconstructedReference[j] == 'N' ||
	   //reconstructedReference[j] == 'I' ||
	   reconstructedReference[j] == 'A' ||
	   reconstructedReference[j] == 'C' ||
	   reconstructedReference[j] == 'G' ||
	   reconstructedReference[j] == 'T' ){
	    lengthMatches++;
	}
    }

    return lengthMatches;
}


inline void transformRef(char * refeBase,char * readBase){
    if( (*refeBase) == 'M'){
	(*refeBase)=(*readBase);
    }
    
}


//checks for an 'R' or 'S' for soft clip
inline bool hasBadCharacter(const string & reconstructedReference){

    for(unsigned int j=0;j<(reconstructedReference.length());j++){
	if(reconstructedReference[j] == 'R'  || 
	   reconstructedReference[j] == 'S' ){
	    return true;
	}
    }
    return false;
}

inline bool skipAlign(const string & reconstructedReference,const BamAlignment  * al,unsigned int * skipped){
    if(hasBadCharacter(reconstructedReference)){
	(*skipped)++;
	return true;
    }
	

    if(hasIinfirstOrLastTwoBases(reconstructedReference)){
	(*skipped)++;
	return true;
    }

    if(deletionsNextToTwo(al)){
	(*skipped)++;
	return true;
    }

    return false;
}


int main (int argc, char *argv[]) {

    int  minBaseQuality = 0;

    string usage=string(""+string(argv[0])+"  [in BAM file] [in VCF file] [chr name] [deam out BAM] [not deam out BAM]"+
			"\nThis program divides aligned single end reads into potentially deaminated\n"+
			"\nreads and the puts the rest into another bam file if the deaminated positions are not called as the alternative base in the VCF.\n"+
			"\nThis is like filterDeaminatedVCF but it loads the VCF before then labels the reads instead of doing it on the fly\n"+
			"\nwhich is good if you have many reads in the bam file.\n"+			
			"\nTip: if you do not need one of them, use /dev/null as your output\n"+
			"\narguments:\n"+
			"\t"+"--bq    [base qual]  : Minimum base quality to flag a deaminated site (Default: "+stringify(minBaseQuality)+")\n"+
			"\t"+"--1000g [vcf file]   : VCF file from 1000g to get the putative A and T positions in modern humans (Default: "+vcf1000g+")\n"+
			"\n");

    if(argc == 1 ||
       argc < 4  ||
       (argc == 2 && (string(argv[0]) == "-h" || string(argv[0]) == "--help") )
       ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }


    for(int i=1;i<(argc-2);i++){ 

	
        if(string(argv[i]) == "--bq"){
	    minBaseQuality=destringify<int>(argv[i+1]);
            i++;
            continue;
	}

        if(string(argv[i]) == "--1000g"){
	    vcf1000g=string(argv[i+1]);
            i++;
            continue;
	}



    }

    unsigned int maxSizeChromosome=250000000;//larger than chr1 hg19
    bool * hasCnoT;
    bool * hasGnoA;
    bool * thousandGenomesHasA;
    bool * thousandGenomesHasT;

    cerr<<"Trying to allocating memory"<<endl;
    try{
	hasCnoT              = new bool[ maxSizeChromosome ];
	hasGnoA              = new bool[ maxSizeChromosome ];
	thousandGenomesHasA  = new bool[ maxSizeChromosome ];
	thousandGenomesHasT  = new bool[ maxSizeChromosome ];
    }catch(bad_alloc& exc){
	cerr<<"ERROR: allocating memory failed"<<endl;
	return 1;
    }
    cerr<<"Success in allocating memory"<<endl;

    for(unsigned int i = 0;i<maxSizeChromosome;i++){
	hasCnoT[i]=false;
	hasGnoA[i]=false;
	thousandGenomesHasA[i]=false;
	thousandGenomesHasT[i]=false;
    }

    string bamfiletopen = string( argv[ argc-5 ] );
    string vcffiletopen = string( argv[ argc-4 ] );
    string chrname      = string( argv[ argc-3 ] );
    string deambam      = string( argv[ argc-2 ] );
    string nondeambam   = string( argv[ argc-1 ] );

    cerr<<"Reading consensus VCF "<<vcffiletopen<<"  ... "<<endl;

    VCFreader vcfr (vcffiletopen,
 		    // vcffiletopen+".tbi",
 		    // chrname,
 		    // 1,
 		    // maxSizeChromosome,
 		    0);

    
    while(vcfr.hasData()){
	SimpleVCF * toprint=vcfr.getData();

	if(toprint->getRef().length() != 1 )
	    continue;
	

	//if the VCF has a at least one G but no A
	if(  toprint->hasAtLeastOneG() && 
	     !toprint->hasAtLeastOneA() ){
	    hasGnoA[ toprint->getPosition() ] =true;
	}

	if(  toprint->hasAtLeastOneC() && 
	     !toprint->hasAtLeastOneT() ){
	    hasCnoT[ toprint->getPosition() ] =true;    
	}

	   
    }
    cerr<<"done reading VCF"<<endl;







    cerr<<"Reading 1000g VCF :"<<vcf1000g<<"  ..."<<endl;
    string line1000g;
    ifstream myFile1000g;
    
    myFile1000g.open(vcf1000g.c_str(), ios::in);

    if (myFile1000g.is_open()){
	while ( getline (myFile1000g,line1000g)){
	    vector<string> fields=allTokens(line1000g,'\t');
	    //0 chr
	    //1 pos
	    //2 id
	    //3 ref
	    //4 alt

	    //check if same chr
	    if(fields[0] != chrname){
		cerr <<"Error, wrong chromosome in 1000g file for line=  "<<line1000g<<endl;
		return 1;
	    }

	    //skip indels
	    if(fields[3].size() != 1 ||
	       fields[4].size() != 1 )
		continue;

	    char         ref=toupper(fields[3][0]);
	    char         alt=toupper(fields[4][0]);
	    unsigned int pos=destringify<unsigned int>( fields[1] );
	    
	    thousandGenomesHasA[ pos ] = ( (ref=='A') || (alt=='A') );
	    thousandGenomesHasT[ pos ] = ( (ref=='T') || (alt=='T') );
	    
	}
	myFile1000g.close();
    }else{
	cerr <<"Unable to open file "<<vcf1000g<<endl;
	return 1;
    }



    cerr<<"done reading 1000g VCF"<<endl;











    BamReader reader;
    
    if ( !reader.Open(bamfiletopen) ) {
    	cerr << "Could not open input BAM file"<< bamfiletopen << endl;
    	return 1;
    }


    //positioning the bam file
    int refid=reader.GetReferenceID(chrname);
    if(refid < 0){
	cerr << "Cannot retrieve the reference ID for "<< chrname << endl;
	return 1;
    }
    //cout<<"redif "<<refid<<endl;	    

    //setting the BAM reader at that position
    reader.SetRegion(refid,
		     0,
		     refid,
		     -1); 	



    vector<RefData>  testRefData=reader.GetReferenceData();
    const SamHeader header = reader.GetHeader();
    const RefVector references = reader.GetReferenceData();

    BamWriter writerDeam;
    if ( !writerDeam.Open(deambam,      header, references) ) {
	cerr << "Could not open output BAM file" << endl;
	return 1;
    }

    BamWriter writerNoDeam;
    if ( !writerNoDeam.Open(nondeambam, header, references) ) {
	cerr << "Could not open output BAM file" << endl;
	return 1;
    }



    unsigned int totalReads      =0;
    unsigned int deaminatedReads =0;
    unsigned int ndeaminatedReads =0;
    unsigned int skipped      =0;



    //iterating over the alignments for these regions
    BamAlignment al;
    int i;

    while ( reader.GetNextAlignment(al) ) {
	// cerr<<al.Name<<endl;

	//skip unmapped
	if(!al.IsMapped()){
	    skipped++;
	    continue;
	}

	//skip paired end !
	if(al.IsPaired() ){  
	    continue;
	    // cerr<<"Paired end not yet coded"<<endl;
	    // return 1;
	}


	string reconstructedReference = reconstructRef(&al);



	char refeBase;
	char readBase;
	bool isDeaminated;
	if(al.Qualities.size() != reconstructedReference.size()){
	    cerr<<"Quality line is not the same size as the reconstructed reference"<<endl;
	    return 1;
	}

	isDeaminated=false;

	if(al.IsReverseStrand()){

	    //first base next to 3'
	    i = 0 ;
	    refeBase=toupper(reconstructedReference[i]);
	    readBase=toupper(         al.QueryBases[i]);

	    if(  readBase  == 'A' && int(al.Qualities[i]-offset) >= minBaseQuality){  //isDeaminated=true; }

		if( skipAlign(reconstructedReference,&al,&skipped) ){ continue; }
		if( hasGnoA[al.Position+1] && 
		    !thousandGenomesHasA[al.Position+1] )
		    isDeaminated=true; 


		// transformRef(&refeBase,&readBase);

		// vcfr.repositionIterator(chrname,al.Position+1,al.Position+1);
		// while(vcfr.hasData()){
		//     SimpleVCF * toprint=vcfr.getData();
		//     // cout<<*toprint<<endl;
		//     //skip deletions in the alt
		//     if(toprint->getRef().length() != 1 )
		// 	continue;

		//     if(toprint->getRef()[0] != refeBase){
		// 	cerr<<reconstructedReference<<endl;
		// 	cerr<<al.Position<<endl;			
		// 	cerr<<numberOfDeletions(&al)<<endl;			
		// 	cerr<<"Problem1 position "<<*toprint<<" does not have a "<<refeBase<<" as reference allele for read "<<al.Name<<endl;
		// 	exit(1);
		//     }
		    
		//     //if the VCF has a at least one G but no A
		//     if(  toprint->hasAtLeastOneG() && 
		// 	!toprint->hasAtLeastOneA() ){
		// 	isDeaminated=true; 
		//     }
		// }

	    }


	    //second base next to 3'
	    i = 1;
	    refeBase=toupper(reconstructedReference[i]);
	    readBase=toupper(         al.QueryBases[i]);

	    //refeBase  == 'G'  &&
	    if( readBase  == 'A' &&  int(al.Qualities[i]-offset) >= minBaseQuality){  //isDeaminated=true; }

		if( skipAlign(reconstructedReference,&al,&skipped) ){ continue; }
		if( hasGnoA[al.Position+2] && 
		   !thousandGenomesHasA[al.Position+2] )
		    isDeaminated=true; 


		// transformRef(&refeBase,&readBase);


		// vcfr.repositionIterator(chrname,al.Position+2,al.Position+2);

		// while(vcfr.hasData()){
		//     SimpleVCF * toprint=vcfr.getData();
		//     // cout<<*toprint<<endl;
		//     //skip deletions in the alt
		//     if(toprint->getRef().length() != 1 )
		// 	continue;

		//     if(toprint->getRef()[0] != refeBase){
		// 	cerr<<reconstructedReference<<endl;
		// 	cerr<<al.Position<<endl;
		// 	cerr<<numberOfDeletions(&al)<<endl;
		// 	cerr<<"Problem2 position "<<*toprint<<" does not have a  "<<refeBase<<" as reference allele for read "<<al.Name<<endl;
		// 	exit(1);
		//     }

		//     //if the VCF has at least one G but no A 
		//     // if(toprint->hasAtLeastOneG() &&
		//     //    toprint->getAlt().find("A") == string::npos){
		//     if(  toprint->hasAtLeastOneG() && 
		// 	!toprint->hasAtLeastOneA() ){
		// 	isDeaminated=true; 
		//     }
		// }
	    }

	    //last  base next to 5'
	    i = (al.QueryBases.length()-1) ;
	    refeBase=toupper(reconstructedReference[i]);
	    readBase=toupper(         al.QueryBases[i]);
	    //refeBase  == 'G'  &&
	    if( readBase  == 'A' &&  int(al.Qualities[i]-offset) >= minBaseQuality){  //isDeaminated=true; }

		if( skipAlign(reconstructedReference,&al,&skipped) ){ continue; }
		transformRef(&refeBase,&readBase);


		int lengthMatches=countMatchesRecons(reconstructedReference,0);		
		int positionJump = al.Position+lengthMatches+numberOfDeletions(&al);
		if(hasGnoA[positionJump] && 
		   !thousandGenomesHasA[positionJump] )
		    isDeaminated=true; 


		// vcfr.repositionIterator(chrname,positionJump,positionJump);
		// while(vcfr.hasData()){
		//     SimpleVCF * toprint=vcfr.getData();

		//     //skip deletions in the alt
		//     if(toprint->getRef().length() != 1 )
		// 	continue;
		    
		//     if(toprint->getRef()[0] != refeBase){
		// 	cerr<<reconstructedReference<<endl;
		// 	cerr<<al.Position<<endl;
		// 	cerr<<lengthMatches<<endl;
		// 	cerr<<numberOfDeletions(&al)<<endl;
		// 	cerr<<positionJump<<endl;
		// 	cerr<<"Problem3 position "<<*toprint<<" does not have a  "<<refeBase<<" as reference allele for read "<<al.Name<<endl;
		// 	exit(1);
		//     }


		//     //if the VCF has at least one G but no A
		//     if(  toprint->hasAtLeastOneG() && 
		// 	!toprint->hasAtLeastOneA() ){
		// 	isDeaminated=true; 
		//     }
		// }

	    }

	}else{

		
	    //first base next to 5'
	    i = 0;
	    refeBase=toupper(reconstructedReference[i]);
	    readBase=toupper(         al.QueryBases[i]);
	    //refeBase  == 'C' 
	    if( readBase  == 'T' &&  int(al.Qualities[i]-offset) >= minBaseQuality){ 

		if( skipAlign(reconstructedReference,&al,&skipped) ){ continue; }
		// transformRef(&refeBase,&readBase);

		if(hasCnoT[al.Position+1]  && 
		   !thousandGenomesHasT[al.Position+1] )
		    isDeaminated=true; 

		// vcfr.repositionIterator(chrname,al.Position+1,al.Position+1);
		// while(vcfr.hasData()){
		//     SimpleVCF * toprint=vcfr.getData();
		//     //cout<<*toprint<<endl;
		//     //skip deletions in the alt
		//     if(toprint->getRef().length() != 1 )
		// 	continue;
		    
		//     if(toprint->getRef()[0] != refeBase){
		// 	cerr<<reconstructedReference<<endl;
		// 	cerr<<al.Position<<endl;
		// 	cerr<<numberOfDeletions(&al)<<endl;			
		// 	cerr<<"Problem4 position "<<*toprint<<" does not have a  "<<refeBase<<" as reference allele for read "<<al.Name<<endl;
		// 	exit(1);
		//     }

		//     //if the VCF has at least one C but no T
		//     if(  toprint->hasAtLeastOneC() && 
		// 	!toprint->hasAtLeastOneT() ){
		// 	isDeaminated=true; 
		//     }

		// }

		//cout<<al.Position+
		 
	    }

	    //second last base next to 3'
	    i = (al.QueryBases.length()-2);
	    refeBase=toupper(reconstructedReference[i]);
	    readBase=toupper(         al.QueryBases[i]);
	    //refeBase  == 'C'  &&
	    if( readBase  == 'T' &&  int(al.Qualities[i]-offset) >= minBaseQuality){  



		if( skipAlign(reconstructedReference,&al,&skipped) ){ continue; }
		//transformRef(&refeBase,&readBase);		
		int lengthMatches=countMatchesRecons(reconstructedReference,1);	
		int positionJump = al.Position+lengthMatches+numberOfDeletions(&al);
		if(hasCnoT[positionJump] && 
		   !thousandGenomesHasT[positionJump] )
		    isDeaminated=true; 

		// vcfr.repositionIterator(chrname,positionJump,positionJump);
		// while(vcfr.hasData()){
		//     SimpleVCF * toprint=vcfr.getData();
		//     //skip deletions in the alt
		//     if(toprint->getRef().length() != 1 )
		// 	continue;

		//     if(toprint->getRef()[0] != refeBase){
		// 	cerr<<reconstructedReference<<endl;
		// 	cerr<<al.Position<<endl;
		// 	cerr<<lengthMatches<<endl;
		// 	cerr<<numberOfDeletions(&al)<<endl;
		// 	cerr<<positionJump<<endl;
		// 	cerr<<"Problem5 position "<<*toprint<<" does not have a  "<<refeBase<<" as reference allele for read "<<al.Name<<endl;
		// 	exit(1);
		//     }

		//     if(  toprint->hasAtLeastOneC() && 
		// 	!toprint->hasAtLeastOneT() ){
		// 	isDeaminated=true; 
		//     }
		// }

		 

	    }

	    //last base next to 3'
	    i = (al.QueryBases.length()-1);
	    refeBase=toupper(reconstructedReference[i]);
	    readBase=toupper(         al.QueryBases[i]);
	    //&& refeBase  == 'C' 
	    if( readBase  == 'T'  && int(al.Qualities[i]-offset) >= minBaseQuality){  
		if( skipAlign(reconstructedReference,&al,&skipped) ){ continue; }
		transformRef(&refeBase,&readBase);		

		int lengthMatches=countMatchesRecons(reconstructedReference,0);	
		int positionJump = al.Position+lengthMatches+numberOfDeletions(&al);
		if(hasCnoT[positionJump] && 
		   !thousandGenomesHasT[positionJump] )
		    isDeaminated=true; 

		// vcfr.repositionIterator(chrname,positionJump,positionJump);
		// while(vcfr.hasData()){
		//     SimpleVCF * toprint=vcfr.getData();
		//     //skip deletions in the alt
		//     if(toprint->getRef().length() != 1 )
		// 	continue;

		//     if(toprint->getRef()[0] != refeBase){
		// 	cerr<<reconstructedReference<<endl;
		// 	cerr<<al.Position<<endl;
		// 	cerr<<lengthMatches<<endl;
		// 	cerr<<numberOfDeletions(&al)<<endl;
		// 	cerr<<positionJump<<endl;
		// 	cerr<<"Problem6 position "<<*toprint<<" does not have a  "<<refeBase<<" as reference allele for read "<<al.Name<<endl;
		// 	exit(1);
		//     }

		//     if(  toprint->hasAtLeastOneC() && 
		// 	!toprint->hasAtLeastOneT() ){
		// 	isDeaminated=true; 
		//     }
		// }

	    }	

	   
	    
	}
		  



	totalReads++;

	if(isDeaminated){
	    deaminatedReads++;
	    writerDeam.SaveAlignment(al);		
	}else{
	    ndeaminatedReads++;
	    writerNoDeam.SaveAlignment(al);		
	}


    
    }//end for each read









    reader.Close();
    writerDeam.Close();
    writerNoDeam.Close();
    delete(hasCnoT);
    delete(hasGnoA);

    cerr<<"Program finished sucessfully, out of "<<totalReads<<" mapped reads (skipped: "<<skipped<<" reads) we flagged "<<deaminatedReads<<" as deaminated and "<<ndeaminatedReads<<" as not deaminated"<<endl;

   
    return 0;
}

