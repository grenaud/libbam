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

#include "libgab.h"
#include "ReconsReferenceBAM.h"
#include "VCFreader.h"
#include "FastQParser.h"

using namespace std;
using namespace BamTools;

const int offset=33;


bool isTransition(char reference,char read){

    if(  (reference == 'A' && read == 'G') ||
	 (reference == 'G' && read == 'A') ||
	 (reference == 'C' && read == 'T') ||
	 (reference == 'T' && read == 'C') ){
	return true;
    }
	 
    return false;
}



bool isTransversions(char reference,char read){
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


bool hasIinfirstOrLastTwoBases(const string & reconstructedReference){
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

bool deletionsNextToTwo(const BamAlignment  * al){
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


void checkFasta(string * seq,int position,char charToHave,bool * deamFlagToChange){
    if(int(seq->length()) < position){
	cerr<<"ERROR cannot retrieve position "<<position<<" in a sequence of length "<<seq->length()<<endl;
	exit(1);
    }

    if( (*seq)[position-1] == charToHave){
	(*deamFlagToChange)=true;
    }// else do nothing
}
 

int countMatchesRecons(const string & reconstructedReference,int minusIndex){
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


void transformRef(char * refeBase,char * readBase){
    if( (*refeBase) == 'M'){
	(*refeBase)=(*readBase);
    }
    
}


//checks for an 'R' or 'S' for soft clip
bool hasBadCharacter(const string & reconstructedReference){

    for(unsigned int j=0;j<(reconstructedReference.length());j++){
	if(reconstructedReference[j] == 'R'  || 
	   reconstructedReference[j] == 'S' ){
	    return true;
	}
    }
    return false;
}

bool skipAlign(const string & reconstructedReference,const BamAlignment  * al,unsigned int * skipped){
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

    string usage=string(""+string(argv[0])+"  [in BAM file] [in fasta file] [chr name] [deam out BAM] [not deam out BAM]"+
			"\nThis program divides aligned read into potentially deaminated\n"+
			"\nreads and the puts the rest into another bam file if the deaminated positions are not called as the alternative base in the fasta.\n"+
			"\nTip: if you do not need one of them, use /dev/null as your output\n"+
			"arguments:\n"+
			"\t"+"--bq  [base qual]   : Minimum base quality to flag a deaminated site (Default: "+stringify(minBaseQuality)+")\n"+
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

    }

    string bamfiletopen    = string( argv[ argc-5 ] );
    string fastafiletopen  = string( argv[ argc-4 ] );
    string chrname         = string( argv[ argc-3 ] );
    string deambam         = string( argv[ argc-2 ] );
    string nondeambam      = string( argv[ argc-1 ] );

    //dummy reader, will need to reposition anyway
    // VCFreader vcfr (vcffiletopen,
    // 		    vcffiletopen+".tbi",
    // 		    chrname,
    // 		    1,
    // 		    1,
    // 		    5);
    FastQParser fqp (fastafiletopen,true);
    if(!fqp.hasData()){
	cerr << "Fasta file "<<fastafiletopen<<" is empty"<<endl;
	return 1; 
    }
    FastQObj * fastaObj	=fqp.getData();
    string * seqFromFasta=fastaObj->getSeq();

    // cout<<"seqFromFasta"<<endl;
    // cout<<*seqFromFasta<<endl;


    


    BamReader reader;
    
    if ( !reader.Open(bamfiletopen) ) {
    	cerr << "Could not open input BAM file"<< bamfiletopen << endl;
    	return 1;
    }

    // if ( !reader.LocateIndex()  ) {
    // 	cerr << "The index for the BAM file cannot be located" << endl;
    // 	return 1;
    // }

    // if ( !reader.HasIndex()  ) {
    // 	cerr << "The BAM file has not been indexed." << endl;
    // 	return 1;
    // }

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

	if(al.IsPaired() ){  
	    cerr<<"Paired end not yet coded"<<endl;
	    return 1;
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
		transformRef(&refeBase,&readBase);
		//if the fasta has a at least one G but no A
		checkFasta(seqFromFasta,al.Position+1,'G',&isDeaminated);


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
		transformRef(&refeBase,&readBase);

		//if the fasta has a at least one G but no A
		checkFasta(seqFromFasta,al.Position+2,'G',&isDeaminated);

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
		//if the fasta has a at least one G but no A
		checkFasta(seqFromFasta,positionJump,'G',&isDeaminated);
		

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
		transformRef(&refeBase,&readBase);

		//if the fasta has a at least one C but no T
		checkFasta(seqFromFasta,al.Position+1,'C',&isDeaminated);
		
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
		transformRef(&refeBase,&readBase);		
		int lengthMatches=countMatchesRecons(reconstructedReference,1);	
		int positionJump = al.Position+lengthMatches+numberOfDeletions(&al);
		//if the fasta has a at least one C but no T
		checkFasta(seqFromFasta,positionJump,'C',&isDeaminated);

		
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
		//if the fasta has a at least one C but no T
		checkFasta(seqFromFasta,positionJump,'C',&isDeaminated);


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






    if(fqp.hasData()){
    	cerr << "Fasta file "<<fastafiletopen<<" has more than 2 records"<<endl;
    	return 1; 
    }



    reader.Close();
    writerDeam.Close();
    writerNoDeam.Close();

    cerr<<"Program finished sucessfully, out of "<<totalReads<<" mapped reads (skipped: "<<skipped<<" reads) we flagged "<<deaminatedReads<<" as deaminated and "<<ndeaminatedReads<<" as not deaminated"<<endl;

   
    return 0;
}

