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

using namespace std;
using namespace BamTools;




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





int main (int argc, char *argv[]) {
    bool pairedEnd=false;
    // unsigned int  binSize        = 1000;
    unsigned int  minBaseQuality = 0;
    // string fastaIndex = "/mnt/solexa/bin/gabriel/index.hg19.fai";

    string usage=string(""+string(argv[0])+"  [in BAM file] [deam out BAM] [not deam out BAM]"+
			"\nThis program divides aligned read into potentially deaminated\n"+
			"\nreads and the puts the rest into another bam file.\n"+
			"\nTip: if you do not need one of them, use /dev/null as your output\n"+
			"arguments:\n"+
			// // "\t"+"--fai [fasta index] : Fasta index (Default: "+stringify(fastaIndex)+")\n"+
			// "\t"+"--bin [bin sizes]   : Bin size to use (Default: "+stringify(binSize)+")\n"+
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
	    minBaseQuality=destringify<unsigned int>(argv[i+1]);
            i++;
            continue;
	}

    }

    string bamfiletopen = string( argv[ argc-3 ] );
    string deambam      = string( argv[ argc-2 ] );
    string nondeambam   = string( argv[ argc-1 ] );


    BamReader reader;
    
    if ( !reader.Open(bamfiletopen) ) {
    	cerr << "Could not open input BAM file"<< bamfiletopen << endl;
    	return 1;
    }

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



    // if ( !reader.LocateIndex()  ) {
    // 	cerr << "The index for the BAM file cannot be located" << endl;
    // 	return 1;
    // }

    // if ( !reader.HasIndex()  ) {
    // 	cerr << "The BAM file has not been indexed." << endl;
    // 	return 1;
    // }







    //iterating over the alignments for these regions
    BamAlignment al;
    //    unsigned int readCounter  =0;
    //    unsigned int sumReadLength=0;


    while ( reader.GetNextAlignment(al) ) {


	//skip unmapped
	if(!al.IsMapped())
	    continue;

	if(al.IsPaired() ){  
	    cerr<<"Paired end not yet coded"<<endl;
	    return 1;
	}

	//	cout<<"1 "<<al.Name<<"\t"<<al.QueryBases<<"\t"<<al.Qualities<<endl;

	string reconstructedReference = reconstructRef(&al);
	char refeBase;
	char readBase;
	bool isDeaminated;
	//	cout<<"2 "<<al.Name<<"\t"<<al.QueryBases<<"\t"<<al.Qualities<<endl;	
	//iterate over each character
	if(al.Qualities.size() != reconstructedReference.size()){
	    cerr<<"Quality line is not the same size as the reconstructed reference"<<endl;
	    return 1;
	}

	for(int i=0;i<al.QueryBases.size();i++){
	    //skip over matches and soft clipped bases, unresolved bases and deletion in the reference 
	    if(reconstructedReference[i] == 'S' ||
	       reconstructedReference[i] == 'M' || 
	       reconstructedReference[i] == 'I' || 
	       al.QueryBases[i]          == 'N' ){
		continue;
	    }

	    //Skip bases with low base quality
	    //this is done specifically to avoid unresolved base pairs
	    // cout<<al.Qualities[i]<<endl;
	    // cout<<int(al.Qualities[i]-33)<<endl;
	    // return 1;
	    if( int(al.Qualities[i]-33) < minBaseQuality ){
		continue;
	    }
		    
	    refeBase=toupper(reconstructedReference[i]);
	    readBase=toupper(al.QueryBases[i]);
	    isDeaminated=false;
	    //	    cout<<"3 "<<al.Name<<"\t"<<al.QueryBases<<"\t"<<al.Qualities<<endl;
	    

	    if(al.IsReverseStrand()){

		//first base next to 3'
		if(i == 0                          && readBase  == 'A' && refeBase  == 'G'){  isDeaminated=true; }
		//seonc base next to 3'
		if(i == 1                          && readBase  == 'A' && refeBase  == 'G'){  isDeaminated=true; }	       		    
		//last  base next to 5'
		if(i == (al.QueryBases.length()-1) && readBase  == 'A' && refeBase  == 'G'){  isDeaminated=true; }

	    }else{
		
		//first base next to 5'
		if(i == 0                          && readBase  == 'T' && refeBase  == 'C'){  isDeaminated=true; }
		//second last base next to 3'
		if(i == (al.QueryBases.length()-2) && readBase  == 'T' && refeBase  == 'C'){  isDeaminated=true; }	       		    
		//last base next to 5'
		if(i == (al.QueryBases.length()-1) && readBase  == 'T' && refeBase  == 'C'){  isDeaminated=true; }
		

	    }
		       
	   
	}//for each base
	
	if(isDeaminated){
	    writerDeam.SaveAlignment(al);		
	}else{
	    writerNoDeam.SaveAlignment(al);		
	}


    
    }//end for each read









    reader.Close();
    writerDeam.Close();
    writerNoDeam.Close();



   
    return 0;
}

