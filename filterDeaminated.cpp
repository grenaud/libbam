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





int main (int argc, char *argv[]) {

    int  minBaseQuality = 0;

    string usage=string(""+string(argv[0])+"  [in BAM file] [deam out BAM] [not deam out BAM]"+
			"\nThis program divides aligned read into potentially deaminated\n"+
			"\nreads and the puts the rest into another bam file.\n"+
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







    //iterating over the alignments for these regions
    BamAlignment al;
    int i;

    while ( reader.GetNextAlignment(al) ) {


	//skip unmapped
	if(!al.IsMapped())
	    continue;

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
	    if(  readBase  == 'A' && refeBase  == 'G' && int(al.Qualities[i]-offset) >= minBaseQuality){  isDeaminated=true; }

	    //second base next to 3'
	    i = 1;
	    refeBase=toupper(reconstructedReference[i]);
	    readBase=toupper(         al.QueryBases[i]);
	    if( readBase  == 'A' && refeBase  == 'G'  && int(al.Qualities[i]-offset) >= minBaseQuality){  isDeaminated=true; }

	    //last  base next to 5'
	    i = (al.QueryBases.length()-1) ;
	    refeBase=toupper(reconstructedReference[i]);
	    readBase=toupper(         al.QueryBases[i]);
	    if( readBase  == 'A' && refeBase  == 'G'  && int(al.Qualities[i]-offset) >= minBaseQuality){  isDeaminated=true; }


	}else{

		
	    //first base next to 5'
	    i = 0;
	    refeBase=toupper(reconstructedReference[i]);
	    readBase=toupper(         al.QueryBases[i]);
	    if( readBase  == 'T' && refeBase  == 'C'  && int(al.Qualities[i]-offset) >= minBaseQuality){  isDeaminated=true; }

	    //second last base next to 3'
	    i = (al.QueryBases.length()-2);
	    refeBase=toupper(reconstructedReference[i]);
	    readBase=toupper(         al.QueryBases[i]);
	    if( readBase  == 'T' && refeBase  == 'C'  && int(al.Qualities[i]-offset) >= minBaseQuality){  isDeaminated=true; }

	    //last base next to 5'
	    i = (al.QueryBases.length()-1);
	    refeBase=toupper(reconstructedReference[i]);
	    readBase=toupper(         al.QueryBases[i]);
	    if( readBase  == 'T' && refeBase  == 'C'  && int(al.Qualities[i]-offset) >= minBaseQuality){  isDeaminated=true; }	
	    
	}
		  



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

