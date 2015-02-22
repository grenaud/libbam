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
#include "PutProgramInHeader.h"
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
    int  bpToDecrease   = 5;

    string usage=string(""+string(argv[0])+"  [in BAM file] [deam out BAM] [not deam out BAM]"+
			"\nThis program divides aligned read into potentially deaminated\n"+
			"\nreads and the puts the rest into another bam file.\n"+
			"\nTip: if you do not need one of them, use /dev/null as your output\n"+
			"arguments:\n"+
			"\t"+"--bq  [base qual]   : Minimum base quality to flag a deaminated site    (Default: "+stringify(minBaseQuality)+")\n"+
			"\t"+"-n    [# bases]     : Consider the nth bases surrounding the 5'/3' ends (Default:"+stringify(bpToDecrease)+") "+"\n"+

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

        if(string(argv[i]) == "-n"){
            bpToDecrease =destringify<int>(argv[i+1]);
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


    const RefVector references = reader.GetReferenceData();

    SamHeader header = reader.GetHeader();
    string pID          = "filterDeaminatedDouble";   
    string pName        = "filterDeaminatedDouble";   
    string pCommandLine = "";
    for(int i=0;i<(argc);i++){
        pCommandLine += (string(argv[i])+" ");
    }
    putProgramInHeader(&header,pID,pName,pCommandLine);

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
    //int i;

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


	if(al.IsPaired() ){  


	    // cerr << "We do not support paired end" << endl;
	    // return 1;

	    if(al.IsReverseStrand()){ //is reverse complemented, we decrease the last As regardless of first or second mate
                    
		for(int indexToCheck=(al.QueryBases.length()-1);indexToCheck>( (al.QueryBases.length()-1)-bpToDecrease);indexToCheck--){

		    refeBase=toupper(reconstructedReference[indexToCheck]);
		    readBase=toupper(         al.QueryBases[indexToCheck]);
		    if(  readBase  == 'A' && refeBase  == 'G' && int(al.Qualities[indexToCheck]-offset) >= minBaseQuality){  isDeaminated=true; }

		}

	    }else{ //if not reverse complemented, we decrease the first 5 Ts regardless of first or second mate

		for(int indexToCheck=0;indexToCheck<bpToDecrease;indexToCheck++){

		    refeBase=toupper(reconstructedReference[indexToCheck]);
		    readBase=toupper(         al.QueryBases[indexToCheck]);
		    if( readBase  == 'T' && refeBase  == 'C'  && int(al.Qualities[indexToCheck]-offset) >= minBaseQuality){  isDeaminated=true; }

		}


	    }
        
	}//end of paired end
	else{//we consider single reads to have been sequenced from 5' to 3'


	    if(al.QueryBases.length() <= 5){
		cerr << "We do not process reads with less than 5bp" << endl;
		return 1;
	    }

	    for(int indexToCheck=0;indexToCheck<bpToDecrease;indexToCheck++){


		refeBase=toupper(reconstructedReference[indexToCheck]);
		readBase=toupper(         al.QueryBases[indexToCheck]);
		if( readBase  == 'T' && refeBase  == 'C'  && int(al.Qualities[indexToCheck]-offset) >= minBaseQuality){  isDeaminated=true; }
		// if(toupper(al.QueryBases[indexToCheck]) == 'T'){
		//     al.Qualities[indexToCheck]=char(offset+baseQualForDeam);
		// }
	    }

	    for(int indexToCheck=(al.QueryBases.length()-1);indexToCheck>( (al.QueryBases.length()-1) -bpToDecrease);indexToCheck--){

		refeBase=toupper(reconstructedReference[indexToCheck]);
		readBase=toupper(         al.QueryBases[indexToCheck]);
		if(  readBase  == 'A' && refeBase  == 'G' && int(al.Qualities[indexToCheck]-offset) >= minBaseQuality){  isDeaminated=true; }
		    
		// if(toupper(al.QueryBases[indexToCheck]) == 'A'){
		//     al.Qualities[indexToCheck]=char(offset+baseQualForDeam);
		// }
	    }
                

	}//end of single end
		  



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

