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
    bool pairedEnd=false;
    // unsigned int  binSize        = 1000;
    unsigned int  minBaseQuality = 1;
    // string fastaIndex = "/mnt/solexa/bin/gabriel/index.hg19.fai";

    string usage=string(""+string(argv[0])+"  [aligned BAM file] "+
			"\nThis program computes the number  of transitions\n"+
			"\nand transversions for each read\n"+
			"arguments:\n"+
			// // "\t"+"--fai [fasta index] : Fasta index (Default: "+stringify(fastaIndex)+")\n"+
			// "\t"+"--bin [bin sizes]   : Bin size to use (Default: "+stringify(binSize)+")\n"+
			"\t"+"--bq  [base qual]   : Minimum base quality (Default: "+stringify(minBaseQuality)+")\n"+
			"\n");

    if(argc == 1 ||
       (argc == 2 && (string(argv[0]) == "-h" || string(argv[0]) == "--help") )
       ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }


    for(int i=1;i<(argc);i++){ 

        //COMMON TO ALL
        //BEGIN loci selection

        // if(string(argv[i]) == "--bin"){
	//     binSize=destringify<unsigned int>(argv[i+1]);
        //     i++;
        //     continue;
	// }

        if(string(argv[i]) == "--bq"){
	    minBaseQuality=destringify<unsigned int>(argv[i+1]);
            i++;
            continue;
	}

    }

    string bamfiletopen = string( argv[ argc-1 ] );


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







    //iterating over the alignments for these regions
    BamAlignment al;
    //    unsigned int readCounter  =0;
    //    unsigned int sumReadLength=0;


    while ( reader.GetNextAlignment(al) ) {


	//skip unmapped
	if(!al.IsMapped())
	    continue;

	unsigned int numberTransitions=0;
	unsigned int numberTransversions=0;

	unsigned int numberA=0;
	unsigned int numberC=0;
	unsigned int numberG=0;
	unsigned int numberT=0;

		
	string reconstructedReference = reconstructRef(&al);
	
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
	    //that were called as 'A's with qual = 0 by Martin's BCL reader 
	    if(int(al.Qualities[i]-offset) < minBaseQuality ){
		continue;
	    }
		    
	    char refeBase=toupper(reconstructedReference[i]);
	    char readBase=toupper(al.QueryBases[i]);

	    if(refeBase == 'A'){
		numberA++;
	    }

	    if(refeBase == 'C'){
		numberC++;
	    }

	    if(refeBase == 'G'){
		numberG++;
	    }

	    if(refeBase == 'T'){
		numberT++;
	    }

		    

			
		    
	    if(isTransition( refeBase, readBase)){
		numberTransitions++;
	    }else{
		if(isTransversions( refeBase, readBase)){
		    numberTransversions++;
		}else{
		    cerr<<"Wrong substitution from "<<refeBase<<" to "<<readBase<<" in "<<al.Name<<endl<<"refe:"<<reconstructedReference<<endl<<"read:"<<al.QueryBases<<endl;
		    return 1;
		}
	    }
		       
		       
	    cout<<al.QueryBases.size()<<"\t"
		<<numberTransitions<<"\t"
		<<numberTransversions<<"\t"
		<<double(numberTransitions)/(double(numberTransitions)+double(numberTransversions) )<<"\t"
		<<double(numberC+numberG)/(double( numberA+numberC+numberG+numberT ) )
		<<endl;
	}

	// cout<<//referenceName<<":"<<":"<<coordinate<<"-"<<(coordinate+binSize)
	//     <<"\t"<<referenceLength
	//     <<"\t"<<readCounter
	//     <<"\t"<<sumReadLength 
	//     <<"\t"<<double(sumReadLength)/double(readCounter)
	//     <<"\t"<<numberTransitions
	//     <<"\t"<<numberTransversions
	//     <<"\t"<<double(numberTransitions)/(double(numberTransitions)+double(numberTransversions) )
	//     <<"\t"<<double(numberC+numberG)/(double( numberA+numberC+numberG+numberT ) )
	//     <<endl;

    
    }//end for each read









    reader.Close();




   
    return 0;
}

