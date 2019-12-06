/*
 * removeTagsMapping
 * Date: Apr-19-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */


#include <iostream>

#include "api/BamMultiReader.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"
#include "utils.h"
#include "PutProgramInHeader.h"

using namespace std;
using namespace BamTools;

const uint32_t flagSingleReads =  4; // 00000100


int main (int argc, char *argv[]) {

    if(argc != 3){
	cerr<<"This program sets unmapped first and second pairs as single-end"<<endl;
	cerr<<"Usage "<<argv[0]<<" [bam file in] [bam file out]"<<endl;
	return 1;
    }

    string bamfiletopen  = string(argv[1]);
    string bamfiletwrite = string(argv[2]);

    cerr<<"Reading "<<bamfiletopen<<" writing to "<<bamfiletwrite<<endl;

    BamReader reader;
    BamWriter writer;

    if ( !reader.Open(bamfiletopen) ) {
    	cerr << "Could not open input BAM files." << endl;
    	return 1;
    }

    SamHeader  myHeader=reader.GetHeader();
    SamProgram sp;

    string pID          = "pairsAsSingle";   
    string pName        = "pairsAsSingle";   
    string pCommandLine = "";
    for(int i=0;i<(argc);i++){
	pCommandLine += (string(argv[i])+" ");
    }
    putProgramInHeader(&myHeader,pID,pName,pCommandLine,returnGitHubVersion(string(argv[0]),"."));

    //no @SQ
    myHeader.Sequences.Clear();
    vector< RefData > 	emptyRefVector;

    if( !writer.Open(bamfiletwrite,myHeader,emptyRefVector ) ) {
    	cerr << "Could not open output BAM file  "<<bamfiletwrite << endl;
    	return 1;	
    }

    
    BamAlignment al;
    while ( reader.GetNextAlignment(al) ) {
	//deleting tag data

	//reset the flag
	if(al.IsPaired()){
	    if(al.IsFirstMate()){
		al.AlignmentFlag =  flagSingleReads;
		al.Name          =  al.Name+"_1";
	    }else{
		al.AlignmentFlag =  flagSingleReads;
		al.Name          =  al.Name+"_2";
	    }
	}else{
	    //fine
	}


	writer.SaveAlignment(al);
    }

    reader.Close();
    writer.Close();

    cerr<<"Program "<<argv[0]<<" terminated gracefully"<<endl;

    return 0;
}

