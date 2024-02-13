/*
 * filterNoNonDNAbase
 * Date: Jan-28-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */


#include <iostream>

#include "api/BamMultiReader.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"
#include "libgab.h"
#include "PutProgramInHeader.h"

using namespace std;
using namespace BamTools;


int main (int argc, char *argv[]) {

    if(argc != 3){
	cerr<<"This program retains only reads with bases A,C,G,T,N"<<endl;
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

    string pID          = "filterNoNonDNAbase";   
    string pName        = "filterNoNonDNAbase";   
    string pCommandLine = "";
    for(int i=0;i<(argc);i++){
	pCommandLine += (string(argv[i])+" ");
    }
    putProgramInHeader(&myHeader,pID,pName,pCommandLine,returnGitHubVersion(string(argv[0]),"."));

    // SamReadGroupDictionary  srgd;
    // myHeader.ReadGroups=srgd;

    if( !writer.Open(bamfiletwrite,myHeader,reader.GetReferenceData() ) ) {
    	cerr << "Could not open output BAM file  "<<bamfiletwrite << endl;
    	return 1;	
    }

    
    BamAlignment al;
    while ( reader.GetNextAlignment(al) ) {
	bool nonDNA=false;
	for(int i=0;i<al.QueryBases.size();i++){
	    if(!isValidDNA(al.QueryBases[i] )){
		nonDNA=true;
		break;
	    }
	}
	
	if(!nonDNA)
	    writer.SaveAlignment(al);
    }

    reader.Close();
    writer.Close();

    cerr<<"Program "<<argv[0]<<" terminated gracefully"<<endl;

    return 0;
}

