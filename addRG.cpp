/*
 * addRG
 * Date: Jan-28-2013 
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


int main (int argc, char *argv[]) {

    if(argc != 4){
	cerr<<"This program adds a RG tag to all reads to make GATK happy"<<endl;
	cerr<<"Usage "<<argv[0]<<" [bam file in] [bam file out] [RG]"<<endl;
	return 1;
    }

    string bamfiletopen  = string(argv[1]);
    string bamfiletwrite = string(argv[2]);
    string RGGROUP       = string(argv[3]);

    cerr<<"Reading "<<bamfiletopen<<" writing to "<<bamfiletwrite<<" with RG:Z:"<<RGGROUP<<endl;

    BamReader reader;
    BamWriter writer;

    if ( !reader.Open(bamfiletopen) ) {
    	cerr << "Could not open input BAM files." << endl;
    	return 1;
    }

    SamHeader  myHeader=reader.GetHeader();
    SamProgram sp;

    string pID          = "addRG";   
    string pName        = "addRG";   
    string pCommandLine = "";
    for(int i=0;i<(argc);i++){
	pCommandLine += (string(argv[i])+" ");
    }
    putProgramInHeader(&myHeader,pID,pName,pCommandLine,returnGitHubVersion(string(argv[0]),"."));



    if( !writer.Open(bamfiletwrite,reader.GetHeader(),reader.GetReferenceData() ) ) {
    	cerr << "Could not open output BAM file  "<<bamfiletwrite << endl;
    	return 1;	
    }

    
    BamAlignment al;
    while ( reader.GetNextAlignment(al) ) {
	//cerr<<al.Name<<endl;
	if(!al.EditTag("RG","Z",RGGROUP)){
	    cerr << "Unable to edit ZQ tag" << endl;
	    exit(1);     
	}
	if(!al.HasTag("RG")){
	    cerr << "Unable to edit ZQ tag" << endl;
	    exit(1);     
	}
	writer.SaveAlignment(al);
    }

    reader.Close();
    writer.Close();

    cerr<<"Program "<<argv[0]<<" terminated gracefully"<<endl;

    return 0;
}

