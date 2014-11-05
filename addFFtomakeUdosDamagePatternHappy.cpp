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

    if(argc != 3){
	cerr<<"This program adds a FF:i:3 tag to all reads to make Udo's damage-pattern happy"<<endl;
	cerr<<"Usage "<<argv[0]<<" [bam file in] [bam file out]"<<endl;
	return 1;
    }

    string bamfiletopen  = string(argv[1]);
    string bamfiletwrite = string(argv[2]);


    cerr<<"Reading "<<bamfiletopen<<" writing to "<<bamfiletwrite<<" with FF:i:3"<<endl;

    BamReader reader;
    BamWriter writer;

    if ( !reader.Open(bamfiletopen) ) {
    	cerr << "Could not open input BAM files." << endl;
    	return 1;
    }

    SamHeader  myHeader=reader.GetHeader();
    SamProgram sp;

    string pID          = "addFF";   
    string pName        = "addFF";   
    string pCommandLine = "";
    for(int i=0;i<(argc);i++){
	pCommandLine += (string(argv[i])+" ");
    }
    putProgramInHeader(&myHeader,pID,pName,pCommandLine,returnGitHubVersion(string(argv[0]),"."));

    // SamReadGroupDictionary  srgd;
    // //for(unsigned int lane=1;lane<=8;lane++){
    // SamReadGroup srg ( RGGROUP );
    // srg.Sample   = RGGROUP;
    // srgd.Add( srg );  
    // //}

    // myHeader.ReadGroups=srgd;




    if( !writer.Open(bamfiletwrite,myHeader,reader.GetReferenceData() ) ) {
    	cerr << "Could not open output BAM file  "<<bamfiletwrite << endl;
    	return 1;	
    }

    
    BamAlignment al;
    while ( reader.GetNextAlignment(al) ) {
	//cerr<<al.Name<<endl;
	if(!al.EditTag("FF","i",3)){
	    cerr << "Unable to edit FF tag" << endl;
	    exit(1);     
	}
	if(!al.HasTag("FF")){
	    cerr << "Unable to edit FF tag" << endl;
	    exit(1);     
	}
	writer.SaveAlignment(al);
    }

    reader.Close();
    writer.Close();

    cerr<<"Program "<<argv[0]<<" terminated gracefully"<<endl;

    return 0;
}

