/*
 * addRG_CTEAM
 * Date: Apr-23-2013
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

    if(argc != 4){
	cerr<<"This program adds a RG tag to all reads depending on the lane"<<endl;
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

    string pID          = "addRG_CTEAM";   
    string pName        = "addRG_CTEAM";   
    string pCommandLine = "";
    for(int i=0;i<(argc);i++){
	pCommandLine += (string(argv[i])+" ");
    }
    putProgramInHeader(&myHeader,pID,pName,pCommandLine,returnGitHubVersion(string(argv[0]),"."));


    //@RG 
    SamReadGroupDictionary  srgd;
    for(unsigned int lane=1;lane<=8;lane++){
	SamReadGroup srg ( RGGROUP+"_"+stringify(lane) );
	srg.Sample   = RGGROUP;
        srgd.Add( srg );  
    }

    myHeader.ReadGroups=srgd;





    if( !writer.Open(bamfiletwrite,myHeader,reader.GetReferenceData() ) ) {
    	cerr << "Could not open output BAM file  "<<bamfiletwrite << endl;
    	return 1;	
    }

    
    BamAlignment al;
    while ( reader.GetNextAlignment(al) ) {
	//cerr<<al.Name<<endl;
	string rgtoadd=RGGROUP;
	vector<string> name = allTokens(al.Name,':');
	rgtoadd+="_"+name[1];
	
	if(!al.EditTag("RG","Z",rgtoadd)){
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

