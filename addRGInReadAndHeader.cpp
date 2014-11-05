/*
 * addRGInReadAndHeader
 * Date: Apr-23-2013
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
	cerr<<"This program adds a RG tag to all reads and the header"<<endl;
	cerr<<"Usage "<<argv[0]<<"  [RG] [bam file in] [bam file out]"<<endl;
	return 1;
    }

    string RGGROUP       = string(argv[1]);
    string bamfiletopen  = string(argv[2]);
    string bamfiletwrite = string(argv[3]);

    cerr<<"Reading "<<bamfiletopen<<" writing to "<<bamfiletwrite<<" with RG:Z:"<<RGGROUP<<endl;

    BamReader reader;
    BamWriter writer;

    if ( !reader.Open(bamfiletopen) ) {
    	cerr << "Could not open input BAM files." << endl;
    	return 1;
    }

    SamHeader  myHeader=reader.GetHeader();
    SamProgram sp;



    //@RG 
    SamReadGroupDictionary  srgd;
    SamReadGroup srg ( RGGROUP );
    srg.Sample       = RGGROUP;
    srgd.Add( srg );  
    

    myHeader.ReadGroups=srgd;



    if( !writer.Open(bamfiletwrite,myHeader,reader.GetReferenceData() ) ) {
    	cerr << "Could not open output BAM file  "<<bamfiletwrite << endl;
    	return 1;	
    }

    
    BamAlignment al;
    while ( reader.GetNextAlignment(al) ) {
	//cerr<<al.Name<<endl;
	// string rgtoadd=RGGROUP;
	// vector<string> name = allTokens(al.Name,':');
	// rgtoadd+="_"+name[1];
	
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

