/*
 * failQualPair
 * Date: Oct-10-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */


#include <iostream>

#include "api/BamMultiReader.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"
#include "libgab.h"

using namespace std;
using namespace BamTools;


int main (int argc, char *argv[]) {

     if( (argc== 1) ||
    	(argc== 2 && string(argv[1]) == "-h") ||
    	(argc== 2 && string(argv[1]) == "-help") ||
    	(argc== 2 && string(argv[1]) == "--help") ){
	 cout<<"Usage: replaceRG [NEWRG] [in bam] [outbam]"<<endl<<"This program removes RG fields and replaces them with [NEWRG]"<<endl;
	 return 1;
    }
     string rgtoadd      = string(argv[1]);
     string bamfiletopen = string(argv[2]);
     string bamFileOUT   = string(argv[3]);

     BamReader reader;
     BamWriter writer;

     if ( !reader.Open(bamfiletopen) ) {
    	cerr << "Could not open input BAM files." << endl;
    	return 1;
     }
    const SamHeader header = reader.GetHeader();
    const RefVector references = reader.GetReferenceData();
    if ( !writer.Open(bamFileOUT,header,references) ) {
    	cerr << "Could not open output BAM file "<<bamFileOUT << endl;
    	return 1;
    }

    BamAlignment al;
 
    while ( reader.GetNextAlignment(al) ) {

	if(al.HasTag("RG")){
	    al.RemoveTag("RG");
	}

	if(!al.EditTag("RG","Z",rgtoadd)){
	    cerr << "Unable to edit ZQ tag" << endl;
	    exit(1);     
	}

	writer.SaveAlignment(al);

    } //while al

    reader.Close();
    writer.Close();

    return 0;
}

