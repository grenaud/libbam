/*
 * cutUMI
 * Date: Dec-06-2019 
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
	 cout<<"Usage: cutUMI [in bam] [outbam]"<<endl<<"This program cuts the UMIs"<<endl;
	 return 1;
    }

     string bamfiletopen = string(argv[1]);
     string bamFileOUT   = string(argv[2]);

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
    string rx;
    string qx;

    while ( reader.GetNextAlignment(al) ) {

	if(al.HasTag("RX")){
	    al.GetTag("RX",rx);
	    al.RemoveTag("RX");
	    rx = rx.substr(0,6);
	    if(!al.EditTag("RX","Z",rx)){
		cerr << "Unable to edit RX tag" << endl;
		exit(1);     
	    }	   
	}

	if(al.HasTag("QX")){
	    al.GetTag("QX",qx);
	    al.RemoveTag("QX");
	    qx = qx.substr(0,6);
	    if(!al.EditTag("QX","Z",qx)){
		cerr << "Unable to edit QX tag" << endl;
		exit(1);     
	    }	   
	}


	writer.SaveAlignment(al);

    } //while al

    reader.Close();
    writer.Close();

    return 0;
}

