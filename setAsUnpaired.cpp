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
#include "utils.h"

using namespace std;
using namespace BamTools;


int main (int argc, char *argv[]) {

     if( (argc== 1) ||
    	(argc== 2 && string(argv[1]) == "-h") ||
    	(argc== 2 && string(argv[1]) == "-help") ||
    	(argc== 2 && string(argv[1]) == "--help") ){
	 cout<<"Usage:setAsUnpaired [in bam] [outbam]"<<endl<<"this program takes flags all paired sequences as singles"<<endl;
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
 
    while ( reader.GetNextAlignment(al) ) {
	if(al.IsMapped()){
	    cerr << "Cannot yet handle mapped reads " << endl;
	    return 1;
	}

	
	al.SetIsPaired (false);
	
	writer.SaveAlignment(al);    

    } //while al

    reader.Close();
    writer.Close();

    return 0;
}

