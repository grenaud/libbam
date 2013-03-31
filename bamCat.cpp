/*
 * bamCat
 * Date: March 31 2013
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
	 cout<<"Usage:"<<argv[0]<<" [out bam] [inbam1] [inbam2] ..."<<endl<<"this program appends all BAM records into a single file using the header from the first file"<<endl;
    	return 1;
    }



     string bamFileOUT   = string(argv[1]);
     BamWriter writer;

     for(int bidx=2;bidx<argc;bidx++){

	 BamReader reader;
	 string bamfiletopen = string(argv[bidx]);
	 cerr<<"Opening "<<bamfiletopen<<endl;
	 if ( !reader.Open(bamfiletopen) ) {
	     cerr << "Could not open input BAM files." << endl;
	     return 1;
	 }
	 if(bidx==2){
	     const SamHeader header     = reader.GetHeader();
	     const RefVector references = reader.GetReferenceData();
    
	     if ( !writer.Open(bamFileOUT,header,references) ) {
		 cerr << "Could not open output BAM file "<<bamFileOUT << endl;
		 return 1;
	     }
	 }

	 BamAlignment al;
 
	 while ( reader.GetNextAlignment(al) ) {
	     writer.SaveAlignment(al);	    
	 } //while al

	 reader.Close();
     }



     writer.Close();
     cerr<<"Program "<<argv[0]<<"finished gracefully"<<endl;
    return 0;
}

