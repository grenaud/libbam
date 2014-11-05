/*
 * removeIndices
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

    bool removeP7=true;
    bool removeP5=true;

     if( (argc== 1) ||
    	(argc== 2 && string(argv[1]) == "-h") ||
    	(argc== 2 && string(argv[1]) == "-help") ||
    	(argc== 2 && string(argv[1]) == "--help") ){
	 cout<<"Usage:removeRG (options) [in bam] [outbam]"<<endl<<"This program removes the indices"<<endl<<endl<<
	     "Options"<<endl<<endl<<
	     "\t-p7\tOnly remove the P7 indices"<<endl<<
	     "\t-p5\tOnly remove the P5 indices"<<endl<<endl;
	 return 1;
    }

    int lastOpt=1;



     for(int i=1;i<(argc);i++){ //all but the last 3 args

	if(string(argv[i])[0] != '-'  ){
	    lastOpt=i;
	    break;
	}


	if(string(argv[i]) == "-p7"  ){
	    removeP7=true;
	    removeP5=false;
	    continue;
	}

	if(string(argv[i]) == "-p5"  ){
	    removeP7=false;
	    removeP5=true;
	    continue;
	}

	cerr<<"Wrong option "<<string(argv[i])<<endl;
	return 1;

     }
     string bamfiletopen = string(argv[lastOpt+0]);
     string bamFileOUT   = string(argv[lastOpt+1]);
     cerr<<"reading from "<<bamfiletopen<<" writing to "<<bamFileOUT<<endl;
     //     return 1;
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
	if(removeP7){
	    if(al.HasTag("XI")){
		al.RemoveTag("XI");
	    }
	    if(al.HasTag("YI")){
		al.RemoveTag("YI");
	    }
	}

	if(removeP5){
	    if(al.HasTag("XJ")){
		al.RemoveTag("XJ");
	    }
	    if(al.HasTag("YJ")){
		al.RemoveTag("YJ");
	    }
	}

	writer.SaveAlignment(al);

    } //while al

    reader.Close();
    writer.Close();

cerr<<"Program "<<argv[0]<<" terminated succesfully"<<endl;
    return 0;
}

