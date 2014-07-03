/*
 * retrieveMapped_single_and_ProperlyPair
 * Date: Oct-10-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */


#include <iostream>

#include "api/BamMultiReader.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"

#include "PutProgramInHeader.h"
#include "utils.h"

using namespace std;
using namespace BamTools;


int main (int argc, char *argv[]) {

    if( (argc== 1) ||
    	(argc== 2 && string(argv[1]) == "-h") ||
    	(argc== 2 && string(argv[1]) == "-help") ||
    	(argc== 2 && string(argv[1]) == "--help") ){
	cout<<"Usage:removeRG [in bam] [outbam]"<<endl<<"This program keeps reads that are mapped and single or properly mapped"<<endl;
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

    SamHeader header = reader.GetHeader();
    string pID          = "retrieveMapped_single_and_ProperlyPair";   
    string pName        = "retrieveMapped_single_and_ProperlyPair";   
    string pCommandLine = "";
    for(int i=0;i<(argc);i++){
        pCommandLine += (string(argv[i])+" ");
    }
    putProgramInHeader(&header,pID,pName,pCommandLine);

    const RefVector references = reader.GetReferenceData();
    if ( !writer.Open(bamFileOUT,header,references) ) {
    	cerr << "Could not open output BAM file "<<bamFileOUT << endl;
    	return 1;
    }

    BamAlignment al;
    uint64_t total=0;
    uint64_t kept=0;

    while ( reader.GetNextAlignment(al) ) {
	total++;

	if(al.IsPaired() ){
	    if(!al.IsProperPair() )
		continue;
	    if(!al.IsMateMapped() )
		continue;
	}

	if(!al.IsMapped()){
	    continue;
	}

	kept++;
	writer.SaveAlignment(al);

    } //while al

    reader.Close();
    writer.Close();
    cerr<<"retrieveMapped_single_and_ProperlyPair: out of "<<total<<" sequences we kept "<<kept<<" sequences "<<endl;

    return 0;
}

