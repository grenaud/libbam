/*
 * filterEditDist
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
	cout<<"Usage:filterEditDist [in bam] [outbam]"<<endl<<"This program removes reads with a high edit distance"<<endl;
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
    unsigned total=0;
    unsigned kept=0;
    int editDist=0;
    while ( reader.GetNextAlignment(al) ) {
	total++;
	
	if(al.HasTag("NM") ){
	    editDist=0;
	    al.GetTag("NM",editDist);
	    // if(al.Name == "SN7001204_0130_AC0M6HACXX_PEdi_SS_L9302_L9303_1:8:2214:4172:88206" ){
	    // 	cerr<<al.Name<<endl;
	    // 	cerr<<double(editDist)<<endl;
	    // 	cerr<<double(al.QueryBases.length())<<endl;
	    // 	cerr<< int(0.2*double(al.QueryBases.length())+0.5 )<<endl;

	    // 	cerr<<((editDist) > int(0.2*double(al.QueryBases.length()) +0.5) )<<endl;
	    // }
	    if((editDist) > int(0.2*double(al.QueryBases.length())+0.5 )){ //This is wrong but to ensure compatibility with python
		continue;
	    }
	}

	kept++;
	writer.SaveAlignment(al);

    } //while al

    reader.Close();
    writer.Close();
    cerr<<"filterEditDist: out of "<<total<<" sequences we kept "<<kept<<" sequences "<<endl;

    return 0;
}

