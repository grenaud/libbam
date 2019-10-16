/*
 * retrieveMappedCertainIsize
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
	cerr<<"Usage:removeRG [min insertsize]  [max insertsize] [in bam] [outbam]"<<endl<<"This program keeps reads that are mapped and single or properly mapped"<<endl;
	return 1;
    }

    string mininsertsizeS   = string(argv[1]);
    string maxinsertsizeS   = string(argv[2]);

    string bamfiletopen  = string(argv[3]);
    string bamFileOUT    = string(argv[4]);

    int mininsertsize = destringify<int>(  mininsertsizeS );
    int maxinsertsize = destringify<int>(  maxinsertsizeS );
    // cerr<<"Min insert size"<< mininsertsize<<endl;
    // cerr<<"Max insert size"<< maxinsertsize<<endl;

    if(maxinsertsizeS <= mininsertsizeS){
	cerr<<"Usage:removeRG [min insertsize]  [max insertsize] [in bam] [outbam]"<<endl<<"This program keeps reads that are mapped and single or properly mapped"<<endl<<"ERROR: the min needs to be smaller than the max"<<endl;
	return 1;
    }
    BamReader reader;
    BamWriter writer;

    if ( !reader.Open(bamfiletopen) ) {
    	cerr << "Could not open input BAM files:" << bamfiletopen<<endl;
    	return 1;
    }

    SamHeader header = reader.GetHeader();
    string pID          = "retrieveMappedCertainIsize";   
    string pName        = "retrieveMappedCertainIsize";   
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

	if(!al.IsMapped()){
	    continue;
	}

	if(al.IsPaired()){

	    if(!al.IsProperPair() )
		continue;
	    if(!al.IsMateMapped() )
		continue;
	    //}

	    //if( al.IsFirstMate()  ){
	    if(al.InsertSize > 0){

		    if( al.InsertSize < mininsertsize)
			continue;
		    
		    if( al.InsertSize > maxinsertsize)
			continue;

	    }else{
		//cout<<(-1.0*al.InsertSize)<<endl;
		
		if( (-1.0*al.InsertSize) < mininsertsize)
		    continue;
		
		if( (-1.0*al.InsertSize) > maxinsertsize)
		    continue;

	    }
	}else{

	    if(al.Length < mininsertsize)
		continue;

	    if(al.Length > maxinsertsize)
		continue;

	 }

	kept++;
	writer.SaveAlignment(al);

    } //while al

    reader.Close();
    writer.Close();
    cerr<<"retrieveMappedCertainIsize: out of "<<total<<" sequences we kept "<<kept<<" sequences "<<endl;

    return 0;
}

