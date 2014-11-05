/*
 * readWhereEitherIsMapped
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
	 cout<<"Usage:readWhereEitherIsMapped [in bam] [outbam]"<<endl<<"This program takes the first bam file and write the pairs where one mapped"<<endl;
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
    BamAlignment al2;
    bool al2Null=true;
 
    while ( reader.GetNextAlignment(al) ) {


	if(al.IsPaired() && 
	   al2Null ){
	    al2=al;
	    al2Null=false;
	    continue;
	}else{
	    if(al.IsPaired() && 
	       !al2Null){
		if(al.Name != al2.Name ){
		    cerr << "Seq#1 "<<al.Name<<" has a different id than seq #2 "<<al2.Name<<", exiting " << endl;
		    return 1;
		} 

		//if at least maps
		if(al.IsMapped() || 
		   al2.IsMapped() ){
		    writer.SaveAlignment(al2);
		    writer.SaveAlignment(al);		    
		}
		
		
	    }else{  //  SINGLE END
		if( al.IsMapped() )
		    writer.SaveAlignment(al);
	    } //end single end

	    al2Null=true;
	}//second pair
		    

    } //while al

    reader.Close();
    writer.Close();

    return 0;
}

