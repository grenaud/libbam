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
	 cout<<"Usage:failQualPair [in bam] [outbam]"<<endl<<"this program takes all pairs in a bam file, if one fails the QC flag, the other will fail as well"<<endl;
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

	// if(al.HasTag("NM") || al.HasTag("MD")  ){
	//     cerr << "Reads should not be aligned" << endl;
	//     return 1;
	// }


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

		//if one fails the QC, both of them fail
		if(al.IsFailedQC() || al2.IsFailedQC() ){
		    al.SetIsFailedQC(true);
		    al2.SetIsFailedQC(true);
		}
		writer.SaveAlignment(al2);
		writer.SaveAlignment(al);

		
		//  SINGLE END, nothing to do
	    }else{ 
		writer.SaveAlignment(al);
	    } //end single end

	    al2Null=true;
	}//second pair
		    

    } //while al

    reader.Close();
    writer.Close();

    return 0;
}

