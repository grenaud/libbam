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

unsigned int total;
unsigned int unmapped;
unsigned int nocigar;

bool passed(BamAlignment * al){
    if(!al->IsMapped()){//leave out unmapped reads
	unmapped++;
	return false;
    }

    if(al->CigarData.size() == 0){//no cigar =leave out
	nocigar++;
	return false;
    }
    return true;
}

int main (int argc, char *argv[]) {

     if( (argc== 1) ||
    	(argc== 2 && string(argv[1]) == "-h") ||
    	(argc== 2 && string(argv[1]) == "-help") ||
    	(argc== 2 && string(argv[1]) == "--help") ){
	 cout<<"Usage:removeUnalignedANDWrongCigar [in bam] [outbam]"<<endl<<"this removes reads that have a flag set as unaligned or have a '*' in their CIGAR line\nuseful prior to calling GATK"<<endl;
    	return 1;
    }

     string bamfiletopen = string(argv[1]);
     string bamFileOUT   = string(argv[2]);

     BamReader reader;
     BamWriter writer;
     total=0;
     unmapped=0;
     nocigar=0;

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

	// if(al.IsPaired() && 
	//    al2Null ){
	//     al2=al;
	//     al2Null=false;
	//     continue;
	// }else{
	//     if(al.IsPaired() && 
	//        !al2Null){
	// 	if(al.Name != al2.Name ){
	// 	    cerr << "removeUnalignedANDWrongCigar: Seq#1 "<<al.Name<<" has a different id than seq #2 "<<al2.Name<<"exiting " << endl;
	// 	    return 1;
	// 	} 
	// 	total++;


	// 	if(passed(&al) &&  passed(&al2) ){
	// 	    writer.SaveAlignment(al2);
	// 	    writer.SaveAlignment(al);
	// 	}
	
		
	//     }else{  //  SINGLE END
	total++;
	if(passed(&al)  )
	    writer.SaveAlignment(al);
	//} //end single end
    
	//     al2Null=true;
	// }
		    

    } //while al


    reader.Close();
    writer.Close();

    cerr<<"removeUnalignedANDWrongCigar: out of "<<total<<" we discarded "<<unmapped<<" due to unmapped "<<nocigar<<" due to no cigar"<<endl;

    return 0;
}

