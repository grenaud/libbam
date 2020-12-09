/*
 * failQualPair
 * Date: Oct-10-2012 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <string.h>
#include <iostream>

#include "api/BamMultiReader.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"
#include "libgab.h"

using namespace std;
using namespace BamTools;

unsigned int discardedLength;
unsigned int discardedEdit  ;
unsigned int total          ;

int    minLength   ;
double fractionEdit;

bool passed(BamAlignment * al){
    if(int(al->QueryBases.length()) < minLength){
	return false;
	discardedLength++;
    }

    if(al->IsMapped()){
	int editDist;
	if(!al->GetTag("NM",editDist)){
	    cerr<<"Cannot get the NM field from "<<al->Name<<endl;
	}
	if(  double(editDist) > (fractionEdit*double(al->QueryBases.length())) ){
	    return false;
	    discardedEdit++;
	}
    }
    return true;
}

int main (int argc, char *argv[]) {

    minLength   =35;
    fractionEdit=0.2;

    const string usage=string(string(argv[0])+" [options] input.bam out.bam"+"\n\n"+
			      "This program takes a BAM file as input and produces\n"+
			      "another without the reads with high edit distance without\n"+
			      "short sequences\n"+
			      "\n"+
			      "Options:\n"+
			      "\t"+"-l , --length" +"\n\t\t"+"Minimum length of sequences (Default: "+stringify(minLength)+")"+"\n"+
			      "\t"+"-f , --frac"   +"\n\t\t"+"Fraction of the length as maximum edit distance (Default: "+stringify(fractionEdit)+")"+"\n");
		


     if( (argc== 1) ||
    	(argc== 2 && string(argv[1]) == "-h") ||
    	(argc== 2 && string(argv[1]) == "-help") ||
    	(argc== 2 && string(argv[1]) == "--help") ){
	 cout<<"Usage: "<<usage<<endl;
	 return 1;
     }

     for(int i=1;i<(argc-2);i++){ //all but the last two args
	 //cout<<"argv["<<i<<"] "<<argv[i]<<endl;
	 if(strcmp(argv[i],"-l") == 0 || strcmp(argv[i],"--length") == 0 ){
	     minLength=destringify<int>(argv[i+1]);
	     i++;
	     continue;
	 }

	 if(strcmp(argv[i],"-f") == 0 || strcmp(argv[i],"--frac") == 0 ){
	     fractionEdit=destringify<double>(argv[i+1]);
	     i++;
	     continue;
	 }

     }

     discardedLength = 0;
     discardedEdit   = 0;
     total           = 0;

     string bamfiletopen = string(argv[argc-2]);
     string bamFileOUT   = string(argv[argc-1]);

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
    // BamAlignment al2;
    // bool al2Null=true;
 
    while ( reader.GetNextAlignment(al) ) {
	total++;
	if(!passed(&al)  )
	    al.SetIsFailedQC(true);
	writer.SaveAlignment(al);

	// if(al.IsPaired() && 
	//    al2Null ){
	//     al2=al;
	//     al2Null=false;
	//     continue;
	// }else{
	//     if(al.IsPaired() && 
	//        !al2Null){
	// 	if(al.Name != al2.Name ){
	// 	    cerr << "filterHighEditDistance: Seq#1 "<<al.Name<<" has a different id than seq #2 "<<al2.Name<<" exiting " << endl;
	// 	    return 1;
	// 	} 
	// 	total++;


	// 	if(passed(&al) &&  passed(&al2) ){
	// 	    writer.SaveAlignment(al2);
	// 	    writer.SaveAlignment(al);
	// 	}
	
		
	//     }else{  //  SINGLE END
	// 	total++;
	// 	if(passed(&al)  )
	// 	    writer.SaveAlignment(al);
	//     } //end single end

	//     al2Null=true;
	// }
		    
    } //while al

    reader.Close();
    writer.Close();

    cerr<<"filterHighEditDistance: out of "<<total<<" sequences we failed "<<discardedLength<<" due to length and "<<discardedEdit<<" due to edit distance"<<endl;
    return 0;
}

