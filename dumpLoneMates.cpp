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
#include "utils.h"

using namespace std;
using namespace BamTools;



int main (int argc, char *argv[]) {

    const string usage=string(string(argv[0])+" [options] input.bam out.bam"+"\n\n"+
			      "This program takes a BAM file as input (sorted by name) and produces\n"+
			      "another as output where the reads marked as paired\n"+
			      "but miss their other pair will be dumped\n");
		
     if( (argc== 1) ||
    	(argc== 2 && string(argv[1]) == "-h") ||
    	(argc== 2 && string(argv[1]) == "-help") ||
    	(argc== 2 && string(argv[1]) == "--help") ){
	 cout<<"Usage: "<<usage<<endl;
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
    unsigned int total=0;
    unsigned int dumped=0;
    string lastName="";
    // bool firstRead=true;

    while ( reader.GetNextAlignment(al) ) {
	total++;

	// if(firstRead){
	//     lastName=al.Name;
	//     firstRead=false;
	// }else{
	//     if(lastName.compare(al.Name) > 0){ //name increases
	// 	cerr << "Read name= "<< al.Name << " and "<< lastName <<" show that the bam file is not sorted by name"<<endl;
	// 	return 1;
	//     }
	//     lastName=al.Name;
	// }

	if(al.IsPaired() && 
	   al2Null ){
	    al2=al;
	    al2Null=false;
	    continue;
	}else{
	    if(al.IsPaired() && !al2Null){

		if(al.Name != al2.Name ){ //this means that the first one is a lone mate 
		    dumped++;
		    
		    //dump al2 

		    //keep al for next iteration
		    al2=al;
		    al2Null=false;
		    continue;
		}else{ //correct pair
		    writer.SaveAlignment(al2);
		    writer.SaveAlignment(al);
		    al2Null=true;
		    continue;
		}
	       	
	    }

	    if(!al.IsPaired()){
		cerr << "Invalid state for "<< al.Name << " read should be not paired"<<endl;
	 	return 1;
	    }

	    //  SINGLE END
	    writer.SaveAlignment(al);
	    // else{  //  SINGLE END
	    // 		total++;
	    // 		writer.SaveAlignment(al);
	    // 	    } //end single end
	    
	    //	    al2Null=true;
	}
		    
    } //while al

    reader.Close();
    writer.Close();

    cerr<<"dumpLoneMates: out of "<<total<<" sequences we dumped "<<dumped<<" sequences "<<endl;
    return 0;
}

