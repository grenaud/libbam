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
			      "but miss their other pair will be marked as single end\n");
		
    cerr<<"To verify"<<endl;
    return 1;
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
    unsigned int fixed=0;

    while ( reader.GetNextAlignment(al) ) {
	total++;

	if(al.IsPaired() && 
	   al2Null ){
	    al2=al;
	    al2Null=false;
	    continue;
	}else{
	    if(al.IsPaired() && 
	       !al2Null){

		if(al.Name != al2.Name ){
		    fixed++;
		    //fix al2
		    al2.SetIsPaired(false);
		    al2.SetIsPrimaryAlignment(true);
		    al2.SetIsProperPair(false);
		    al2.SetIsMateMapped(false);

		    writer.SaveAlignment(al2);


		    //keep al for next iteration
		    al2=al;
		    al2Null=false;
		}else{ //correct pair
		    writer.SaveAlignment(al2);
		    writer.SaveAlignment(al);
		}
	
		
	    }else{  //  SINGLE END
		total++;
		if(passed(&al)  )
		    writer.SaveAlignment(al);
	    } //end single end

	    al2Null=true;
	}
		    
    } //while al

    reader.Close();
    writer.Close();

    cerr<<"filterHighEditDistance: out of "<<total<<" sequences we fixed "<<fixed<<" sequences "<<endl;
    return 0;
}

