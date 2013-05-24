#include <iostream>
#include <string>
#include <cstring>

#include "api/BamMultiReader.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"

#include "utils.h"
#include "PutProgramInHeader.h"

using namespace std;
using namespace BamTools;


const int baseQualForDeam=2;
const int offset=33;



int main (int argc, char *argv[]) {

    // bool mapped  =false;
    // bool unmapped=false;
    int bpToDecrease=5;

    const string usage=string(string(argv[0])+" [options] input.bam out.bam"+"\n\n"+
			      "This program takes a BAM file as input and produces\n"+
			      "another where the putative deaminated bases have\n"+
			      "a base quality score of "+stringify(baseQualForDeam)+"\n"+
			      "given an "+stringify(offset)+" offset \n"+
			      "\n"+
			      "Options:\n"+
			      "\t"+"-n" +"\t\t\t"+"Decrease the nth bases surrounding the 5'/3' ends (Default:"+stringify(bpToDecrease)+") "+"\n");
			      
			      

    if( (argc== 1) ||
	(argc== 2 && string(argv[1]) == "-h") ||
	(argc== 2 && string(argv[1]) == "-help") ||
	(argc== 2 && string(argv[1]) == "--help") ){
	cout<<"Usage:"<<endl;
	cout<<usage<<endl;
	cout<<""<<endl;
	return 1;
    }

    //all but last 2
    for(int i=1;i<(argc-2);i++){ 

        if(string(argv[i]) == "-n"){
	    bpToDecrease =destringify<int>(argv[i+1]);
            i++;
            continue;
	}

    }



    if(argc < 3){
	cerr<<"Error: Must specify the input and output BAM files";
	return 1;
    }

    string inbamFile =argv[argc-2];
    string outbamFile=argv[argc-1];


    BamReader reader;

    if ( !reader.Open(inbamFile) ) {
    	cerr << "Could not open input BAM files." << endl;
    	return 1;
    }

    SamHeader header = reader.GetHeader();
    string pID          = "decrQualDeaminatedDoubleStranded";   
    string pName        = "decrQualDeaminatedDoubleStranded";   
    string pCommandLine = "";
    for(int i=0;i<(argc);i++){
        pCommandLine += (string(argv[i])+" ");
    }
    putProgramInHeader(&header,pID,pName,pCommandLine);


    vector<RefData>  testRefData=reader.GetReferenceData();
    // const SamHeader header = reader.GetHeader();
    const RefVector references = reader.GetReferenceData();

    BamWriter writer;
    if ( !writer.Open(outbamFile, header, references) ) {
	cerr << "Could not open output BAM file" << endl;
	return 1;
    }

    BamAlignment al;
    // BamAlignment al2;
    // bool al2Null=true;
    
    while ( reader.GetNextAlignment(al) ) {

	    if(al.IsPaired() ){  


		// cerr << "We do not support paired end" << endl;
		// return 1;

		if(al.IsReverseStrand()){ //is reverse complemented, we decrease the last As regardless of first or second mate
		    
		    for(int indexToCheck=(al.QueryBases.length()-1);indexToCheck>( (al.QueryBases.length()-1)-bpToDecrease);indexToCheck--){
			if(toupper(al.QueryBases[indexToCheck]) == 'A'){
			    al.Qualities[indexToCheck]=char(offset+baseQualForDeam);
			}
		    }

		}else{ //if not reverse complemented, we decrease the first 5 Ts regardless of first or second mate

		    for(int indexToCheck=0;indexToCheck<bpToDecrease;indexToCheck++){
			if(toupper(al.QueryBases[indexToCheck]) == 'T'){
			    al.Qualities[indexToCheck]=char(offset+baseQualForDeam);
			}
		    }


		}
	
	    }//end of paired end
	    else{//we consider single reads to have been sequenced from 5' to 3'


		if(al.QueryBases.length() <= 5){
		    cerr << "We do not process reads with less than 5bp" << endl;
		    return 1;
		}
		for(int indexToCheck=0;indexToCheck<bpToDecrease;indexToCheck++){
		    if(toupper(al.QueryBases[indexToCheck]) == 'T'){
			al.Qualities[indexToCheck]=char(offset+baseQualForDeam);
		    }
		}

		for(int indexToCheck=(al.QueryBases.length()-1);indexToCheck>( (al.QueryBases.length()-1) -bpToDecrease);indexToCheck--){
		    if(toupper(al.QueryBases[indexToCheck]) == 'A'){
			al.Qualities[indexToCheck]=char(offset+baseQualForDeam);
		    }
		}
		

	    }//end of single end

	    writer.SaveAlignment(al);		


    }//    while ( reader.GetNextAlignment(al) ) {

    reader.Close();
    writer.Close();
   
    return 0;
}

