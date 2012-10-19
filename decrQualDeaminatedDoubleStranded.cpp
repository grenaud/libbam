#include <iostream>
#include <string>
#include <cstring>

#include "api/BamMultiReader.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"

using namespace std;
using namespace BamTools;

string intStringify(int i) {
    stringstream s;
    s << i;
    return s.str();
}

const int baseQualForDeam=2;
const int offset=33;



int main (int argc, char *argv[]) {

    bool mapped  =false;
    bool unmapped=false;

    const string usage=string(string(argv[0])+" [options] input.bam out.bam"+"\n\n"+
			      "This program takes a BAM file as input and produces\n"+
			      "another where the putative deaminated bases have\n"+
			      "a base quality score of "+intStringify(baseQualForDeam)+"\n"+
			      "given an "+intStringify(offset)+" offset \n"+
			      "\n"+
			      "Options:\n");
			      // "\t"+"-u , --unmapped" +"\n\t\t"+"For an unmapped bam file"+"\n"+
			      // "\t"+"-m , --mapped"   +"\n\t\t"+"For an mapped bam file"+"\n");
			      
			      

    if( (argc== 1) ||
	(argc== 2 && string(argv[1]) == "-h") ||
	(argc== 2 && string(argv[1]) == "-help") ||
	(argc== 2 && string(argv[1]) == "--help") ){
	cout<<"Usage:"<<endl;
	cout<<usage<<endl;
	cout<<""<<endl;
	return 1;
    }


    if(argc != 3){
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


    vector<RefData>  testRefData=reader.GetReferenceData();
    const SamHeader header = reader.GetHeader();
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
		    
		    for(int indexToCheck=(al.QueryBases.length()-1);indexToCheck>(al.QueryBases.length()-6);indexToCheck--){
			if(toupper(al.QueryBases[indexToCheck]) == 'A'){
			    al.Qualities[indexToCheck]=char(offset+baseQualForDeam);
			}
		    }

		}else{ //if not reverse complemented, we decrease the first 5 Ts regardless of first or second mate

		    for(int indexToCheck=0;indexToCheck<5;indexToCheck++){
			if(toupper(al.QueryBases[indexToCheck]) == 'T'){
			    al.Qualities[indexToCheck]=char(offset+baseQualForDeam);
			}
		    }


		}
	
	    }//end of paired end
	    else{//we consider single reads to have been sequenced from 5' to 3'


		if(al.QueryBases.length() <= 5){
		    cerr << "We do not reads with less than 5bp" << endl;
		    return 1;
		}
		for(int indexToCheck=0;indexToCheck<5;indexToCheck++){
		    if(toupper(al.QueryBases[indexToCheck]) == 'T'){
			al.Qualities[indexToCheck]=char(offset+baseQualForDeam);
		    }
		}

		for(int indexToCheck=(al.QueryBases.length()-1);indexToCheck>(al.QueryBases.length()-6);indexToCheck--){
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

