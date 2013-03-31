#include <iostream>
#include <cstring>
#include <string>

#include "api/BamMultiReader.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"

using namespace std;
using namespace BamTools;

int main (int argc, char *argv[]) {

    if( (argc== 1) ||
	(argc== 2 && string(argv[1]) == "-h") ||
	(argc== 2 && string(argv[1]) == "-help") ||
	(argc== 2 && string(argv[1]) == "--help") ){
	cout<<"To check if RG changes between two BAM file (except unknown,conflict,wrong)"<<endl<<"Usage:"<<endl;
	cout<<""<<endl;
	cout<<"retrieveRG old.bam new.bam"<<endl;
	cout<<endl;
	return 1;
    }


    string bamfiletopen1=string(argv[1]);
    string bamfiletopen2=string(argv[2]);

    BamReader reader1;
    BamReader reader2;

    if ( !reader1.Open(bamfiletopen1) ) {
    	cerr << "Could not open input BAM file: "<<    bamfiletopen1 << endl;
    	return 1;
    }

    if ( !reader2.Open(bamfiletopen2) ) {
    	cerr << "Could not open input BAM file: "<<    bamfiletopen2 << endl;
    	return 1;
    }

    
    BamAlignment al1;
    BamAlignment al2;

    while ( reader1.GetNextAlignment(al1) ) {	
	reader2.GetNextAlignment(al2);

	if(al1.Name != al2.Name){
	    cerr<<"reads differ in name "<<al1.Name<<"\t"<<al2.Name << endl;
	    return 1;
	}

	
	string rgTag1;
	if(!al1.GetTag("RG",rgTag1)){
	    cerr << "Could not get tag from: "<<    al1.Name << endl;
	    return 1;
	}
	string rgTag2;
	if(!al2.GetTag("RG",rgTag2)){
	    cerr << "Could not get tag from: "<<    al2.Name << endl;
	    return 1;
	}
	
	if(rgTag1 == "unknown"  ||
	   rgTag1 == "conflict" ||
	   rgTag1 == "wrong"    ||
	   rgTag2 == "unknown"  ||
	   rgTag2 == "conflict" ||
	   rgTag2 == "wrong"    ){
	    //cerr<<"rgTag1 "<<rgTag1<<" rgTag2 "<<rgTag2<<endl;
	    //skip
	    continue;
	}else{
	    if(rgTag1 != rgTag2){
		cerr<<"reads  "<<al1.Name<<"\t"<<al2.Name << "changed read group from "<<rgTag1<<" to "<<rgTag2<<endl;		
	    }else{
		//cerr<<"rgTag1 "<<rgTag1<<" rgTag2 "<<rgTag2<<endl;
	    }
	}
	   
    }

    reader1.Close();
    reader2.Close();

    return 0;
}
