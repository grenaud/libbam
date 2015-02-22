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
	cout<<endl<<"To print stats about the same reads mapped to more than one sample"<<endl<<endl<<"Usage:"<<endl;
	cout<<""<<endl;
	cout<<"\tcrossSamplestat file1.bam file2.bam ..."<<endl;
	cout<<endl;
	return 1;
    }
    // cout<<argc<<endl;
    // return 1;
    
    BamReader readers[argc-1];
    //open
    cout<<"ReadID\tReadSeq\tReadQuality\t";
    for(int i=1;i<(argc);i++){
	// cout<<i<<endl;
	// cout<<string(argv[i])<<endl;
	cout<<"\tLocationRef"<<i<<"\tRef"<<i<<"AlignScore";
	if ( !readers[i-1].Open( string(argv[i]) ) ) {
	    cerr << "Could not open input BAM file: "<<  string(argv[i]) << endl;
	    return 1;
	}
    }

    cout<<endl;    

    BamAlignment al0;
    BamAlignment alN;

    while ( readers[0].GetNextAlignment(al0) ) {
	cout<<al0.Name<<"\t"<<al0.QueryBases<<"\t"<<al0.Qualities;
	cout<< "\t"<<al0.Position <<"\t"<<al0.MapQuality;
	 for(int i=1;i<(argc-1);i++){
	     if(!readers[i].GetNextAlignment(alN)){
		 cerr << "Could not read alignment from : "<<  string(argv[i+1]) << endl;
		 return 1;
	     }


	     if(al0.Name != alN.Name){
		 cerr<<"reads differ in name "<<al0.Name<<"\t"<<alN.Name << endl;
		 return 1;
	     }
	     //cout<<"name "<<i<<"\t"<<alN.Name<<endl;
	     cout<< "\t"<<alN.Position <<"\t"<<alN.MapQuality;
	 }
	 cout<<endl;
    }

    //close 
    for(int i=1;i<(argc);i++){
	if ( !readers[i-1].Close(  ) ) {
	    cerr << "Could not close BAM file: "<< string(argv[i-1]) << endl;
	    return 1;
	}

    }
    

    return 0;
}
