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
	cout<<"To return the RG for each record in a BAM file"<<endl<<"Usage:"<<endl;
	cout<<""<<endl;
	cout<<"retrieveRG [options] in.bam"<<endl;
	cout<<" Options:"<<endl;
	cout<<" "<<endl;
	cout<<" -m : Use only mapped reads"<<endl;
	cout<<" -u : Use only unmapped reads"<<endl;
	cout<<endl;
	return 1;
    }


    string bamfiletopen;
    bool useOnlyMapped=false;
    bool useOnlyUnmapped=false;
    for(int i=1;i<argc;i++){
	if(strcmp(argv[i],"-m")==0)
	    useOnlyMapped   = true;	
	if(strcmp(argv[i],"-u")==0)
	    useOnlyUnmapped = true;	       
    }
    bamfiletopen=string(argv[argc-1]);

    BamReader reader;

    if(useOnlyUnmapped && useOnlyMapped){
	cout<<"Cannot specify both unmapped and mapped reads"<< endl;
    	return 1;
    }
    if ( !reader.Open(bamfiletopen) ) {
    	cerr << "Could not open input BAM file: "<<    bamfiletopen << endl;
    	return 1;
    }

    if ( !reader.LocateIndex() ){
	cerr << "warning: cannot locate index for file " << bamfiletopen<<endl;
	//return 1;
    }


    BamAlignment al;
    while ( reader.GetNextAlignment(al) ) {	
	if( useOnlyMapped){
	    if(al.IsMapped()){
		string rgTag;
		al.GetTag("RG",rgTag);
		cout<<rgTag<<endl;
	    }
	}else{
	    if( useOnlyUnmapped){
		if(!al.IsMapped()){
		    string rgTag;
		    al.GetTag("RG",rgTag);
		    cout<<rgTag<<endl;
		}
	    }else{
		string rgTag;
		al.GetTag("RG",rgTag);
		cout<<rgTag<<endl;
	    }
	}
    }

    reader.Close();

    return 0;
}
