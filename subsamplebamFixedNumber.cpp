#include <iostream>
#include <vector>
#include <climits> 

#include "api/BamMultiReader.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"

#include "libgab.h"

using namespace std;
using namespace BamTools;

int main (int argc, char *argv[]) {
    const string usage=string("\t"+string(argv[0])+
			      " [options] [BAM file in] [BAM file out] [number to subsample] "+"\n\n"+
			      "This program subsamples (without duplicates) a BAM file to have an exact number of sequences\n"+
			      "It does NOT shuffle the sequences"+
			      "Options:\n"+
			      "none");
			      
			      
    if( (argc== 1) ||
	(argc== 2 && string(argv[1]) == "-h") ||
	(argc== 2 && string(argv[1]) == "-help") ||
	(argc== 2 && string(argv[1]) == "--help") ){
	cerr<<usage<<endl;
	cerr<<endl;
	return 1;
    }


    string       bamfiletopen                                = string(argv[argc-3]);
    string       outfile                                     = string(argv[argc-2]);
    unsigned int numberToSubsample = destringify<unsigned int>(string(argv[argc-1]));

    unsigned int totalSeq=0;
    unsigned int countSeq=0;
    vector<bool> * flagUse=new vector<bool>();
    


    BamReader reader;

    if ( !reader.Open(bamfiletopen) ) {
    	cerr << "Could not open input BAM files." << endl;
    	return 1;
    }
    cerr<<"Reading the bam file for the first time ... "<<endl;
    BamAlignment al;
    while ( reader.GetNextAlignment(al) ) {
	flagUse->push_back(false);
	if(totalSeq==UINT_MAX){
	    cerr<<"Errror: your input BAM file contains too many sequences limit: "<<UINT_MAX<<endl;
	}

	//TODO check for unsigned int limit
	totalSeq++;
    }


    if(totalSeq<=numberToSubsample){
	cerr<<"Cannot subsample "<<numberToSubsample<<" from a BAM file with "<<totalSeq<<" sequences"<<endl;
	return 1;
    }
    cerr<<"done!"<<endl;

    cerr<<"Selecting random reads ... "<<endl;

    unsigned int numberFlagged=0;
    while(numberFlagged!=numberToSubsample){
	unsigned int randomUnsInt=randomUint();
	if(randomUnsInt>=totalSeq)
	    continue;

	if(! flagUse->at(randomUnsInt) ){
	    flagUse->at(randomUnsInt)=true;
	    numberFlagged++;
	}else{
	    continue;
	}
    }
    cerr<<"done!"<<endl;



    cerr<<"Writing to output ... "<<endl;
    
    //re-open
    BamWriter writer;
    if( !writer.Open(outfile,
		     reader.GetHeader(),
		     reader.GetReferenceData() 
		     ) ) {
    	cerr << "Could not open output BAM file  "<<outfile << endl;
    	return 1;	
    }

    //re-open
    if ( !reader.Open(bamfiletopen) ) {
    	cerr << "Could not open input BAM files." << endl;
    	return 1;
    }

    while ( reader.GetNextAlignment(al) ) {
	if( flagUse->at(countSeq) )
	    writer.SaveAlignment(al);
	countSeq++;
    }
    

    reader.Close();
    writer.Close();

    cerr<<"Wrote succesfully "<<numberToSubsample<<" reads out of "<<totalSeq<<" sequences"<<endl;

    return 0;
}



