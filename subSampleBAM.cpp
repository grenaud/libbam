#include <iostream>
#include <string>
#include <cstdlib>
#include <sys/time.h>

#include "api/BamMultiReader.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"

using namespace std;
using namespace BamTools;

bool   initialized       = false;

static inline double randomGen(){
    if(initialized == false){
	struct timeval time;
	gettimeofday(&time,NULL);
	srand(time.tv_usec);
	// srand((unsigned)time(0));
	initialized=true;
	return  (double(rand())/double(RAND_MAX));       
    }else{
	return  (double(rand())/double(RAND_MAX));       
    }
}

int main (int argc, char *argv[]) {

    if( (argc== 1) ||
	(argc== 2 && string(argv[1]) == "-h") ||
	(argc== 2 && string(argv[1]) == "-help") ||
	(argc== 2 && string(argv[1]) == "--help") ){
	cout<<"Usage:"<<endl;
	cout<<""<<endl;
	cout<<"subSampleBAM fraction in.bam out.bam"<<endl;
	return 1;
    }

    double  fraction      = atof  (argv[1]);
    string bamfiletopen   = string(argv[2]);
    string outputFilename = string(argv[3]);
    BamReader reader;

    if(fraction < 0.0 || fraction > 1.0){
    	cerr << "The fraction must be between 0 and 1" << endl;
    	return 1;
    }

    if ( !reader.Open(bamfiletopen) ) {
    	cerr << "Could not open input BAM files." << endl;
    	return 1;
    }



    const SamHeader header = reader.GetHeader();
    const RefVector references = reader.GetReferenceData();

    BamWriter writer;
    if ( !writer.Open(outputFilename, header, references) ) {
	cerr << "Could not open output BAM file" << endl;
	return 1;
    }

    BamAlignment al;
    while ( reader.GetNextAlignment(al) ) {
	if(al.IsPaired()){
	    if(al.IsFirstMate() && randomGen() < fraction){
		writer.SaveAlignment(al);
		reader.GetNextAlignment(al);
		writer.SaveAlignment(al);
	    }
	}else{
	    if( randomGen() < fraction){
		writer.SaveAlignment(al);
	    }
	    
	}
    }
    reader.Close();
    writer.Close();
   
    return 0;
}

