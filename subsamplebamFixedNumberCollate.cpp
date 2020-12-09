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
			      "This program subsamples (without duplicates) a BAM file which has been processed through samtools collate to have an exact number of sequences\n"+
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

    bool sawFirst=false;
    while ( reader.GetNextAlignment(al) ) {
	if( ((al.AlignmentFlag & 0x0400) != 0) || //PCR optical dup
	    ((al.AlignmentFlag & 0x0800) != 0)){    //supp. alignment
	    //cerr<<"sec "<<al.Name<<" "<<al.AlignmentFlag<<endl;
	    continue; }//skip secondary

	if(totalSeq==UINT_MAX){
	    cerr<<"Errror: your input BAM file contains too many sequences limit: "<<UINT_MAX<<endl;
	}

	if( al.IsPaired() ){
	    if(al.IsFirstMate()){
		flagUse->push_back(false);
		sawFirst=true;
		totalSeq++;
	    }else{
		if(!sawFirst){
		    cerr<<"ERROR: read"<<al.Name<<" is a second mate but was not preceded by a first mate, did you use \"samtools collate\"?"<<endl;
		    return 1;
		}
	    }
	}else{//single
	    flagUse->push_back(false);
	    sawFirst=false;
	    totalSeq++;	    
	}
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

    //BamAlignment al2;
    string al2n;
    sawFirst=false;
    bool flushSecond=false;

    while ( reader.GetNextAlignment(al) ) {
	if( ((al.AlignmentFlag & 0x0400) != 0) || //PCR optical dup
	    ((al.AlignmentFlag & 0x0800) != 0)){    //supp. alignment
	    //cerr<<"sec "<<al.Name<<endl;
	    continue; }//skip secondary

	if( al.IsPaired() ){
	    if(al.IsFirstMate()){
		//flagUse->push_back(false);
		if( flagUse->at(countSeq) ){
		    writer.SaveAlignment(al);	    
		    flushSecond=true;
		    al2n    =al.Name;
		}
		
		sawFirst=true;	     
		countSeq++;
	    }else{
		if(!sawFirst){  cerr<<"ERROR: read "<<al.Name<<" is a second mate but was not preceded by a first mate, did you use \"samtools collate\"?"<<endl;   return 1; }

		if(flushSecond){
		    if(al2n != al.Name){cerr<<"ERROR: read name "<<al.Name<<" does not match the previous "<<al2n<<" "<<endl;   return 1; }
		    
		    writer.SaveAlignment(al);	    
		    flushSecond=false;
		}
	    }
	}else{//single
	    if( flagUse->at(countSeq) )
		writer.SaveAlignment(al);	    
	    sawFirst   =false;
	    flushSecond=false;
	    countSeq++;
	}

    }
    

    reader.Close();
    writer.Close();

    cerr<<"Wrote succesfully "<<numberToSubsample<<" reads out of "<<totalSeq<<" sequences"<<endl;

    return 0;
}



