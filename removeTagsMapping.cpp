/*
 * removeTagsMapping
 * Date: Apr-19-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */


#include <iostream>

#include "api/BamMultiReader.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"
#include "libgab.h"
#include "PutProgramInHeader.h"

using namespace std;
using namespace BamTools;

const uint32_t flagSingleReads =  4; // 00000100
const uint32_t flagFirstPair   = 77; // 01001101
const uint32_t flagSecondPair  =141; // 10001101


int main (int argc, char *argv[]) {

    bool revc1=false;
    bool revc2=false;

    if(argc < 3){
	cerr<<"This program strips the mapping information and tags"<<endl;
	cerr<<"Usage "<<argv[0]<<" [bam file in] [bam file out]"<<endl;
	cerr<<"\tOptions:\n"<<endl;
	cerr<<"\t\t--revc1\t\t\treverse complement the reverse reads, ignore flags (Default:"<<booleanAsString(revc1)<<") "<<endl;
	cerr<<"\t\t--revc2\t\t\treverse complement the reverse reads, ignore flags (Default:"<<booleanAsString(revc2)<<") "<<endl;
	return 1;
    }


    for(int i=1;i<(argc-2);i++){ //all but the last arg
	
	if(string(argv[i]) == "--revc1"){
	    revc1=true;
	    continue;
	}

	if(string(argv[i]) == "--revc2"){
	    revc2=true;
	    continue;
	}
       
	cerr<<"Unknown option "<<argv[i] <<" exiting"<<endl;
	return 1;
    }


    string bamfiletopen  = string(argv[argc-2]);
    string bamfiletwrite = string(argv[argc-1]);

    cerr<<"Reading "<<bamfiletopen<<" writing to "<<bamfiletwrite<<endl;

    BamReader reader;
    BamWriter writer;

    if ( !reader.Open(bamfiletopen) ) {
    	cerr << "Could not open input BAM files." << endl;
    	return 1;
    }

    SamHeader  myHeader=reader.GetHeader();
    SamProgram sp;

    string pID          = "removeTagsMapping";   
    string pName        = "removeTagsMapping";   
    string pCommandLine = "";
    for(int i=0;i<(argc);i++){
	pCommandLine += (string(argv[i])+" ");
    }
    putProgramInHeader(&myHeader,pID,pName,pCommandLine,returnGitHubVersion(string(argv[0]),"."));

    //no @SQ
    myHeader.Sequences.Clear();
    vector< RefData > 	emptyRefVector;

    if( !writer.Open(bamfiletwrite,myHeader,emptyRefVector ) ) {
    	cerr << "Could not open output BAM file  "<<bamfiletwrite << endl;
    	return 1;	
    }

    
    BamAlignment al;
    while ( reader.GetNextAlignment(al) ) {
	//deleting tag data
	al.TagData="";

	//reset the flag
	if(al.IsPaired()){
	    if(al.IsFirstMate()){

		if(revc1){
		    al.QueryBases = reverseComplement(al.QueryBases);
		    string st     = al.Qualities;
		    reverse(st.begin(), st.end()); 
		    al.Qualities    = st;
		}else{
		    if(al.IsReverseStrand()){//reverse the first mate
			al.QueryBases = reverseComplement(al.QueryBases);
			string st     = al.Qualities;
			reverse(st.begin(), st.end()); 
			al.Qualities    = st;
		    }
		}
		al.AlignmentFlag =  flagFirstPair;
	    }else{

		if(revc2){
		    al.QueryBases = reverseComplement(al.QueryBases);
		    string st     = al.Qualities;
		    reverse(st.begin(), st.end()); 
		    al.Qualities    = st;
		}else{
		    //cerr<<al.Name<<" "<<al.AlignmentFlag<<" "<<al.IsReverseStrand()<<endl;
		    if(!al.IsReverseStrand()){//reverse the second mate if maps to the first strand
			al.QueryBases = reverseComplement(al.QueryBases);
			string st     = al.Qualities;
			reverse(st.begin(), st.end()); 
			al.Qualities    = st;
		    }
		}
		al.AlignmentFlag =  flagSecondPair;
	    }
	}else{


	    if(al.IsReverseStrand()){//reverse single-end
		al.QueryBases = reverseComplement(al.QueryBases);
		string st     = al.Qualities;
		reverse(st.begin(), st.end()); 
		al.Qualities    = st;
	    }
	    al.AlignmentFlag =  flagSingleReads;
	}

	//no ref or positon
	al.RefID=-1;
	al.MateRefID=-1;
	al.Position=-1;
	al.MatePosition=-1;
	//no insert size
	al.InsertSize=0;
	//no cigar
	al.CigarData.clear();
	//no mapping quality 
	al.MapQuality=0;

	writer.SaveAlignment(al);
    }

    reader.Close();
    writer.Close();

    cerr<<"Program "<<argv[0]<<" terminated gracefully"<<endl;

    return 0;
}
