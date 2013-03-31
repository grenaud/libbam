#include <iostream>
#include <string>
#include <cstring>

#include "api/BamMultiReader.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"
#include "PutProgramInHeader.h"

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

    // bool mapped  =false;
    // bool unmapped=false;

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

    // for(int i=1;i<(argc-1);i++){ //all but the last arg

    // 	if(strcmp(argv[i],"-m") == 0 || strcmp(argv[i],"--mapped") == 0 ){
    // 	    mapped=true;
    // 	    continue;
    // 	}

    // 	if(strcmp(argv[i],"-u") == 0 || strcmp(argv[i],"--unmapped") == 0 ){
    // 	    unmapped=true;
    // 	    continue;
    // 	}
       
    // 	cerr<<"Unknown option "<<argv[i] <<" exiting"<<endl;
    // 	return 1;
    // }

    if(argc != 3){
	cerr<<"Error: Must specify the input and output BAM files";
	return 1;
    }

    string inbamFile =argv[argc-2];
    string outbamFile=argv[argc-1];

    if(inbamFile == outbamFile){
	cerr<<"Input and output files are the same"<<endl;
	return 1;    
    }
    // if(!mapped && !unmapped){
    // 	cerr << "Please specify whether you reads are mapped or unmapped" << endl;
    // 	return 1;
    // }
    // if(mapped && unmapped){
    // 	cerr << "Please specify either mapped or unmapped but not both" << endl;
    // 	return 1;
    // }

    BamReader reader;

    if ( !reader.Open(inbamFile) ) {
    	cerr << "Could not open input BAM files." << endl;
    	return 1;
    }


    vector<RefData>  testRefData=reader.GetReferenceData();
    SamHeader header = reader.GetHeader();
    string pID          = "decrQualDeaminated";   
    string pName        = "decrQualDeaminated";   
    string pCommandLine = "";
    for(int i=0;i<(argc);i++){
        pCommandLine += (string(argv[i])+" ");
    }
    putProgramInHeader(&header,pID,pName,pCommandLine);


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

		if(al.IsFirstMate() ){ //5' end, need to check first base only
		    if(al.IsReverseStrand()){ //
			if(!al.IsMapped()){
			    cerr << "Cannot have reverse complemented unmapped reads :" <<al.Name<< endl;
			    //return 1;
			}
			int indexToCheck;
			//last
			indexToCheck=al.QueryBases.length()-1;
			if(toupper(al.QueryBases[indexToCheck]) == 'A'){
			    al.Qualities[indexToCheck]=char(offset+baseQualForDeam);
			}
		    }else{
			int indexToCheck;
			//first base
			indexToCheck=0;
			if(toupper(al.QueryBases[indexToCheck]) == 'T'){
			    al.Qualities[indexToCheck]=char(offset+baseQualForDeam);
			}
		    }
		}else{ //3' end, need to check last two bases only
		    if( al.IsSecondMate() ){
			if(al.IsReverseStrand()){ //
			    if(!al.IsMapped()){
				cerr << "Cannot have reverse complemented unmapped reads :" <<al.Name<< endl;
				//return 1;
			    }
			    int indexToCheck;

			    //second to last
			    indexToCheck=al.QueryBases.length()-2;
			    if(toupper(al.QueryBases[indexToCheck]) == 'T'){
				al.Qualities[indexToCheck]=char(offset+baseQualForDeam);
			    }

			    //last
			    indexToCheck=al.QueryBases.length()-1;
			    if(toupper(al.QueryBases[indexToCheck]) == 'T'){
				al.Qualities[indexToCheck]=char(offset+baseQualForDeam);
			    }

			}else{
			    int indexToCheck;
			    //first base
			    indexToCheck=0;
			    if(toupper(al.QueryBases[indexToCheck]) == 'A'){
				al.Qualities[indexToCheck]=char(offset+baseQualForDeam);
			    }

			    //second base
			    indexToCheck=1;
			    if(toupper(al.QueryBases[indexToCheck]) == 'A'){
				al.Qualities[indexToCheck]=char(offset+baseQualForDeam);
			    }
			}
		    }else{
			cerr << "Wrong state" << endl;
			return 1;
		    }
		}

	    }//end of paired end
	    else{//we consider single reads to have been sequenced from 5' to 3'

		if(al.IsReverseStrand()){ //need to consider 
		    if(!al.IsMapped()){
			cerr << "Cannot have reverse complemented unmapped reads :" <<al.Name<< endl;
			//return 1;
		    }

		    int indexToCheck;
		    //first base
		    indexToCheck=0;
		    if(toupper(al.QueryBases[indexToCheck]) == 'A'){
			al.Qualities[indexToCheck]=char(offset+baseQualForDeam);
		    }

		    //second base
		    indexToCheck=1;
		    if(toupper(al.QueryBases[indexToCheck]) == 'A'){
			al.Qualities[indexToCheck]=char(offset+baseQualForDeam);
		    }

		    //last
		    indexToCheck=al.QueryBases.length()-1;
		    if(toupper(al.QueryBases[indexToCheck]) == 'A'){
			al.Qualities[indexToCheck]=char(offset+baseQualForDeam);
		    }
		}else{

		    int indexToCheck;
		    //first base
		    indexToCheck=0;
		    if(toupper(al.QueryBases[indexToCheck]) == 'T'){
			al.Qualities[indexToCheck]=char(offset+baseQualForDeam);
		    }

		    //second to last
		    indexToCheck=al.QueryBases.length()-2;
		    if(toupper(al.QueryBases[indexToCheck]) == 'T'){
			al.Qualities[indexToCheck]=char(offset+baseQualForDeam);
		    }

		    //last
		    indexToCheck=al.QueryBases.length()-1;
		    if(toupper(al.QueryBases[indexToCheck]) == 'T'){
			al.Qualities[indexToCheck]=char(offset+baseQualForDeam);
		    }

		}
	    }//end of single end

	    writer.SaveAlignment(al);		


    }//    while ( reader.GetNextAlignment(al) ) {

    reader.Close();
    writer.Close();

    cerr<<"Program terminated gracefully"<<endl;   
    return 0;
}

