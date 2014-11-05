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


#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

int main (int argc, char *argv[]) {

    int length=10;
    bool onlyRev=false;
    bool onlyFwd=false;

    const string usage=string(string(argv[0])+" [options] input.bam out.bam"+"\n\n"+
			      "This program takes a BAM file as input and produces\n"+
			      "another where the first l bases have been kept\n"+
			      "\n"+
			      "Options:\n"+
			      "\t"+"-l" +"\t\t"+"Keep this length (Default "+stringify(length)+")"+"\n"
			      "\t"+"-r" +"\t\t"+"Only trim reverse reads (Default "+booleanAsString(onlyRev)+")"+"\n"
			      "\t"+"-f" +"\t\t"+"Only trim forward reads (Default "+booleanAsString(onlyFwd)+")"+"\n"
			      );
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

    for(int i=1;i<(argc-2);i++){ //all but the last arg

	// 	if(strcmp(argv[i],"-m") == 0 || strcmp(argv[i],"--mapped") == 0 ){
	// 	    mapped=true;
	// 	    continue;
	// 	}

    	if(string(argv[i]) == "-l" ){
	    length=destringify<int>(string(argv[i+1]));
	    i++;
	    continue;
	}

    	if(string(argv[i]) == "-r" ){
	    onlyRev=true;
	    continue;
	}

    	if(string(argv[i]) == "-f" ){
	    onlyFwd=true;
	    continue;
	}
       


	cerr<<"Unknown option "<<argv[i] <<" exiting"<<endl;
	return 1;
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


    vector<RefData>  testRefData=reader.GetReferenceData();
    SamHeader myHeader = reader.GetHeader();
    const RefVector references = reader.GetReferenceData();

    string pID          = "cutEndKeepBeginning";   
    string pName        = "cutEndKeepBeginning";   
    string pCommandLine = "";
    for(int i=0;i<(argc);i++){
	pCommandLine += (string(argv[i])+" ");
    }
    putProgramInHeader(&myHeader,pID,pName,pCommandLine,returnGitHubVersion(string(argv[0]),"."));

    BamWriter writer;
    if ( !writer.Open(outbamFile, myHeader, references) ) {
	cerr << "Could not open output BAM file" << endl;
	return 1;
    }

    BamAlignment al;
    // BamAlignment al2;
    // bool al2Null=true;
    
    while ( reader.GetNextAlignment(al) ) {
	if(onlyRev && al.IsFirstMate() ){//do not care about first mate
	    writer.SaveAlignment(al);		
	    continue;
	}

	if(onlyFwd && !al.IsFirstMate() ){
	    writer.SaveAlignment(al);		
	    continue;
	}


	int lengthAl=MIN(al.QueryBases.size(),length);
			 
	al.QueryBases = al.QueryBases.substr(0,lengthAl);
	al.Qualities  = al.Qualities.substr(0,lengthAl);

	writer.SaveAlignment(al);		

    }//    while ( reader.GetNextAlignment(al) ) {

    reader.Close();
    writer.Close();
    cerr<<"Program "<<argv[0]<<" terminated succesfully"<<endl;

    return 0;
}

