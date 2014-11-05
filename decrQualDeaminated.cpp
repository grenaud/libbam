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
    int bpToDecrease5=1;
    int bpToDecrease3=2;

    const string usage=string(string(argv[0])+" [options] input.bam out.bam"+"\n\n"+
			      "\tThis program takes a BAM file as input and produces\n"+
			      "\tanother where the putative deaminated bases have\n"+
			      "\ta base quality score of "+intStringify(baseQualForDeam)+"\n"+
			      "\tgiven an "+intStringify(offset)+" offset \n"+
			      "\n"+
			      "\tOptions:\n"+
			      "\t\t"+"-n5" +"\t\t\t"+"Decrease the nth bases surrounding the 5' ends (Default:"+stringify(bpToDecrease5)+") "+"\n"+
			      "\t\t"+"-n3" +"\t\t\t"+"Decrease the nth bases surrounding the 3' ends (Default:"+stringify(bpToDecrease3)+") "+"\n"
			      );
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

    for(int i=1;i<(argc-2);i++){ //all but the last arg

    	if( string(argv[i]) == "-n5" ){
	    bpToDecrease5 = destringify<int>(argv[i+1]);
            i++;
    	    continue;
    	}

    	if( string(argv[i]) == "-n3" ){
	    bpToDecrease3 = destringify<int>(argv[i+1]);
            i++;
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

			//5' of first mate reversed
			indexToCheck=al.QueryBases.length()-1;
			for(int i=0;i<bpToDecrease5;i++){			
			    if(toupper(al.QueryBases[indexToCheck]) == 'A'){
				al.Qualities[indexToCheck]=char(offset+baseQualForDeam);
			    }
			    indexToCheck=max(indexToCheck-1,0);
			}

			
		    }else{
			int indexToCheck;
			//5' of first mate
			indexToCheck=0;
			for(int i=0;i<bpToDecrease5;i++){ //first base			
			    if(toupper(al.QueryBases[indexToCheck]) == 'T'){
				al.Qualities[indexToCheck]=char(offset+baseQualForDeam);
			    }
			    indexToCheck=min(indexToCheck+1,int(al.Qualities.size()));
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

			    //3' of second mate reversed
			    indexToCheck=al.QueryBases.length()-1;
			    for(int i=0;i<bpToDecrease3;i++){			
				if(toupper(al.QueryBases[indexToCheck]) == 'T'){
				    al.Qualities[indexToCheck]=char(offset+baseQualForDeam);
				}
				indexToCheck=max(indexToCheck-1,0);
			    }
			    

			}else{
			    int indexToCheck;
			    
			    //3' of second mate forward
			    indexToCheck=0;
			    for(int i=0;i<bpToDecrease3;i++){ //first base			
				if(toupper(al.QueryBases[indexToCheck]) == 'A'){
				    al.Qualities[indexToCheck]=char(offset+baseQualForDeam);
				}
				indexToCheck=min(indexToCheck+1,int(al.Qualities.size()));
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

		    //5' of single read reversed
		    indexToCheck=al.QueryBases.length()-1;
		    for(int i=0;i<bpToDecrease5;i++){			
			if(toupper(al.QueryBases[indexToCheck]) == 'A'){
			    al.Qualities[indexToCheck]=char(offset+baseQualForDeam);
			}
			indexToCheck=max(indexToCheck-1,0);
		    }

		    //3' of single read reversed
		    indexToCheck=0;
		    for(int i=0;i<bpToDecrease3;i++){ //first base			
			if(toupper(al.QueryBases[indexToCheck]) == 'A'){
			    al.Qualities[indexToCheck]=char(offset+baseQualForDeam);
			}
			indexToCheck=min(indexToCheck+1,int(al.Qualities.size()));
		    }

		    

		}else{

		    int indexToCheck;
		    
		    //5' of single read
		    indexToCheck=0;
		    for(int i=0;i<bpToDecrease5;i++){ //first base			
			if(toupper(al.QueryBases[indexToCheck]) == 'T'){
			    al.Qualities[indexToCheck]=char(offset+baseQualForDeam);
			}
			indexToCheck=min(indexToCheck+1,int(al.Qualities.size()));
		    }


		    //3' of single read
		    indexToCheck=al.QueryBases.length()-1;
		    for(int i=0;i<bpToDecrease3;i++){			
			if(toupper(al.QueryBases[indexToCheck]) == 'T'){
			    al.Qualities[indexToCheck]=char(offset+baseQualForDeam);
			}
			indexToCheck=max(indexToCheck-1,0);
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

