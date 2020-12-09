/*
 * removeTagsMapping
 * Date: Apr-19-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */


#include <iostream>
#include <map>

#include "api/BamMultiReader.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"
#include "libgab.h"
#include "PutProgramInHeader.h"

using namespace std;
using namespace BamTools;

int main (int argc, char *argv[]) {

    if(argc != 4){
	cerr<<"This program reads a list of readnames and an input bam file and prints the reads matching those"<<endl;
	cerr<<"Usage "<<argv[0]<<" [list names] [bam file in] [bam file out]"<<endl;
	return 1;
    }

    string namefile      = string(argv[1]);
    string bamfiletopen  = string(argv[2]);
    string bamfiletwrite = string(argv[3]);
    map<string,bool> names2found;
    igzstream nameFile;
    string line;
    nameFile.open(namefile.c_str(), ios::in);
    
    if (nameFile.good()){
	while ( getline (nameFile,line)){
	    vector<string> token=allTokens(line,' ');
	    if(token.size() != 1){
		cerr << "Line does not have a single field "<<line<<endl;
		return 1;
	    }
	    //cout<<"#"<<line<<"#"<<endl;
	    names2found[line]=false;
	}
	nameFile.close();
    }else{
	cerr << "Unable to open file "<<namefile<<endl;
	return 1;
    }


    cerr<<"Reading "<<bamfiletopen<<" writing to "<<bamfiletwrite<<endl;

    BamReader reader;
    BamWriter writer;

    if ( !reader.Open(bamfiletopen) ) {
    	cerr << "Could not open input BAM files." << endl;
    	return 1;
    }

    SamHeader  myHeader=reader.GetHeader();
    const RefVector references = reader.GetReferenceData();

    SamProgram sp;

    string pID          = "retrieveReadsWithName";   
    string pName        = "retrieveReadsWithName";   
    string pCommandLine = "";
    for(int i=0;i<(argc);i++){
	pCommandLine += (string(argv[i])+" ");
    }
    putProgramInHeader(&myHeader,pID,pName,pCommandLine,returnGitHubVersion(string(argv[0]),"."));


    if( !writer.Open(bamfiletwrite,myHeader,references ) ) {
    	cerr << "Could not open output BAM file  "<<bamfiletwrite << endl;
    	return 1;	
    }

    
    BamAlignment al;
    while ( reader.GetNextAlignment(al) ) {
	map<string,bool>::iterator it=names2found.find(al.Name);

	if(it != names2found.end()){
	    //names2found.find(al.Name)
	    // cout<<"    found"<<al.Name<<endl;
	    it->second=true;
	    writer.SaveAlignment(al);
	}else{
	    //cout<<
	    // cout<<"not found"<<al.Name<<endl;
	}


    }

    reader.Close();
    writer.Close();

    cerr<<"Program "<<argv[0]<<" terminated gracefully"<<endl;

    vector<string> notfound;
    for(map<string,bool>::iterator it = names2found.begin();it != names2found.end(); ++it) {
	if(!it->second)
	    notfound.push_back(it->first);
    }


    for(unsigned int i=0;i<notfound.size();i++){
	cerr<<"Name "<<notfound[i]<<" not found"<<endl;
    }

    return 0;
}

