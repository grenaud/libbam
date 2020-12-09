/*
 * addRGinHeaderHack
 * Date: Feb-27-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */


#include <iostream>
#include <set>

#include "api/BamMultiReader.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"
#include "libgab.h"

using namespace std;
using namespace BamTools;


int main (int argc, char *argv[]) {

    if( (argc <= 2) ||
    	(argc== 2 && string(argv[1]) == "-h") ||
    	(argc== 2 && string(argv[1]) == "-help") ||
    	(argc== 2 && string(argv[1]) == "--help") ){
	cout<<"Usage:"<<argv[0]<<" [in bam] [outbam] [sample]"<<endl<<"This program looks at the BAM file, remembers the RG, writes them in the header\nThe sample name is optional but will use the RG::ID as the sample name\nThe input CANNOT be a file descriptor because it is read twice"<<endl;
	return 1;
    }

    string bamfiletopen = string(argv[1]);
    
    string sampleName;
    if(argc == 4){
    	sampleName = string(argv[3]);
    }

    set<string> allRG;
    BamReader  reader;

    if ( !reader.Open(bamfiletopen) ) {
    	cerr << "Could not open input BAM file:" << bamfiletopen<<endl;
    	return 1;
    }
    
    SamHeader myHeader = reader.GetHeader();
    const RefVector references = reader.GetReferenceData();

    BamAlignment al;
    unsigned totalRead   =0;
    unsigned totalWritten=0;

    while ( reader.GetNextAlignment(al) ) {

	if(al.HasTag("RG")){
	    //al.RemoveTag("RG");
	    string rgTag;
	    al.GetTag("RG",rgTag);
	    //cout<<rgTag<<endl;
	    allRG.insert(rgTag);
	}


	//writer.SaveAlignment(al);
	totalRead++;
    } //while al

    reader.Close();

    SamReadGroupDictionary  srgd;
    set<string>::iterator itRG;                                                                                                                                                   
    for ( itRG=allRG.begin(); itRG != allRG.end(); itRG++ ){                                                                                                         
	//cout<<*itRG <<endl;
        SamReadGroup srg ( *itRG );
        if( sampleName.empty() ){
	    srg.Sample   =(*itRG);
	}else{
	    srg.Sample   = sampleName;
	}
	srg.SequencingTechnology = "ILLUMINA";
        srgd.Add( srg );  
    }

    myHeader.ReadGroups=srgd;

    string bamFileOUT   = string(argv[2]);
    BamWriter writer;

    if ( !writer.Open(bamFileOUT,myHeader,references) ) {
    	cerr << "Could not open output BAM file:"<<bamFileOUT << endl;
    	return 1;
    }

    if ( !reader.Open(bamfiletopen) ) {
    	cerr << "Could not open input BAM file:" << bamfiletopen<<endl;
    	return 1;
    }

    while ( reader.GetNextAlignment(al) ) {
	writer.SaveAlignment(al);
	totalWritten++;
    } 
    reader.Close();
    writer.Close();

    if(totalWritten != totalRead){
	cerr << "ERROR, read "<<totalRead<<" but wrote "<<totalWritten<<" reads"<< endl;
    	return 1;
    }else{
	cerr << "Success "<<totalRead<<" and wrote "<<totalWritten<<" reads"<< endl;
    }

    return 0;
}

