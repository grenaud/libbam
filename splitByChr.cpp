/*
 * splitByChr
 * Date: Feb-28-2013 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */


#include <iostream>

#include "api/BamMultiReader.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"
#include "utils.h"

using namespace std;
using namespace BamTools;


int main (int argc, char *argv[]) {

     if( (argc!= 3) ||
    	(argc== 2 && string(argv[1]) == "-h") ||
    	(argc== 2 && string(argv[1]) == "-help") ||
    	(argc== 2 && string(argv[1]) == "--help") ){
	 cerr<<"Usage:splitByChr [in bam] [out prefix]"<<endl<<"this program creates one bam file per chr in the with the outprefix\nFor example splitByChr in.bam out will create\nout.chr1.bam\nout.chr2.bam\n"<<endl;
    	return 1;
    }


     string bamfiletopen = string(argv[1]);
     // if(!strEndsWith(bamfiletopen,".bam")){

     // }
     string bamDirOutPrefix    = string(argv[2]);
     map<string,BamWriter *> chr2BamWriter;
     
     // if(!isDirectory(bamDirOut)){
     // 	 cerr<<"ERROR: the out directory does not exist"<<endl;
     // 	return 1;
     // }

     BamReader reader;

     if ( !reader.Open(bamfiletopen) ) {
    	cerr << "Could not open input BAM files." << endl;
    	return 1;
     }
    const SamHeader header = reader.GetHeader();
    const RefVector references = reader.GetReferenceData();
    vector<RefData>  refData=reader.GetReferenceData();


    BamWriter unmapped;

    // cout<<header.ToString()<<endl;
    // return 1;
    if ( !unmapped.Open(bamDirOutPrefix+".unmapped.bam",header,references) ) {
    	cerr << "Could not open output BAM file "<< bamDirOutPrefix+".unmapped.bam" << endl;
    	return 1;
    }



    BamAlignment al;
    unsigned int total=0;
    while ( reader.GetNextAlignment(al) ) {

	// al.SetIsFailedQC(false);
	// writer.SaveAlignment(al);
	if(al.IsMapped () ){
	    if(chr2BamWriter.find(refData[al.RefID].RefName) == chr2BamWriter.end()){ //new
		chr2BamWriter[refData[al.RefID].RefName] = new  BamWriter();
		if ( !chr2BamWriter[refData[al.RefID].RefName]->Open(bamDirOutPrefix+"."+refData[al.RefID].RefName+".bam",header,references) ) {
		    cerr     << "Could not open output BAM file "<< bamDirOutPrefix<<"."<<refData[al.RefID].RefName<<".bam" << endl;
		    return 1;
		}
		chr2BamWriter[refData[al.RefID].RefName]->SaveAlignment(al);
	    }else{
		chr2BamWriter[refData[al.RefID].RefName]->SaveAlignment(al);
	    }
	}else{
	    unmapped.SaveAlignment(al);
	}
		    
	total++;
    } //while al

    reader.Close();
    // writer.Close();
    
    unmapped.Close();

    map<string,BamWriter *>::iterator chr2BamWriterIt;
    for (chr2BamWriterIt =chr2BamWriter.begin(); 
	 chr2BamWriterIt!=chr2BamWriter.end(); 
	 chr2BamWriterIt++){
	chr2BamWriterIt->second->Close();
    }
    cerr<<"Wrote succesfully "<<total<<" reads"<<endl;


    return 0;
}

