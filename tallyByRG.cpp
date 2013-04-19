/*
 * tallyByRG
 * Date: Apr-16-2013 
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

     if( (argc!= 2) ||
    	(argc== 2 && string(argv[1]) == "-h") ||
    	(argc== 2 && string(argv[1]) == "-help") ||
    	(argc== 2 && string(argv[1]) == "--help") ){
	 cerr<<"Usage:tallyByRG [in bam]"<<endl<<"this program computes how many reads mapped\nto various chromosomes on a per RG basis"<<endl;
    	return 1;
    }


     string bamfiletopen = string(argv[1]);
     // if(!strEndsWith(bamfiletopen,".bam")){

     // }
     //     string bamDirOutPrefix    = string(argv[2]);
     map<string,vector<unsigned int> > rg2Tally;
     
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


    SamReadGroupDictionary 	srgd=header.ReadGroups;
    for(SamReadGroupConstIterator srgci=srgd.ConstBegin();
	srgci<srgd.ConstEnd();
	srgci++){
	//cout<<*srgci<<endl;
	const SamReadGroup rg = (*srgci);
	//cout<<rg.ID<<endl;
	rg2Tally[rg.ID] =  vector<unsigned int> (references.size(),0);	
    }
    //    return 1;

    //    BamWriter unmapped;

    // cout<<header.ToString()<<endl;
    // return 1;

    // if ( !unmapped.Open(bamDirOutPrefix+".unmapped.bam",header,references) ) {
    // 	cerr << "Could not open output BAM file "<< bamDirOutPrefix+".unmapped.bam" << endl;
    // 	return 1;
    // }



    BamAlignment al;
    unsigned int total=0;
    while ( reader.GetNextAlignment(al) ) {

	// al.SetIsFailedQC(false);
	// writer.SaveAlignment(al);
	// if(al.IsMapped () ){
	//     if(rg2Tally.find(refData[al.RefID].RefName) == rg2Tally.end()){ //new
	// 	rg2Tally[refData[al.RefID].RefName] = new  BamWriter();
	// 	if ( !rg2Tally[refData[al.RefID].RefName]->Open(bamDirOutPrefix+"."+refData[al.RefID].RefName+".bam",header,references) ) {
	// 	    cerr     << "Could not open output BAM file "<< bamDirOutPrefix<<"."<<refData[al.RefID].RefName<<".bam" << endl;
	// 	    return 1;
	// 	}
	
	//     }else{
	// 	rg2Tally[refData[al.RefID].RefName]->SaveAlignment(al);
	//     }
	// }else{
	//     unmapped.SaveAlignment(al);
	// }
	if(al.HasTag("RG")){
	    string rgTag;
	    al.GetTag("RG",rgTag);
	    //cout<<rgTag<<endl;
	    if(rg2Tally.find(rgTag) == rg2Tally.end()){ //new
		cerr<<"Unfound new RG "<<rgTag<<endl;
		return 1;
	    }else{
		rg2Tally[rgTag][al.RefID]++;
	    }
	}else{
	    cerr << "Cannot get RG tag for " << al.Name<<endl;
	    //return 1;
	}

	total++;
    } //while al

    reader.Close();
    // writer.Close();
    
    // unmapped.Close();
    cout<<"\t";
    vector<string> toprint;
    for(unsigned i =0;i<references.size();i++){
	toprint.push_back(refData[i].RefName);       
    }
    cout<<vectorToString(toprint,"\t")<<endl;


    map<string, vector<unsigned int> >::iterator rg2TallyIt;
    for (rg2TallyIt =rg2Tally.begin(); 
	 rg2TallyIt!=rg2Tally.end(); 
	 rg2TallyIt++){

	cout<<rg2TallyIt->first<<"\t"<<vectorToString(rg2Tally[rg2TallyIt->first],"\t")<<endl;
	//rg2TallyIt->second->Close();
    }

    cerr<<"Finished succesfully"<<endl;


    return 0;
}

