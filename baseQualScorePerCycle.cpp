#include <iostream>
#include <vector>
#include <set>
#include <ctype.h>
#include <stdlib.h>

#include "api/BamMultiReader.h"
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"
#include "api/SamSequenceDictionary.h"

#include "libgab.h"
#include "ReconsReferenceBAM.h"
#define maxqual 64
#define maxbase 16

using namespace std;
using namespace BamTools;

const int offset=33;


int numberOfCycles;


// bool isTransition(char reference,char read){

//     if(  (reference == 'A' && read == 'G') ||
// 	 (reference == 'G' && read == 'A') ||
// 	 (reference == 'C' && read == 'T') ||
// 	 (reference == 'T' && read == 'C') ){
// 	return true;
//     }
	 
//     return false;
// }



// bool isTransversions(char reference,char read){
//     if(  (reference == 'A' && read == 'C') ||
// 	 (reference == 'A' && read == 'T') ||
	 
// 	 (reference == 'G' && read == 'C') ||
// 	 (reference == 'G' && read == 'T') ||
	 
// 	 (reference == 'C' && read == 'A') ||
// 	 (reference == 'C' && read == 'G') ||

// 	 (reference == 'T' && read == 'A') ||
// 	 (reference == 'T' && read == 'G') ){
// 	return true;
//     }

//     return false;	
// }

// vector<unsigned int>    counterOfOccurences;

// vector<unsigned int> mismatches;

// typedef unsigned int [maxqual] qualScoresCount;
// typedef struct qualScoresCount { unsigned int [ maxqual ]  } qualScoresCount;
// typedef unsigned int  qualScoresCount [maxqual];
// typedef unsigned int  pairBP [ maxbase ];

vector<          vector<unsigned int>   >     counterOfOccurencesPerCycle; // first dim = cycle, second is allelePair2Int
vector< vector<  vector<unsigned int> > >     counterOfQScoresPerCycle;    // first dim = cycle, second is allelePair2Int, third is the quality score


//increases the counters mismatches and typesOfMismatches of a given BamAlignment object
inline void increaseCounters(BamAlignment & al,string & reconstructedReference,int firstCycleRead,int increment){

    char refeBase;
    char readBase;
    unsigned int qualScore;
    int cycleToUse=firstCycleRead;
    // cout<<"name "<<al.Name<<endl;
    // cout<<"firstCycleRead "<<firstCycleRead<<endl;
    // cout<<"increment      "<<increment<<endl;

    for(int i=0;i<numberOfCycles;i++,cycleToUse+=increment){
        // cout<<"i = "<<i<<" cyc "<<cycleToUse<<endl;
	refeBase  = toupper(reconstructedReference[i]);
	readBase  = toupper(         al.QueryBases[i]);
	qualScore = int(al.Qualities[i]-offset);

	if(refeBase == 'S' ||refeBase == 'I'){ //don't care about soft clipped or indels
	    continue;
	}


	// //match
	 if(refeBase == 'M'){
	     //matches[cycleToUse]++;

	     if(al.IsReverseStrand()){ //need to take the complement
		readBase=complement(readBase);
	    }

	     //cout<<readBase<<endl;
	     counterOfOccurencesPerCycle[cycleToUse][allelePair2Int(readBase,readBase)]++;
	     counterOfQScoresPerCycle[cycleToUse][allelePair2Int(readBase,readBase)][qualScore]++;
	     continue;
	}
	    
	//mismatch
	if( isResolvedDNA(refeBase)  && 
	    isResolvedDNA(readBase) ){
	    
	    if(al.IsReverseStrand()){ //need to take the complement
		refeBase=complement(refeBase);
		readBase=complement(readBase);
	    }
	 
	    if(readBase == refeBase){
		cerr<<"Internal error in reconstruction of read "<<al.Name<<", contact developer"<<endl;
		exit(1);;
	    }
					
	    counterOfOccurencesPerCycle[cycleToUse][allelePair2Int(refeBase,readBase)]++;
	    counterOfQScoresPerCycle[cycleToUse][allelePair2Int(refeBase,readBase)][qualScore]++;

	    continue;
	}	    		     
    }
}


int main (int argc, char *argv[]) {

    string usage=string(""+string(argv[0])+"  [in BAM file]"+
			"\nThis program reads a BAM file and computes the error rate for each cycle\n"+
			// "\nreads and the puts the rest into another bam file.\n"+
			// "\nTip: if you do not need one of them, use /dev/null as your output\n"+
			// "arguments:\n"+
			// "\t"+"--bq  [base qual]   : Minimum base quality to flag a deaminated site (Default: "+stringify(minBaseQuality)+")\n"+
			"\n");

    if(argc == 1 ||
       (argc == 2 && (string(argv[0]) == "-h" || string(argv[0]) == "--help") )
       ){
	cerr << "Usage "<<usage<<endl;
	return 1;       
    }


    // for(int i=1;i<(argc-2);i++){ 

	
    //     if(string(argv[i]) == "--bq"){
    // 	    minBaseQuality=destringify<int>(argv[i+1]);
    //         i++;
    //         continue;
    // 	}

    // }

    string bamfiletopen = string( argv[ argc-1 ] );
    // string deambam      = string( argv[ argc-2 ] );
    // string nondeambam   = string( argv[ argc-1 ] );


    BamReader reader;
    
    if ( !reader.Open(bamfiletopen) ) {
    	cerr << "Could not open input BAM file"<< bamfiletopen << endl;
    	return 1;
    }








    //iterating over the alignments for these regions
    BamAlignment al;
    bool pairedEnd=false;
    bool firstRead=true;

    while ( reader.GetNextAlignment(al) ) {

	if(firstRead){ //reads are either all paired end or single end, I don't allow a mix
	    numberOfCycles=al.QueryBases.size();
	    // cout<<"numberOfCycles "<<numberOfCycles<<endl;
	    if(al.IsPaired() ){  
		pairedEnd=true;		

		counterOfOccurencesPerCycle = vector<  vector<unsigned int>  >             (2*numberOfCycles, vector<unsigned int> (maxbase,0) );
		counterOfQScoresPerCycle    = vector<  vector< vector<unsigned int> >  >   (2*numberOfCycles, vector< vector<unsigned int> > (maxbase,vector<unsigned int>(maxqual,0)));
	    }else{

		counterOfOccurencesPerCycle = vector<  vector<unsigned int>  >             (numberOfCycles, vector<unsigned int> (maxbase,0) );
		counterOfQScoresPerCycle    = vector<  vector< vector<unsigned int> >  >   (numberOfCycles, vector< vector<unsigned int> > (maxbase,vector<unsigned int>(maxqual,0)));
		
	    }
	    firstRead=false;
	}

	if( (  pairedEnd && !al.IsPaired()) ||  
	    ( !pairedEnd &&  al.IsPaired())   ){
	    cerr<<"Read "<<al.Name<<" is wrong, cannot have a mixture of paired and unpaired read for this program"<<endl;
	    return 1;
	}


	//skip unmapped
	if(!al.IsMapped())
	    continue;

	if(numberOfCycles!=int(al.QueryBases.size())){
	    cerr<<"The length of read "<<al.Name<<" is wrong, should be "<<numberOfCycles<<"bp"<<endl;
	    return 1;
	}


	string reconstructedReference = reconstructRef(&al);
	if(al.Qualities.size() != reconstructedReference.size()){
	    cerr<<"Quality line is not the same size as the reconstructed reference"<<endl;
	    return 1;
	}


	if( pairedEnd ){
	    if( al.IsFirstMate() ){ //start cycle 0

		if( al.IsReverseStrand()  ){ 
		    increaseCounters(al,reconstructedReference,numberOfCycles-1,-1); //start cycle numberOfCycles-1
		}else{
		    increaseCounters(al,reconstructedReference,0               , 1); //start cycle 0
		}
       
	    }else{

		if( al.IsSecondMate()  ){ 
		    if( al.IsReverseStrand()  ){ 
			increaseCounters(al,reconstructedReference,2*numberOfCycles-1,-1); //start cycle 2*numberOfCycles-1
		    }else{
			increaseCounters(al,reconstructedReference,numberOfCycles    , 1); //start cycle numberOfCycles
		    }
		}else{
	    	    cerr<<"Reads "<<al.Name<<" must be either first or second mate"<<endl;
	    	    return 1;
	    	}
	    }	
	}else{ //single end 

	    if( al.IsReverseStrand()  ){ 
		increaseCounters(al,reconstructedReference,numberOfCycles-1,-1); //start cycle numberOfCycles-1
	    }else{
		increaseCounters(al,reconstructedReference,0               , 1); //start cycle 0
	    }
	    
	}


    }//end while  each read
	


    reader.Close();

//     cout<<"cycle\tmatches\tmismatches\tmismatches%\tA>C\tA>C%\tA>G\tA>G%\tA>T\tA>T%\tC>A\tC>A%\tC>G\tC>G%\tC>T\tC>T%\tG>A\tG>A%\tG>C\tG>C%\tG>T\tG>T%\tT>A\tT>A%\tT>C\tT>C%\tT>G\tT>G%"<<endl;


    cout<<"cycle\tA>A\tA>C\tA>G\tA>T\tC>A\tC>C\tC>G\tC>T\tG>A\tG>C\tG>G\tG>T\tT>A\tT>C\tT>G\tT>T\tA>A(qc)\tA>C(qc)\tA>G(qc)\tA>T(qc)\tC>A(qc)\tC>C(qc)\tC>G(qc)\tC>T(qc)\tG>A(qc)\tG>C(qc)\tG>G(qc)\tG>T(qc)\tT>A(qc)\tT>C(qc)\tT>G(qc)\tT>T(qc)"<<endl;


    
    for(unsigned int i=0;i<counterOfQScoresPerCycle.size();i++){
	cout<<(i+1)<<"\t"<<vectorToString(counterOfOccurencesPerCycle[i],"\t");
	for(unsigned int j=0;j<maxbase;j++){	    	   
	    cout<<"\t"<<vectorToString(counterOfQScoresPerCycle[i][j],",");
	}
	cout<<endl;
    }

    // for(unsigned int j=0;j<maxbase;j++){	    
    // }    }
    
    //     for(unsigned int i=0;i<matches.size();i++){
    // cout<<(i+1);
    // 	if( (matches[i]+mismatches[i]!=0) )
    // 	    cout<<"\t"<<matches[i]<<"\t"<<mismatches[i]<<"\t"<< 100.0*(double(mismatches[i])/double(matches[i]+mismatches[i])) ;
    // 	else
    // 	    cout<<"\t"<<matches[i]<<"\t"<<mismatches[i]<<"\tNA";

    // 	for(int j=0;j<12;j++){   
    // 	    cout<<"\t"<<typesOfMismatches[j][i];
    // 	    if( (matches[i]+mismatches[i]!=0) )
    // 		cout<<"\t"<<100.0*double(typesOfMismatches[j][i])/double(matches[i]+mismatches[i]);
    // 	    else
    // 		cout<<"\tNA";
    // 	}

    // 	cout<<endl;
    //     }




   
    return 0;
}

