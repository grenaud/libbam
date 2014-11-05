/*
 * failQualPair
 * Date: Oct-10-2012 
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

     if( (argc== 1) ||
    	(argc== 2 && string(argv[1]) == "-h") ||
    	(argc== 2 && string(argv[1]) == "-help") ||
    	(argc== 2 && string(argv[1]) == "--help") ){
	 cout<<"Usage:sumBases [in bam]"<<endl<<"this program returns the sum of all bases for all reads"<<endl;
	 return 1;
     }

     string bamfiletopen = string(argv[1]);

     BamReader reader;
     // cout<<"ok"<<endl;
     if ( !reader.Open(bamfiletopen) ) {
	 cerr << "Could not open input BAM files." << endl;
	 return 1;
     }

     BamAlignment al;
     unsigned int totalBases=0;
     while ( reader.GetNextAlignment(al) ) {
	 totalBases+=al.QueryBases.size();		   
     } //while al

     reader.Close();

     cout<<totalBases<<endl;

     return 0;
}

