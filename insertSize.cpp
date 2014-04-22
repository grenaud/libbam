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
	 cout<<"Usage:editDist [in bam]"<<endl<<"this program returns the NM field of all aligned reads"<<endl;
	 return 1;
     }

     string bamfiletopen = string(argv[1]);
     // cout<<bamfiletopen<<endl;
     BamReader reader;
     // cout<<"ok"<<endl;
     if ( !reader.Open(bamfiletopen) ) {
	 cerr << "Could not open input BAM files." << endl;
	 return 1;
     }

     BamAlignment al;
     // cout<<"ok"<<endl;
     while ( reader.GetNextAlignment(al) ) {
	 // cout<<al.Name<<endl;
	 if(!al.IsMapped())
	     continue;

	 if(al.IsPaired()){
	     if( al.IsFirstMate()  ){
		 cout<<al.InsertSize<<endl;
	     }
	 }else{
	     cout<<al.Length<<endl;
	 }

		    

     } //while al

     reader.Close();

     return 0;
}

