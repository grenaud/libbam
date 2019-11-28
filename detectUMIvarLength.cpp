/*
 * detectUMIvarLength
 * Date: Nov-28-2019
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

//#define DEBUG2


//likelihood variables
double likeMatch[64];
double likeMismatch[64];
double likeMatchProb[64];
double likeMismatchProb[64];
    

double probForQual[64];
double likeRandomMatch;    // 1/4
double likeRandomMisMatch; // 3/4
double likeRandomMatchProb;    // 1/4
double likeRandomMisMatchProb; // 3/4
double max_prob_N ;

void init(){

    max_prob_N =0.25;

    likeRandomMatchProb          =  1.0/4.0 ;
    likeRandomMisMatchProb       =  3.0/4.0 ;

    likeRandomMatch       = log( 1.0/4.0 )/log(10);
    likeRandomMisMatch    = log( 3.0/4.0 )/log(10);

#ifdef DEBUG2
    cout<<"likeRandomMisMatch "<<likeRandomMisMatch<<endl;
    cout<<"likeRandomMatch    "<<likeRandomMatch<<endl;
#endif

    //probability scores for qscores less than 2 make no sense
    //since random DNA has by default a probability of error of 0.75
    //anyway, we set the min for qual scores at 2
    for(int i=0;i<2;i++){

        likeMatch[i]        = log1p(    -pow(10.0,2.0/-10.0) )    /log(10);         
        likeMismatch[i]     = log  (     pow(10.0,2.0/-10.0)/3.0 )/log(10);

        likeMatchProb[i]           = 1.0-pow(10.0,2.0/-10.0) ;
        likeMismatchProb[i]        =     pow(10.0,2.0/-10.0)/3.0 ;
 
        probForQual[i]      = max(double(1.0)-pow(double(10.0),double(2.0)/double(-10.0)),
                                  max_prob_N);
#ifdef DEBUG2
        cout<<endl<<"qual = "<<i<<endl;
        cout<<"match        " <<likeMatch[i]<<endl;
        cout<<"mismatch     " <<likeMismatch[i]<<endl;
        cout<<"match        " <<likeMatchProb[i]<<endl;
        cout<<"mismatch     " <<likeMismatchProb[i]<<endl;
	cout<<"probForQual  "<<probForQual[i]<<endl;
#endif

    }
    //Computing for quality scores 2 and up
    for(int i=2;i<64;i++){
        likeMatch[i]        = log1p(    -pow(10.0,i/-10.0) )     /log(10);          
        likeMismatch[i]     = log  (     pow(10.0,i/-10.0)/3.0  )/log(10); //i/(-10.0*1.0);  


        likeMatchProb[i]           = 1.0-pow(10.0,i/-10.0);
        likeMismatchProb[i]        =     pow(10.0,i/-10.0)/3.0;

        probForQual[i] = max(double(1.0)-pow(double(10.0),double(i)/double(-10.0)),
                             max_prob_N);

#ifdef DEBUG2
        cout<<endl<<"qual = "<<i<<endl;
        cout<<"match    " <<likeMatch[i]<<endl;
        cout<<"mismatch " <<likeMismatch[i]<<endl;
        cout<<"match    " <<likeMatchProb[i]<<endl;
        cout<<"mismatch " <<likeMismatchProb[i]<<endl;
	cout<<"probForQual  "<<probForQual[i]<<endl;
#endif

    }

}

int main (int argc, char *argv[]) {
    init();

    const int qualOffset=33;
    string keyF="TGACT";
    string keyR=reverseComplement(keyF);
    vector<unsigned int> lengthsUMI;
    lengthsUMI.push_back(6);
    lengthsUMI.push_back(12);
    vector<double> ll5p;
    vector<double> ll3p;
    double mininf= - std::numeric_limits<double>::infinity();
    

    unsigned int maxLengthUMI=0;
    for(unsigned int i=0;i<lengthsUMI.size();i++){
	if(maxLengthUMI < lengthsUMI[i]){  maxLengthUMI = lengthsUMI[i];  }
	ll5p.push_back(mininf);
	ll3p.push_back(mininf);
    }

    double randomLL=(maxLengthUMI+keyF.size())*likeRandomMisMatch;

    ll5p.push_back(randomLL);
    ll3p.push_back(randomLL);

     if( (argc== 1) ||
    	(argc== 2 && string(argv[1]) == "-h") ||
    	(argc== 2 && string(argv[1]) == "-help") ||
    	(argc== 2 && string(argv[1]) == "--help") ){
	 cout<<"Usage: detectUMIvarLength [in bam] [outbam]"<<endl<<"This program detect variable length UMI (not customized!)"<<endl;
    	return 1;
     }

     string bamfiletopen = string(argv[1]);
     string bamFileOUT   = string(argv[2]);

     BamReader reader;
     BamWriter writer;

     if ( !reader.Open(bamfiletopen) ) {
	 cerr << "Could not open input BAM files." << endl;
	 return 1;
     }
     const SamHeader header = reader.GetHeader();
     const RefVector references = reader.GetReferenceData();
     if ( !writer.Open(bamFileOUT,header,references) ) {
	 cerr << "Could not open output BAM file "<<bamFileOUT << endl;
	 return 1;
     }

     BamAlignment al;
     BamAlignment al2;
     bool al2Null=true;

     unsigned int trimmedSingles=0;
     unsigned int failedSingles=0;
     unsigned int trimmedPairs=0;
     unsigned int failedPairs=0;

     while ( reader.GetNextAlignment(al) ) {

        if(al.IsPaired() && al2Null ){
            al2=al;
            al2Null=false;
            continue;
        }else{
            if(al.IsPaired() && !al2Null){
                if(al.Name != al2.Name ){
                    cerr << "Seq#1 "<<al.Name<<" has a different id than seq #2 "<<al2.Name<<", exiting " << endl;
                    return 1;
                } 
		//al2=first 
		//al =second

		int ll5i=lengthsUMI.size();
		int ll3i=lengthsUMI.size();

		for(unsigned int i=0;i<lengthsUMI.size();i++){//each UMI length
		    double ll=0;
		    if(al2.QueryBases.length()<lengthsUMI[i]){
			ll      = mininf;//dummy value
			ll5p[i] = mininf;
		    }else{		
			for(unsigned int j=0;j<keyF.size();j++){
			    int idx = lengthsUMI[i] + j;
			    int q   = max(int(char( al2.Qualities[ idx ] ))-qualOffset,2);			    
			    if( al2.QueryBases[ idx ] == keyF[j]){
				ll += likeMatch[    q ];
			    }else{
				ll += likeMismatch[ q ];
			    }
			}
			//add other bases
			ll += (maxLengthUMI)*likeRandomMisMatch;
			ll5p[i]=ll;//dummy value
		    } 
		    //cerr<<i<<" "<<ll<<" "<<ll5p[ll5i]<<endl;
		    if( ll>ll5p[ll5i] ){
			ll5i=i;
		    }
		}

		for(unsigned int i=0;i<lengthsUMI.size();i++){//each UMI length
		    double ll=0;
		    if(al.QueryBases.length()<lengthsUMI[i]){
			ll      = mininf;//dummy value
			ll3p[i] = mininf;
		    }else{			
			for(unsigned int j=0;j<keyR.size();j++){
			    int idx = lengthsUMI[i] + j;
			    int q   = max(int(char( al.Qualities[ idx ] ))-qualOffset,2);			    

			    if( al.QueryBases[ idx ] == keyF[j]){
				ll += likeMatch[    q ];
			    }else{
				ll += likeMismatch[ q ];
			    }
			    
			}
			//add other bases
			ll += (maxLengthUMI)*likeRandomMisMatch;
			ll3p[i]=ll;//dummy value
		    }
		    //cerr<<i<<" "<<ll<<" "<<ll3p[ll3i]<<endl;
		    if( ll>ll3p[ll3i] ){
			ll3i=i;
		    }

		}

		// cerr<<al2.QueryBases<<endl;
		// cerr<<al.QueryBases<<endl;
		// cerr<<"5p "<<vectorToString(ll5p)<<" "<<ll5i<<endl;
		// cerr<<"3p "<<vectorToString(ll3p)<<" "<<ll3i<<endl;

		if(ll5i != int(lengthsUMI.size()) && ll3i != int(lengthsUMI.size()) ){//trim
		    int idxstr=lengthsUMI[ll5i];
		    int idxend=lengthsUMI[ll3i];
		    
		    string umiS         = al2.QueryBases.substr(0,idxstr)+"-"+al.QueryBases.substr(0,idxend);
		    string umiQ         = al2.Qualities.substr( 0,idxstr)+" "+al.Qualities.substr( 0,idxend);
		    //cerr<<"UMI "<<umiS<<endl;
		    if(!al2.AddTag("RX","Z",umiS)){ cerr << "Error while editing RX tags for UMI:"<<umiS<<"#"<< endl; exit(1);   } 
		    if(!al.AddTag("RX","Z",umiS)){  cerr << "Error while editing RX tags for UMI:"<<umiS<<"#"<< endl; exit(1);   } 
		    if(!al2.AddTag("QX","Z",umiQ)){ cerr << "Error while editing QX tags for UMI:"<<umiQ<<"#"<< endl; exit(1);   }
		    if(!al.AddTag("QX","Z",umiQ)){  cerr << "Error while editing QX tags for UMI:"<<umiQ<<"#"<< endl; exit(1);   }
		    
		    int idxstrl =lengthsUMI[ll5i]+keyF.size();
		    int idxendll=lengthsUMI[ll3i]+keyR.size();
		    
		    al2.QueryBases = al2.QueryBases.substr(idxstrl);
		    al2.Qualities  = al2.Qualities.substr( idxstrl);
		    al.QueryBases = al.QueryBases.substr(idxendll);
		    al.Qualities  = al.Qualities.substr( idxendll);

		    // cerr<<al2.QueryBases<<endl;
		    // cerr<<al.QueryBases<<endl;
		    trimmedPairs++;
		}else{
		    failedPairs++;
		    al2.SetIsFailedQC(true);
		    al.SetIsFailedQC(true);
		}

		writer.SaveAlignment(al2);
		writer.SaveAlignment(al);               
                               
            }else{  //  SINGLE END
                //if( al.IsMapped() )
		int ll5i=lengthsUMI.size();
		int ll3i=lengthsUMI.size();

		for(unsigned int i=0;i<lengthsUMI.size();i++){//each UMI length
		    double ll=0;
		    if(al.QueryBases.length()<lengthsUMI[i]){
			ll      = mininf;//dummy value
			ll5p[i] = mininf;
		    }else{		
			for(unsigned int j=0;j<keyF.size();j++){
			    int idx = lengthsUMI[i] + j;
			    int q   = max(int(char( al.Qualities[ idx ] ))-qualOffset,2);			    
			    if( al.QueryBases[ idx ] == keyF[j]){
				ll += likeMatch[    q ];
			    }else{
				ll += likeMismatch[ q ];
			    }
			}
			//add other bases
			ll += (maxLengthUMI)*likeRandomMisMatch;
			ll5p[i]=ll;//dummy value
		    } 
		    //cerr<<i<<" "<<ll<<" "<<ll5p[ll5i]<<endl;
		    if( ll>ll5p[ll5i] ){
			ll5i=i;
		    }
		}

		for(unsigned int i=0;i<lengthsUMI.size();i++){//each UMI length
		    double ll=0;
		    if(al.QueryBases.length()<lengthsUMI[i]){
			ll      = mininf;//dummy value
			ll3p[i] = mininf;
		    }else{			
			for(unsigned int j=0;j<keyR.size();j++){
			    int idx = (al.QueryBases.length()-lengthsUMI[i]-keyR.size()) + j;
			    int q   = max(int(char( al.Qualities[ idx ] ))-qualOffset,2);			    

			    if( al.QueryBases[ idx ] == keyR[j]){
				ll += likeMatch[    q ];
			    }else{
				ll += likeMismatch[ q ];
			    }
			    
			}
			//add other bases
			ll += (maxLengthUMI)*likeRandomMisMatch;
			ll3p[i]=ll;//dummy value
		    }
		    //cerr<<i<<" "<<ll<<" "<<ll3p[ll3i]<<endl;
		    if( ll>ll3p[ll3i] ){
			ll3i=i;
		    }

		}

		//cerr<<al.QueryBases<<endl;
		if(ll5i != int(lengthsUMI.size()) && ll3i != int(lengthsUMI.size()) ){//trim
		    int idxstr=lengthsUMI[ll5i];
		    int idxend=al.QueryBases.size()-lengthsUMI[ll3i];
		    
		    string umiS         = al.QueryBases.substr(0,idxstr)+"-"+al.QueryBases.substr(idxend);
		    string umiQ         = al.Qualities.substr( 0,idxstr)+" "+al.Qualities.substr( idxend);
		    //cerr<<"UMI "<<umiS<<endl;
		    if(!al.AddTag("RX","Z",umiS)){ cerr << "Error while editing RX tags for UMI:"<<umiS<<"#"<< endl; exit(1);   } 
		    if(!al.AddTag("QX","Z",umiQ)){ cerr << "Error while editing QX tags for UMI:"<<umiQ<<"#"<< endl; exit(1);   }
		    
		    int idxstrl =lengthsUMI[ll5i]+keyF.size();
		    int idxendll=(lengthsUMI[ll3i]+keyR.size());
		    
		    al.QueryBases = al.QueryBases.substr(idxstrl,al.QueryBases.size()-idxstrl-idxendll);
		    al.Qualities  = al.Qualities.substr( idxstrl,al.Qualities.size() -idxstrl-idxendll);
		    //cerr<<al.QueryBases<<endl;
		    trimmedSingles++;
		}else{
		    failedSingles++;
		    al.SetIsFailedQC(true);
		}

		// cerr<<"5p "<<vectorToString(ll5p)<<" "<<ll5i<<endl;
		// cerr<<"3p "<<vectorToString(ll3p)<<" "<<ll3i<<endl;

		writer.SaveAlignment(al);
            } //end single end

            al2Null=true;
        }//second pair

	 // al.SetIsFailedQC(true);
	 // writer.SaveAlignment(al);
		 
     } //while al

     reader.Close();
     writer.Close();


     cerr<<"Done: trimmedSingles "<<trimmedSingles<<" "<<
	 "failedSingles "<<failedSingles<<" "<<   
	 "trimmedPairs "<<trimmedPairs<<" "<<
	 "failedPairs "<< failedPairs<<endl;

     return 0;
}

