/*
 * algo_init.cpp
 *
 *  Created on: Aug 13, 2014
 *      Author: Qi
 */

#include "algo_init.h"
#include "globals.h"
#include "profiling.h"
#include "utils.h"
#include "hash.h"
#include <cstring>
#include <string>
#include <fstream>
#include <algorithm>
#include <string>
#include <bitset>

using namespace std;

typedef pair<string, int> PAIR;  

struct pinfo{
	unsigned int did;
	// float s;
	float s1;
	float s2;
};

struct sinfo{
	unsigned int did;
	// unsigned int freq;
	float score;
};

struct scores{
	float s1;
	float s2;
	int f1;
	int f2;
};

struct fullinfo{
	int did;
	float score;
	short kbits;
};

bool myfunc (const sinfo& a, const sinfo& b){
	return (a.score > b.score);
}

bool sortcdt (const fullinfo& a, const fullinfo& b){
	return (a.score > b.score);
}	

bool cmp_by_value(const PAIR& lhs, const PAIR& rhs) {  
  return lhs.second < rhs.second;  
}  

void algo_init::operator() (CluewebReader* Reader, int qn, pairlists& pls, lptrArray& lps, const int topK, QpResult* res, profilerC& p) {


	const int p_l = pls.lengths.size();
	const int t_l = lps.size();
	const int a_l = pls.lengths.size() + lps.size();
	cout<<"origin singleterms: "<<lps.size()<<endl;
	cout<<"origin pairs: "<<pls.lengths.size()<<endl;

	string queryline;

	short filter = (1<<t_l) - 1;
	bitset<8> y(filter);
	cout<<y<<endl;

	/*total # of docids, used for hashtable size*/
	size_t totaldids = 0;

	/*used to parse term1+term2 in pair*/
	string dem = "+";
	string term1, term2;
	size_t position = 0;


	/*each term has a # with a binary bits representation, 0 indicates existing, sorted by impact, the lower the bit is, the higher the impact*/
	map<string, int> termmapping;

	/*each term with the sum of its top1 scores, will be used to assign the binary bits*/
	vector<PAIR> tllv;

	/*pairmap has the pinfo vector stores all the dids in slices seq, while in each slice ordered by inmpact*/
	/*psize map has #=slice elem, each elem is the # of dids within each slice*/
	map<string, vector<pinfo>> pairmap;
	map<string, vector<size_t>> psizemap;

	/*singlemap has the pinfo vector stores all the dids in slices seq, while in each slice ordered by inmpact*/
	/*ssize map has #=slice elem, each elem is the # of dids within each slice*/
	map<string, vector<sinfo>> singlemap;
	map<string, vector<size_t>> ssizemap;

	/*# of slices, build the boundries for each slice*/
	const int slices = 6;
	int bs[slices+1];
	bs[0] = 0;
	bs[slices] = CONSTS::MAXD;
	int jump = CONSTS::MAXD/slices;

	for(int i=1; i<slices; i++){
		bs[i] = bs[i-1] + jump;
	}

	for(int i=0; i<=slices; i++){
		cout<<"slice # "<<i<<": "<<bs[i]<<endl;
	}

	vector<pinfo> tempp_vu;
	vector<pinfo> tempp_vo;
	vector<pinfo> tempp_v;
	vector<size_t> psize;

	vector<sinfo> temps_v;
	vector<size_t> ssize;


	/*load the singlelist into map*/
	for (int i=0; i<lps.size(); ++i){



		const string term = lps[i]->term;
		const int length = lps[i]->unpadded_list_length;

		// termmapping[term] = 0; //for termmapping, just insert the term in the map, the value 0 doesn't matter
		tllv.push_back(make_pair(term, length)); //for termmapping, pair value is the listlength


		const int l_sqrt = sqrt(lps[i]->unpadded_list_length);
		// cout<<term<<": "<<length<<endl;

		if( lps.size()>2 && length >= 10000000){
			cout<<"Query size greater than 3, throw away singles larger than 10m: "<< term <<endl;
			continue;
		}

		int term_id = Reader->term_map[term];
		// cout<<term<<" "<<term_id<<endl;
	 	RawIndexList Rlist = Reader->load_raw_list(term,term_id);

	 	vector<sinfo> t_list;

	 	for(int h = 0; h < length; h++){
	 		sinfo t_p;
	 		t_p.did  = Rlist.doc_ids.at(h);
	 		// t_p.freq = Rlist.freq_s.at(h);
	 		t_p.score = Rlist.scores.at(h);
	 		t_list.push_back(t_p);
	 	}

	 	sort(t_list.begin(), t_list.end(), myfunc);

	 	size_t cut = 1000;
	 	if(length > cut)
	 	t_list.resize(cut);// resize according to Configuration file
	 	totaldids += cut;

	 	// cout<<t_list.size()<<endl;

		ssize.push_back(0); //the first should be 0;

		for(int h=0; h<slices; h++){
			int ct = 0;
	 	  for(int j = 0; j < t_list.size(); j++){


	 		if(t_list.at(j).did > bs[h] && t_list.at(j).did <= bs[h + 1]){
	 			temps_v.push_back(t_list.at(j));
	 			ct ++;
	 		}

	 	  }
	 	  	// cout<<"ct: "<<ct<<endl;
	 	    ssize.push_back(ct+ssize.back());
		}



	 	singlemap[term] = temps_v;
	 	ssizemap[term] = ssize;
	 	cout<<term<<": "<<temps_v.size()<<endl;
	 	cout<<term<<": "<<ssize.size() - 1<<endl;

  		temps_v.clear();
  		ssize.clear();

	}

		unsigned short offset = 0;
		unsigned short elem = 0;
		hashTable *ht;
		if( 5*totaldids/slices < 8191)
		ht = initHash(totaldids/slices, 0); 
		else
		ht = initHash(8191, 1); //table size 8191
		// vector<int> didresults;
		// vector<float> scoreresults;
		// vector<short> kbitsresults;
		// didresults.reserve(15000);
		// scoreresults.reserve(15000);
		// kbitsresults.reserve(15000);
		vector<fullinfo> fresults;
		fresults.reserve(15000);


		/*term mapping according to listlen*/
		int mapping = 0;

		sort(tllv.begin(), tllv.end(), cmp_by_value);

		for(int i=0; i<tllv.size(); i++){
			cout<<tllv.at(i).first<<": "<<tllv.at(i).second<<endl;
			termmapping[tllv.at(i).first] = ~ (1<<mapping++);
		}

		// for(map<string, int>::iterator it = termmapping.begin(); it!=termmapping.end(); ++it){
		// 	it->second = ~ (1<<mapping++) ;
		// }

		for(map<string, int>::iterator it = termmapping.begin(); it!=termmapping.end(); ++it){

			bitset<8> x(it->second);
			cout<<it->first<<": "<<x<<endl;

		}
		/*term mapping according to listlen*/


	/*---load the pair lists into map----*/
  	string did_s;
  	int did;
  	string sc1_s;
  	float sc1;
  	string sc2_s;
  	float sc2;
  	int count = 0;
  	// string ts_s;
  	// float ts;

  	if(pls.lengths.size()>0){/*if pair list size is not 0*/

  	for(int k = 0; k < pls.lengths.size(); k++){

  	position = pls.pairnames.at(k).find(dem);
	term1 = pls.pairnames.at(k).substr(0, position);
	term2 = pls.pairnames.at(k).substr(position+1,pls.pairnames.at(k).size());
	// bitset<8> x(termmapping[term1]);
	// bitset<8> y(termmapping[term2]);
	// bitset<8> z(((~termmapping[term1])>>3) | ((~termmapping[term2])>>3));
	// cout<<"in loading pairs stage: "<<term1<<": "<<x<<", "<<term2<<": "<<y<<": "<<z<<endl;

	if( ( ((~termmapping[term1])>>3) | ((~termmapping[term2])>>3) ) > 0)
		continue;

  	ifstream pair_stream;
  	string pair_dir = CONSTS::pair_index + pls.pairnames.at(k);
  	// string pair_dir = "/data/qw376/pair_index/" + pls.pairnames.at(k);
  	pair_stream.open(pair_dir.c_str());

  	count = 0;

  	// cout<<pls.pairnames.at(k)<<endl;

  	while(getline(pair_stream, queryline)){

  		pinfo temp;

  		// cout<<queryline<<endl;
  		string::iterator itr = queryline.begin();
		string::iterator start = itr;

   		while(itr != queryline.end() && !isspace(*itr)){
			++itr;
		}

		did_s = string(start, itr);
		did = atoi(did_s.c_str());

		//take the total score
		start = itr+1;
  	  	itr++;
    	while(itr != queryline.end() && !isspace(*itr)){
			++itr;
		}
		// ts_s = string(start, itr);
		// ts = atof(ts_s.c_str());

		//ignore the first freq
		start = itr+1;
  	  	itr++;
    	while(itr != queryline.end() && !isspace(*itr)){
			++itr;
		}

		//take the first score
		start = itr+1;
  	  	itr++;
    	while(itr != queryline.end() && !isspace(*itr)){
			++itr;
		}

		sc1_s = string(start, itr);
		sc1 = atof(sc1_s.c_str());

		//ignore the second freq
		start = itr+1;
  	  	itr++;
    	while(itr != queryline.end() && !isspace(*itr)){
			++itr;
		}

		//take the second score
		start = itr+1;
  	  	itr++;
    	while(itr != queryline.end() && !isspace(*itr)){
			++itr;
		}

		sc2_s = string(start, itr);
		sc2 = atof(sc2_s.c_str());


		temp.did = did;
		// temp.s = ts;
		temp.s1 = sc1;
		temp.s2 = sc2;

		tempp_vo.push_back(temp);

  		count ++;		
		// t_map[didl] = count;
  		if (count == pls.lengths.at(k)){
  			totaldids += count;
  			break;
  		}

  	}


  		psize.push_back(0); //the first should be 0;

  		// cout<<pls.pairnames.at(k)<<": "<<tempp_vo.size()<<endl;

  		for(int j=0; j<slices; j++){
  			int ct = 0;
	 		for(int i=0; i<tempp_vo.size(); i++){	
	 		  if(tempp_vo.at(i).did > bs[j] && tempp_vo.at(i).did <= bs[j + 1]){
	 			tempp_v.push_back(tempp_vo.at(i));
	 			ct ++;
	 		  }
	 	  }
	 	  	// cout<<"ct: "<<ct<<endl;
	 	    psize.push_back(ct+psize.back());
	 	}


	 	pairmap[pls.pairnames.at(k)] = tempp_v;
	 	psizemap[pls.pairnames.at(k)] = psize;
	 	cout<<pls.pairnames.at(k)<<": "<<tempp_v.size()<<endl;
	 	cout<<pls.pairnames.at(k)<<": "<<psize.size() - 1<<endl;

	 	tempp_vo.clear();
  		tempp_v.clear();
  		psize.clear();

  	  	pair_stream.close();
	}

	}//if pair list size is not 0
	/*-----------------------*/

		cout<<"singleterms after purge: "<<singlemap.size()<<endl;
	    cout<<"pairs after purge: "<<pairmap.size()<<endl;

		int hit = 0;

		p.start(CONSTS::ALLQS);


		for(int k=0; k<slices; k++){

		/*for the pair lists hashing*/
		/*the reason put pairs before singles is that pairs is loaded later than singles in loading process, kind of cache warm-up*/
	
		for(map<string, vector<pinfo>>::iterator it = pairmap.begin(); it!=pairmap.end(); ++it){
			// cout<<k+1<<" pair try: "<<it->first<<" "<<psizemap[it->first].at(k+1)-psizemap[it->first].at(k)<<endl;
			// cout<<psizemap[it->first].at(k+1)<<" - "<<psizemap[it->first].at(k)<<endl;
			position = it->first.find(dem);
			term1 = it->first.substr(0, position);
			term2 = it->first.substr(position+1,it->first.size());
			// cout<<term1<<" "<<term2<<endl; 
			// for(int i=0; i<it->second.size(); i++){
			for(int i=psizemap[it->first].at(k); i<psizemap[it->first].at(k+1); i++){

				// cout<<it->second.at(i).did<<": "<<it->second.at(i).s1<<" "<<it->second.at(i).s2<<endl;
				// int pos = insertHash(ht, it->second.at(i).did, elem, 0, didresults);
				// int pos = insertHash(ht, it->second.at(i).did, elem, 0, didresults, hit);
				int pos = insertHash(ht, it->second.at(i).did, elem, 0, fresults);
    			if (pos){ //key already existed
     				   // printf("Key %d already exists!\n", GETKEY(ht->table[pos-1]-1, didresults));

    				   // scoreresults.at(ht->table[pos-1]-1) = scoreresults.at(ht->table[pos-1]-1) + it->second.at(i).s1 * ( kbitsresults.at(ht->table[pos-1]-1) & (termmapping[term1]) != filter) + it->second.at(i).s2 * ( kbitsresults.at(ht->table[pos-1]-1) & (termmapping[term2]) != filter);
    				   // kbitsresults.at(ht->table[pos-1]-1) = (kbitsresults.at(ht->table[pos-1]-1)) & (termmapping[term1]) & (termmapping[term2]);

    					fresults.at(ht->table[pos-1]-1).score = fresults.at(ht->table[pos-1]-1).score + it->second.at(i).s1 * ( fresults.at(ht->table[pos-1]-1).kbits & (termmapping[term1]) != filter) + it->second.at(i).s2 * ( fresults.at(ht->table[pos-1]-1).kbits & (termmapping[term2]) != filter);
    					fresults.at(ht->table[pos-1]-1).kbits = (fresults.at(ht->table[pos-1]-1).kbits) & (termmapping[term1]) & (termmapping[term2]);

   					 }else{ //successfully inserted

   					 	// didresults.push_back(it->second.at(i).did);
   					 	// scoreresults.push_back(it->second.at(i).s1 + it->second.at(i).s2);
   					 	// kbitsresults.push_back( (termmapping[term1]) & (termmapping[term2]) );


   					 	fullinfo ftemp;
   					 	ftemp.did = it->second.at(i).did;
   					 	ftemp.score = it->second.at(i).s1 + it->second.at(i).s2;
   					 	ftemp.kbits = (termmapping[term1]) & (termmapping[term2]);
   					 	fresults.push_back(ftemp);


     			 		elem = ++ offset;
   					 	// printf("successfully inserted!\n");
   				 }
			}
		}/*pairs finished*/

		/*for the single lists hashing*/
		for(map<string, vector<sinfo>>::iterator it = singlemap.begin(); it!=singlemap.end(); ++it){

			// cout<<k+1<<" single try: "<<it->first<<" "<<ssizemap[it->first].at(k+1)-ssizemap[it->first].at(k)<<endl;
			// cout<<ssizemap[it->first].at(k+1)<<" - "<<ssizemap[it->first].at(k)<<endl;
			// for(int i=0; i<it->second.size(); i++){
			   for(int i=ssizemap[it->first].at(k); i<ssizemap[it->first].at(k+1); i++){

				// cout<<it->second.at(i).did<<": "<<it->second.at(i).score<<endl;
				// int pos = insertHash(ht, it->second.at(i).did, elem, 0, didresults);
				// int pos = insertHash(ht, it->second.at(i).did, elem, 0, didresults, hit);
				int pos = insertHash(ht, it->second.at(i).did, elem, 0, fresults);
    			if (pos){ //key already existed
     				    // printf("Key %d already exists!\n", GETKEY(ht->table[pos-1]-1, didresults));
    					// cout<<it->second.at(i).did<<" key already existed at: "<<ht->table[pos-1]<<" pos: "<<pos<<endl;

    					// scoreresults.at(ht->table[pos-1]-1) = scoreresults.at(ht->table[pos-1]-1) + it->second.at(i).score;
    					// kbitsresults.at(ht->table[pos-1]-1) = (kbitsresults.at(ht->table[pos-1]-1)) & (termmapping[it->first]);

    					fresults.at(ht->table[pos-1]-1).score = fresults.at(ht->table[pos-1]-1).score + it->second.at(i).score;
    					fresults.at(ht->table[pos-1]-1).kbits = fresults.at(ht->table[pos-1]-1).kbits & (termmapping[it->first]);

   					 }else{ //successfully inserted

   					 	// didresults.push_back(it->second.at(i).did);
   					 	// scoreresults.push_back(it->second.at(i).score);
   					 	// kbitsresults.push_back(termmapping[it->first]);

   					 	fullinfo ftemp;
   					 	ftemp.did = it->second.at(i).did;
   					 	ftemp.score = it->second.at(i).score;
   					 	ftemp.kbits = termmapping[it->first];
   					 	fresults.push_back(ftemp);

     			 		elem = ++ offset;
   					 	// printf("successfully inserted!\n");
   				 }
			}
		}/*singles finished*/

			clearHash(ht);
		}
		/*pairlist*/

		sort(fresults.begin(), fresults.end(), sortcdt);
		fresults.resize(200);

		p.end(CONSTS::ALLQS);



		// cout<<"final did #: "<<didresults.size()<<endl;
		cout<<"final did #: "<<fresults.size()<<endl;
		cout<<"final did #: "<<offset<<endl;
		// cout<<"attempts: "<<hit<<endl;
		// cout<<"rate: "<<(float)hit/(float)offset<<endl;
		// cout<<"final did #: "<<fresults.size()<<endl;

		// for(int i=0; i<didresults.size(); i++){
		// 	bitset<8> x(kbitsresults.at(i));
		// 	cout<<didresults.at(i)<<", "<<scoreresults.at(i)<<", "<<x<<endl;
		// }

		int lookups = 0;

		ofstream out_stream;
		out_stream.open(CONSTS::Candidates_200.c_str(), ofstream::app);

		
		for(int i=0; i<fresults.size(); i++){

			int n = ~(fresults.at(i).kbits);
			count = 0;
			while (n>0) { 
				count = count + (n&1);
				n=n>>1; //Right shift by 1 
			}

    		lookups += t_l - count; 

    		// bitset<8> x(fresults.at(i).kbits);
			// cout<<fresults.at(i).did<<", "<<fresults.at(i).score<<", "<<x<<", "<<t_l - count<<endl;
			out_stream<<fresults.at(i).did<<" ";
		}
		out_stream<<endl;
		out_stream.close();

		cout<<"Total lookups needed: "<<lookups<<endl;

}
