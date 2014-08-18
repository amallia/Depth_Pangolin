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

using namespace std;

struct pinfo{
	unsigned int did;
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


bool myfunc (sinfo a, sinfo b){
	return (a.score > b.score);
}

void algo_init::operator() (CluewebReader* Reader, int qn, pairlists& pls, lptrArray& lps, const int topK, QpResult* res, profilerC& p) {


	const int p_l = pls.lengths.size();
	const int t_l = lps.size();
	const int a_l = pls.lengths.size() + lps.size();
	cout<<"singleterms: "<<lps.size()<<endl;
	cout<<"pairs: "<<pls.lengths.size()<<endl;

	string queryline;

	// map<string, vector<pinfo>> pairmap1;
	// map<string, vector<pinfo>> pairmap2;
	// map<string, vector<pinfo>> pairmap3;
	// map<string, vector<pinfo>> pairmap4;
	// map<string, vector<pinfo>> pairmap5;
	map<string, vector<pinfo>> pairmap;

	vector<pinfo> tempp_v1;
	vector<pinfo> tempp_v2;
	vector<pinfo> tempp_v3;
	vector<pinfo> tempp_v4;
	vector<pinfo> tempp_v5;
	vector<pinfo> tempp_v;

	// map<string, vector<sinfo>> singlemap1;
	// map<string, vector<sinfo>> singlemap2;
	// map<string, vector<sinfo>> singlemap3;
	// map<string, vector<sinfo>> singlemap4;
	// map<string, vector<sinfo>> singlemap5;
	map<string, vector<sinfo>> singlemap;

	vector<sinfo> temps_v1;
	vector<sinfo> temps_v2;
	vector<sinfo> temps_v3;
	vector<sinfo> temps_v4;
	vector<sinfo> temps_v5;
	vector<sinfo> temps_v;

	// const int tablesize = 3079;
	// int didtable[tablesize];
	// float scoretable[tablesize];
	// bool kbittable[tablesize][t_l];

	// for(int i = 0; i< tablesize; i++){
	// 	didtable[i] = -1;
	// }


    char * docid_s;
	// int docid;
	int count = 0;
	int test;


  	/*---load the pair lists into map----*/
  	string did_s;
  	int did;
  	string sc1_s;
  	float sc1;
  	string sc2_s;
  	float sc2;

  	if(pls.lengths.size()>0){

  	for(int k = 0; k < pls.lengths.size(); k++){


  	ifstream pair_stream;
  	string pair_dir = CONSTS::pair_index + pls.pairnames.at(k);
  	// string pair_dir = "/data/qw376/pair_index/" + pls.pairnames.at(k);
  	pair_stream.open(pair_dir.c_str());

  	count = 0;

  	cout<<pls.pairnames.at(k)<<endl;

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

		//ignore the total score
		start = itr+1;
  	  	itr++;
    	while(itr != queryline.end() && !isspace(*itr)){
			++itr;
		}

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
		temp.s1 = sc1;
		temp.s2 = sc2;


		if(did <= 10000000)
			tempp_v1.push_back(temp);
		if(did > 10000000 && did <= 20000000)
			tempp_v2.push_back(temp);
		if(did > 20000000 && did <= 30000000)
			tempp_v3.push_back(temp);
		if(did > 30000000 && did <= 40000000)
			tempp_v4.push_back(temp);
		if(did > 40000000 && did <= CONSTS::MAXD)
			tempp_v5.push_back(temp);

  		count ++;		
		// t_map[didl] = count;
  		if (count == pls.lengths.at(k)){
  			// pairmap1[pls.pairnames.at(k)] = tempp_v1;
  			// pairmap2[pls.pairnames.at(k)] = tempp_v2;
  			// pairmap3[pls.pairnames.at(k)] = tempp_v3;
  			// pairmap4[pls.pairnames.at(k)] = tempp_v4;
  			// pairmap5[pls.pairnames.at(k)] = tempp_v5;
  			break;
  		}

  	}

  		tempp_v.insert(tempp_v.end(),tempp_v1.begin(),tempp_v1.end());
	 	tempp_v.insert(tempp_v.end(),tempp_v2.begin(),tempp_v2.end());
	 	tempp_v.insert(tempp_v.end(),tempp_v3.begin(),tempp_v3.end());
	 	tempp_v.insert(tempp_v.end(),tempp_v4.begin(),tempp_v4.end());
	 	tempp_v.insert(tempp_v.end(),tempp_v5.begin(),tempp_v5.end());


	 	// for(int i = 0; i<tempp_v.size(); i++){
	 	// 	cout<<tempp_v.at(i).did<<endl;
	 	// }

	 	pairmap[pls.pairnames.at(k)] = tempp_v;

	 	tempp_v1.clear();
  		tempp_v2.clear();
  		tempp_v3.clear();
  		tempp_v4.clear();
  		tempp_v5.clear();
  		tempp_v.clear();

  	  	pair_stream.close();
	}

	}//if pair list size is not 0
	/*-----------------------*/

	// cout<<"map1 size: "<<pairmap1.size()<<endl;
	// cout<<"map2 size: "<<pairmap2.size()<<endl;
	// cout<<"map3 size: "<<pairmap3.size()<<endl;
	// cout<<"map4 size: "<<pairmap4.size()<<endl;
	// cout<<"map5 size: "<<pairmap5.size()<<endl;

	// for(map<string, vector<pinfo>>::iterator it = pairmap1.begin(); it!=pairmap1.end(); ++it){
	// 	cout<< it->first << ": " <<it->second.size()<<endl;
	// 	for(int i=0; i<it->second.size(); i++){
	// 		cout<<it->second.at(i).did<<": "<<it->second.at(i).s1<<" "<<it->second.at(i).s2<<endl;
	// 	}
	// }

	// for(map<string, vector<pinfo>>::iterator it = pairmap2.begin(); it!=pairmap2.end(); ++it){
	// 	cout<< it->first << ": " <<it->second.size()<<endl;
	// 	for(int i=0; i<it->second.size(); i++){
	// 		cout<<it->second.at(i).did<<": "<<it->second.at(i).s1<<" "<<it->second.at(i).s2<<endl;
	// 	}
	// }

	// for(map<string, vector<pinfo>>::iterator it = pairmap3.begin(); it!=pairmap3.end(); ++it){
	// 	cout<< it->first << ": " <<it->second.size()<<endl;
	// 	for(int i=0; i<it->second.size(); i++){
	// 		cout<<it->second.at(i).did<<": "<<it->second.at(i).s1<<" "<<it->second.at(i).s2<<endl;
	// 	}
	// }

	// for(map<string, vector<pinfo>>::iterator it = pairmap4.begin(); it!=pairmap4.end(); ++it){
	// 	cout<< it->first << ": " <<it->second.size()<<endl;
	// 	for(int i=0; i<it->second.size(); i++){
	// 		cout<<it->second.at(i).did<<": "<<it->second.at(i).s1<<" "<<it->second.at(i).s2<<endl;
	// 	}
	// }

	// for(map<string, vector<pinfo>>::iterator it = pairmap5.begin(); it!=pairmap5.end(); ++it){
	// 	cout<< it->first << ": " <<it->second.size()<<endl;
	// 	for(int i=0; i<it->second.size(); i++){
	// 		cout<<it->second.at(i).did<<": "<<it->second.at(i).s1<<" "<<it->second.at(i).s2<<endl;
	// 	}
	// }


	/*load the singlelist into map*/
	for (int i=0; i<lps.size(); ++i){

		// const int l_sqrt = sqrt(lps[i]->unpadded_list_length);
		cout<<lps[i]->term<<": "<<lps[i]->unpadded_list_length<<endl;
		// cout<<lps[i]->term<<" "<<lps[i]->unpadded_list_length<<endl;
		const string term = lps[i]->term;
		int term_id = Reader->term_map[term];
		// cout<<term<<" "<<term_id<<endl;
	 	RawIndexList Rlist = Reader->load_raw_list(term,term_id);

	 	vector<sinfo> t_list;

	 	for(int h = 0; h < lps[i]->unpadded_list_length; h++){
	 		sinfo t_p;
	 		t_p.did  = Rlist.doc_ids.at(h);
	 		// t_p.freq = Rlist.freq_s.at(h);
	 		t_p.score = Rlist.scores.at(h);
	 		t_list.push_back(t_p);
	 	}

	 	sort(t_list.begin(), t_list.end(), myfunc);

	 	t_list.resize(1000);// resize according to Configuration file

	 	// cout<<t_list.size()<<endl;

	 	for(int j = 0; j < t_list.size(); j++){

	 		if(t_list.at(j).did <= 10000000)
	 			temps_v1.push_back(t_list.at(j));
	 		if(t_list.at(j).did > 10000000 && t_list.at(j).did <= 20000000)
	 			temps_v2.push_back(t_list.at(j));
	 		if(t_list.at(j).did > 20000000 && t_list.at(j).did <= 30000000)
	 			temps_v3.push_back(t_list.at(j));
	 		if(t_list.at(j).did > 30000000 && t_list.at(j).did <= 40000000)
	 			temps_v4.push_back(t_list.at(j));
	 		if(t_list.at(j).did > 40000000 && t_list.at(j).did < CONSTS::MAXD)
	 			temps_v5.push_back(t_list.at(j));

	 	}

	 	temps_v.insert(temps_v.end(),temps_v1.begin(),temps_v1.end());
	 	temps_v.insert(temps_v.end(),temps_v2.begin(),temps_v2.end());
	 	temps_v.insert(temps_v.end(),temps_v3.begin(),temps_v3.end());
	 	temps_v.insert(temps_v.end(),temps_v4.begin(),temps_v4.end());
	 	temps_v.insert(temps_v.end(),temps_v5.begin(),temps_v5.end());

	 	// for(int i = 0; i<temps_v.size(); i++){
	 	// 	cout<<temps_v.at(i).did<<endl;
	 	// }

	 // singlemap1[term] = temps_v1;
  	// 	singlemap2[term] = temps_v2;
  	// 	singlemap3[term] = temps_v3;
  	// 	singlemap4[term] = temps_v4;
  	// 	singlemap5[term] = temps_v5;

	 	singlemap[term] = temps_v;
	 	temps_v1.clear();
  		temps_v2.clear();
  		temps_v3.clear();
  		temps_v4.clear();
  		temps_v5.clear();
  		temps_v.clear();

	}

	 //  	for(map<string, vector<sinfo>>::iterator it = singlemap1.begin(); it!=singlemap1.end(); ++it){
		// 	cout<< it->first << ": " <<it->second.size()<<endl;
		// 	// for(int i=0; i<it->second.size(); i++){
		// 	// 	cout<<it->second.at(i).did<<": "<<it->second.at(i).score<<endl;
		// 	// }
		// }

		// for(map<string, vector<sinfo>>::iterator it = singlemap2.begin(); it!=singlemap2.end(); ++it){
		// 	cout<< it->first << ": " <<it->second.size()<<endl;
		// 	// for(int i=0; i<it->second.size(); i++){
		// 	// 	cout<<it->second.at(i).did<<": "<<it->second.at(i).score<<endl;
		// 	// }
		// }

		// for(map<string, vector<sinfo>>::iterator it = singlemap3.begin(); it!=singlemap3.end(); ++it){
		// 	cout<< it->first << ": " <<it->second.size()<<endl;
		// 	// for(int i=0; i<it->second.size(); i++){
		// 	// 	cout<<it->second.at(i).did<<": "<<it->second.at(i).score<<endl;
		// 	// }
		// }

		// for(map<string, vector<sinfo>>::iterator it = singlemap4.begin(); it!=singlemap4.end(); ++it){
		// 	cout<< it->first << ": " <<it->second.size()<<endl;
		// 	// for(int i=0; i<it->second.size(); i++){
		// 	// 	cout<<it->second.at(i).did<<": "<<it->second.at(i).score<<endl;
		// 	// }
		// }

		// for(map<string, vector<sinfo>>::iterator it = singlemap5.begin(); it!=singlemap5.end(); ++it){
		// 	cout<< it->first << ": " <<it->second.size()<<endl;
		// 	// for(int i=0; i<it->second.size(); i++){
		// 	// 	cout<<it->second.at(i).did<<": "<<it->second.at(i).score<<endl;
		// 	// }
		// }

		// int size = 0;
		// for(map<string, vector<sinfo>>::iterator it = singlemap1.begin(); it!=singlemap1.end(); ++it){
		// 	cout<< it->first << ": " <<it->second.size()<<endl;
		// 	size = size + it->second.size();
		// }

		short offset = 0;
		short elem = 0;
		hashTable *ht;
		ht = initHash(100, 0); //table size 3079
		vector<int> didresults;
		vector<float> scoreresults;
		vector<string> kbitsresults;


		p.start(CONSTS::ALLQS);

		/*for the single lists hashing*/
		for(map<string, vector<sinfo>>::iterator it = singlemap.begin(); it!=singlemap.end(); ++it){
			for(int i=0; i<it->second.size(); i++){
				// cout<<it->second.at(i).did<<": "<<it->second.at(i).score<<endl;
				int pos = insertHash(ht, it->second.at(i).did, elem, 0, didresults);
    			if (pos){ //key already existed
     				 printf("Key %d already exists!\n", GETKEY(ht->table[pos-1], didresults));
      				 // score[ht->table[pos-1]]++;

   					 }else{ //successfully inserted
      					 // array[offset] = list[i];
   					 	// didresults.at[offset] = it->second.at(i).did;
   					 	didresults.push_back(it->second.at(i).did);
     					 // score[offset]++;
     			 		elem = ++ offset;
   					 	// printf("successfully inserted!\n");
   				 }
			}
		}

		// for(int i=0; i<didresults.size(); i++){

		// 	cout<<didresults.at(i)<<endl;
		// }
		/*singlelist*/

		offset = didresults.size();
		elem = offset;

		/*for the pair lists hashing*/
		for(map<string, vector<pinfo>>::iterator it = pairmap.begin(); it!=pairmap.end(); ++it){
			for(int i=0; i<it->second.size(); i++){
				// cout<<it->second.at(i).did<<": "<<it->second.at(i).s1<<" "<<it->second.at(i).s2<<endl;
				int pos = insertHash(ht, it->second.at(i).did, elem, 0, didresults);
    			if (pos){ //key already existed
     				 printf("Key %d already exists!\n", GETKEY(ht->table[pos-1], didresults));
      				 // score[ht->table[pos-1]]++;

   					 }else{ //successfully inserted
      					 // array[offset] = list[i];
   					 	// didresults.at[offset] = it->second.at(i).did;
   					 	didresults.push_back(it->second.at(i).did);
     					 // score[offset]++;
     			 		elem = ++ offset;
   					 	// printf("successfully inserted!\n");
   				 }
			}
		}
		/*pairlist*/

		p.end(CONSTS::ALLQS);


		cout<<"final did #: "<<didresults.size()<<endl;
		
		for(int i=0; i<didresults.size(); i++){
			cout<<didresults.at(i)<<endl;
		}

}
