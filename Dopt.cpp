/*
 * Dopt.cpp
 *
 *  Created on: July 29, 2014
 *      Author: Qi
 */

#include "Dopt.h"
#include "globals.h"
#include "profiling.h"
#include "utils.h"
#include <cstring>
#include <string>
#include <fstream>
#include <algorithm>
#include <set>
using namespace std;

void Dopt::operator() (int qn, pairlists& pls, lptrArray& lps, const int topK, QpResult* res) {

	// vector<map<ulong, int>> pimv;
	// vector<map<ulong, int>> timv;
	const int p_l = pls.lengths.size();
	const int t_l = lps.size();
	const int a_l = pls.lengths.size() + lps.size();
	// cout<<"singleterms: "<<lps.size()<<endl;
	// cout<<"pairs: "<<pls.lengths.size()<<endl;
	// cout<<"a_l: "<<a_l<<endl;
	int p_e[10][p_l]; 
	int t_e[10][t_l]; 

	const int buffer_size = 1024*100;
	string queryline;
	ulong results[10];
	int known_l[10];
	for(int i = 0; i < 10; i++){
		known_l[i] = 0;
	}

    char * docid_s;
	// int docid;
	int count = 0;
	int test;

	/*---Get the top 10 results from Complicated ranking func----*/

	ifstream input_stream;
	// input_stream.open("../../../Dropbox/WSDM_Index_Script/result_log");  //for or
	// input_stream.open("/data/qw376/WSDM_Index_Script/500_cr");  //for cr
	input_stream.open(CONSTS::cr_results.c_str());  //for cr


	for(int i=0; i < qn-1; ++i){
        input_stream.ignore(numeric_limits<streamsize>::max(),'\n');
    }

	getline(input_stream, queryline);

	char * ql = new char [queryline.length()+1];
    std::strcpy (ql, queryline.c_str());

	// cout<<queryline<<endl;


	docid_s = strtok (ql," ");

	while (docid_s != NULL)
  	{
    	// printf ("%s\n",docid_s);
 		results[count] = atol(docid_s);   	
 		count++;
    	docid_s = strtok (NULL, " ,.-");
    	if(count == 10)
    		break;
  	}

  	input_stream.close();

  	// for(int b=0; b<count; b++){
  	// 	cout<<results[b]<<endl;
  	// }
  	/*-----------------------*/


  	/*---load the pair lists into map and outputs the depth----*/
  	string did_s;
  	ulong didl;
  	// string sc_s;
  	// float sc;

  	if(pls.lengths.size()>0){

  	for(int k = 0; k < pls.lengths.size(); k++){

  	map<ulong, int> t_map;

  	ifstream pair_stream;
  	string pair_dir = CONSTS::pair_index + pls.pairnames.at(k);
  	// string pair_dir = "/data/qw376/pair_index/" + pls.pairnames.at(k);
  	pair_stream.open(pair_dir.c_str());

  	count = 0;

  	// cout<<pls.pairnames.at(k)<<endl;

  	while(getline(pair_stream, queryline)){

  		// cout<<queryline<<endl;
  		string::iterator itr = queryline.begin();
		string::iterator start = itr;

   		while(itr != queryline.end() && !isspace(*itr)){
			++itr;
		}

		did_s = string(start, itr);
		didl = atol(did_s.c_str());


		// start = itr+1;
  // 	  	itr++;
  //   	while(itr != queryline.end() && !isspace(*itr)){
		// 	++itr;
		// }

		// sc_s = string(start, itr);
		// sc = atof(s1_s.c_str());

  		count ++;		
		t_map[didl] = count;
  		if (count == pls.lengths.at(k))
  			break;

  	}
  		// pimv.push_back(t_map);
		// for (map<ulong,scores>::iterator it=t_map.begin(); it!=t_map.end(); ++it)
  //  			 cout << it->first << " => " << it->second.s1 << endl;

  		// cout<<t_map.size()<<endl;

  		for(int i = 0; i < 10; i++){

  			map<ulong, int>:: iterator it;
			it = t_map.find(results[i]);
			if(it != t_map.end()){
				// cout<<results[i]<<" is at depth: "<<it->second<<" of "<<pls.pairnames.at(k)<<endl;
				p_e[i][k] = it->second;

				
			}else{
				// cout<<results[i]<<" is at depth: - of "<<pls.pairnames.at(k)<<endl;
				p_e[i][k] = 0;
			}

  		}


  	  	pair_stream.close();
	}

	}//if pair list size is not 0
	/*-----------------------*/

	/*---load the single term lists into map and outputs the depth (From ASCII lists)----*/
	// for(int k = 0; k<lps.size(); k++){

	// 	map<ulong, int> t_map;
	// 	ifstream index_score;
	// 	string score_ind_dir = "/home/qi/Dropbox/score_index/" + lps[k]->term;
	// 	index_score.open(score_ind_dir.c_str(), ofstream::app);


	// 	count = 0;
	// 	cout<<lps[k]->term<<endl;
	// 	while(getline(index_score, queryline)){

	// 		string::iterator itr = queryline.begin();
	// 		string::iterator start = itr;

 //   			while(itr != queryline.end() && !isspace(*itr)){
	// 			++itr;
	// 		}

	// 		did_s = string(start, itr);
	// 		didl = atol(did_s.c_str());

	// 		// t_a[count] = didl;
	// 		count++;
	// 		t_map[didl] = count;
	// 	}

	// 	cout<<count<<endl;
	// 	cout<<t_map.size()<<endl;
	// 	cout<<lps[k]->unpadded_list_length<<endl;

 //  		for(int i = 0; i < 10; i++){

 //  			map<ulong, int>:: iterator it;
	// 		it = t_map.find(results[i]);
	// 		if(it != t_map.end()){
	// 			cout<<results[i]<<" is at depth: "<<it->second<<" of "<<lps[k]->term<<endl;
				
	// 		}else{
	// 			cout<<results[i]<<" is at depth: -"<<" of "<<lps[k]->term<<endl;
	// 		}

 //  		}
	// }
	/*-----*/

	/*---load the single term lists into map and outputs the depth (From Binary lists)----*/
	for(int k = 0; k<lps.size(); k++){

		ifstream index_score;
		// string score_ind_dir = "/home/qi/Dropbox/score_index_binary_m/" + lps[k]->term;
		// string score_ind_dir = "/data/qw376/score_index_binary/" + lps[k]->term;
		string score_ind_dir = CONSTS::score_index_binary.c_str() + lps[k]->term;

		const int size = lps[k] -> unpadded_list_length;

		// cout<<lps[k] -> term<<": "<<size<<endl;

		/*if list length is less than a block*/
		if (size <= buffer_size){
			index_score.open(score_ind_dir.c_str(), ios::binary);
			unsigned int t_a[size];
			index_score.read((char*)t_a, sizeof(unsigned int)*size);

  		for(int i = 0; i < 10; i++){
  				int j;
  				for(j = 0; j<size; j++){
  					if (results[i] == t_a[j])
  						break;
  				}

				if(j!=size){
					// cout<<results[i]<<" is at depth: "<<j+1<<" of "<<lps[k]->term<<endl;
					t_e[i][k] = j+1;
				
				}else{
					// cout<<results[i]<<" is at depth: -"<<" of "<<lps[k]->term<<endl;
					t_e[i][k] = 0;
				}

  			}
  			index_score.close();
  		}
  		/*----*/

  		/*if the list length is more than a block*/
  		if (size >= buffer_size){

  			const int block =  size / buffer_size;
			const int rem = size % buffer_size;
			// cout<<"block: "<<block<<endl;
			// cout<<"rem: "<<rem<<endl;
			int depths[10];
			for(int m = 0; m < 10; m++){
				depths[m] = -1;
			}

			//interate the 10 results
  			for(int i = 0; i < 10; i++){

  				index_score.open(score_ind_dir.c_str(), ios::binary);

  				int depth = -1;

  				//iterate all the blocks
  				for (int h = 0; h < block; h++){

  					// cout<<"p1"<<endl;
  					unsigned int t_a[buffer_size];
  					// cout<<"p2"<<endl;
  					index_score.read((char*)t_a, sizeof(unsigned int)*buffer_size);

  					int j;
  					for (j = 0; j < buffer_size; j++){
  						// cout<<t_a[j]<<": "<<depth<<endl;
  						if(results[i] == t_a[j]){
  							depth = j + buffer_size * h + 1;
  							break;
  						}
  					}

  					if(depth > 0) {
  						depths[i] = depth;
  						break;
  					}
  					// int test;
  					// cin>>test;
  				}

  				//iterate remaining parts
  				if(depth == -1){

  				   unsigned int t_a[rem];

  				   index_score.read((char*)t_a, sizeof(unsigned int)*rem);

  				   int j;
  					for (j = 0; j < rem; j++){
  						if(results[i] == t_a[j]){
  							depth = j + buffer_size * block;
  							break;
  						}
  					}

  					if(depth > 0) {
  						depths[i] = depth;
  					}
  				}

  				index_score.close();
  			}

  			for(int m = 0; m < 10; m++){
				if(depths[m] == -1){
					// cout<<results[m]<<" is at depth: -"<<" of "<<lps[k]->term<<endl;
					t_e[m][k] = 0;
				}else{
					// cout<<results[m]<<" is at depth: "<<depths[m]<<" of "<<lps[k]->term<<endl;
					t_e[m][k] = depths[m] + 1;
				}
			}
  		

  		}

  		/*----*/



	}
	/*---------*/


	/*verify the lists info are correct*/
	// ifstream index_score;
	// ofstream index_veri;
	// for(int k = 0; k<lps.size(); k++){
	// 	if(lps[k]->unpadded_list_length<buffer_size){
	// 	const int size = lps[k]->unpadded_list_length;
	// 	unsigned int t_a[size];

	// 	cout<<"term: "<<lps[k]->term<<endl;

	// 	const string score_ind_dir = "/home/qi/score_index_binary/" + lps[k]->term;
	// 	FILE * index_score = fopen ( score_ind_dir.c_str() , "rb" );
	// 	if(index_score==NULL){
	// 		cout<<"file broken"<<endl;
	// 	}
	// 	cout<<fread (t_a, sizeof(unsigned int), size, index_score)<<endl;

	// 	const string score_ind_dir_v = "/home/qi/Dropbox/score_index_veri/" + lps[k]->term;
	// 	index_veri.open(score_ind_dir_v.c_str(), ofstream::app);
	// 	for(int n = 0; n < size; n++){
	// 		index_veri<<t_a[n]<<endl;
	// 	}
	// 	index_veri.close();
	//   }
	// }
	/*---------------------*/


	ofstream index_depth;
	// const string dir = "/data/qw376/experiments/depth_index";
	index_depth.open(CONSTS::index_depth.c_str(), ofstream::app);
	for(int i = 0; i < 10; i++){
		for(int n = 0; n < t_l; n++){
			index_depth<<lps[n]->term<<" "<<t_e[i][n]<<" ";
		}
		for(int n = 0; n < p_l; n++){
			index_depth<<pls.pairnames.at(n)<<" "<<p_e[i][n]<<" ";
		}
		index_depth<<endl;
	}
	index_depth.close();

}
