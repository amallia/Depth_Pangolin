/*
 * algo_init.h
 *
 *  Created on: Aug 13, 2014
 *      Author: Qi
 */

#ifndef ALGO_INIT_H_
#define ALGO_INIT_H_

#include "PostingOriented_BMW.h"
#include "CluewebReader.h"

class algo_init{
private:
	unsigned int* pages;

public:
	algo_init(unsigned int* pgs) : pages(pgs) {}
	void operator()(CluewebReader* Reader, int qn, pairlists& pls, lptrArray& lps, const int topK, QpResult* res, profilerC& p);
};

#endif /* ALGO_INIT_H_ */
