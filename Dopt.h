/*
 * Dopt.h
 *
 *  Created on: July 29, 2014
 *      Author: Qi
 */

#ifndef DOPT_H_
#define DOPT_H_

#include "PostingOriented_BMW.h"

class Dopt{
private:
	unsigned int* pages;

public:
	Dopt(unsigned int* pgs) : pages(pgs) {}
	void operator()(int qn, pairlists& pls, lptrArray& lps, const int topK, QpResult* res);
};

#endif /* DOPT_H_ */
