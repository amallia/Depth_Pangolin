/*
 * DepthCal.h
 *
 *  Created on: July 22, 2014
 *      Author: Qi
 */

#ifndef DEPTHCAL_H_
#define DEPTHCAL_H_

#include "PostingOriented_BMW.h"

class DepthCal{
private:
	unsigned int* pages;

public:
	DepthCal(unsigned int* pgs) : pages(pgs) {}
	void operator()(int qn, pairlists& pls, lptrArray& lps, const int topK, QpResult* res);
};

#endif /* DEPTHCAL_H_ */
