#ifndef __CMeanFilter_
#define __CMeanFilter_

#include "CMatrix.h"
#include <cmath>
#define PI acos(-1.0)

class CMeanFilter
{
public:
	CMatrix* FilteredFreq;
	CMatrix* FilteredAmp;
	double ndMean;
private:
	CMatrix* PrevAveFreq;
	CMatrix* PrevAveAmp;
	double Counter; // double type or cast to double to avoid being converted to int in RunMeanFilter()
	int NumberOfInput;

public:
	CMeanFilter(double spindleSpeed , double samplingPeriod, int numOfInput);
	~CMeanFilter();
	void RunMeanFilter(CMatrix* freqInput, CMatrix* ampInput);
};

#endif
