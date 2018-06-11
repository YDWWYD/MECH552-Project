#include "CMeanFilter.h"

CMeanFilter::CMeanFilter(double spindleSpeed, double samplingPeriod, int numOfInput) 
{
	Counter = 0;
	ndMean = round(1 / (spindleSpeed / 2 / PI) / samplingPeriod);

	FilteredFreq = new CMatrix(numOfInput, 1);
	FilteredAmp = new CMatrix(numOfInput, 1);
	PrevAveAmp = new CMatrix(numOfInput, 1);
	PrevAveFreq = new CMatrix(numOfInput, 1);

	NumberOfInput = numOfInput;
}

CMeanFilter::~CMeanFilter(void)
{
	delete FilteredFreq;
	delete FilteredAmp;
	delete PrevAveFreq;
	delete PrevAveAmp;
}

void CMeanFilter::RunMeanFilter(CMatrix* freqIn, CMatrix* ampIn)
{
	if (Counter < ndMean)
	{
		if (Counter != 0)
		{
			for (int i = 0; i < NumberOfInput; i++)
			{
				FilteredFreq->Content[i] = (Counter - 1) / (Counter)*PrevAveFreq->Content[i] + 1 / (Counter)*freqIn->Content[i];
				FilteredAmp->Content[i] = (Counter - 1) / (Counter)*PrevAveAmp->Content[i] + 1 / (Counter)*ampIn->Content[i];
			}
		}
	}

	else
	{
		for (int i = 0; i < NumberOfInput; i++)
		{
			FilteredFreq->Content[i] = (ndMean - 1) / ndMean*PrevAveFreq->Content[i] + 1 / ndMean*freqIn->Content[i]; // In [rad/s]
			FilteredAmp->Content[i] = (ndMean - 1) / ndMean*PrevAveAmp->Content[i] + 1 / ndMean*ampIn->Content[i];
		}
	}

	// Update previous average frequency and amplitude
	for (int i = 0; i < NumberOfInput; i++)
	{
		PrevAveFreq->Content[i] = FilteredFreq->Content[i]; // in [rad/s]
		PrevAveAmp->Content[i] = FilteredAmp->Content[i];
	}

	//Update counter;
	Counter++;
}