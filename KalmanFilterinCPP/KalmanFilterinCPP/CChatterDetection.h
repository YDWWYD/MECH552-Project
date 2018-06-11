#ifndef __CChatterDetection_
#define __CChatterDetection_

#include "CMatrix.h"
#include <cmath>
#define PI acos(-1.0)

class CChatterDetection
{
public:
	int ChatterDetected;
	double PeriodicEnergy;
	double ChatterEnergy;
	double EnergyRatio;
	double EnergyThreshold;
	double EnergyRatioLimit;

private:
	int PrevChatterDetected;
	int DetectionDelay;
	double IntegrationFactor;
	CMatrix* DelayedPeriodicEnergy;

public:
	CChatterDetection(double ndMean, double energyThreshold, double energyRatioLimit, double integrationFactor);
	~CChatterDetection();
	double CalculatePeriodicEnergy(CMatrix* amplitude, double spindleSpeed);
	double CalculateChatterEnergy(CMatrix* amplitude, CMatrix* frequency);
	void RunChatterDetection(CMatrix* periodicAmplitude, double spindleSpeed, CMatrix* chatterAmp, CMatrix* chatterFreq);
};

#endif
