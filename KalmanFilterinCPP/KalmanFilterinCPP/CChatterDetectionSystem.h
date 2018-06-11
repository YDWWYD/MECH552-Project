#ifndef __CChatterDetectionSystem_
#define __CChatterDetectionSystem_

#include "CKalmanFilter.h"
#include "CBandpassFilters.h"
#include "CNonlinearEnergyOperator.h"
#include "CMeanFilter.h"
#include "CChatterDetection.h"

class CChatterDetectionSystem
{
//public:
//	double PeriodicEnergy;
//	double ChatterEnergy;
//	double EnergyRatio;
//	double EnergyThreshold;
//	double EnergyRatioLimit;

public:
	CKalmanFilter* KalmanFilter;
	CBandpassFilters* BandpassFilters;
	CNonlinearEnergyOperator* NEO;
	CMeanFilter* MeanFilter;
	CChatterDetection* ChatterDetection;

public:
	CChatterDetectionSystem(int N, double spindleSpeed, double samplingPeriod, double lamda, double R, double startFreq, double stopFreq, double TPE, double dt, double stopBandMag, double rippleMag, double energyThreshold, double energyRatioLimit, double integrationFactor);
	~CChatterDetectionSystem();
	void Run(double measuremnt, double spindleSpeed);
};

#endif
