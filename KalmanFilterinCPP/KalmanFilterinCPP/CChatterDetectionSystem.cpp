#include "CChatterDetectionSystem.h"

CChatterDetectionSystem::CChatterDetectionSystem(int N, double spindleSpeed, double samplingPeriod, double lamda, double R, 
												 double startFreq, double stopFreq, double TPE, double dt, double stopBandMag, double rippleMag, 
												 double energyThreshold, double energyRatioLimit, double integrationFactor)
{
	KalmanFilter = new CKalmanFilter(N, spindleSpeed, samplingPeriod, lamda, R);
	BandpassFilters = new CBandpassFilters(startFreq, stopFreq, TPE, dt, stopBandMag, rippleMag, samplingPeriod);

	NEO = new CNonlinearEnergyOperator(samplingPeriod, TPE, *BandpassFilters);
	MeanFilter = new CMeanFilter(spindleSpeed, samplingPeriod, BandpassFilters->NumberOfFilters);
	ChatterDetection = new CChatterDetection(MeanFilter->ndMean, energyThreshold, energyRatioLimit, integrationFactor);
}

CChatterDetectionSystem::~CChatterDetectionSystem()
{
	delete KalmanFilter;
	delete BandpassFilters;
	delete NEO;
	delete MeanFilter;
	delete ChatterDetection;
}

void CChatterDetectionSystem::Run(double measurement, double spindleSpeed)
{
	double estimation = KalmanFilter->RunKalman(measurement);
	BandpassFilters->RunBandpassFilters(BandpassFilters->CalculateInput(measurement, estimation));
	NEO->RunNEO(BandpassFilters->BandpassOutputs);
	MeanFilter->RunMeanFilter(NEO->freq, NEO->amp);
	ChatterDetection->RunChatterDetection(KalmanFilter->PeriodicAmp, spindleSpeed, MeanFilter->FilteredAmp, MeanFilter->FilteredFreq);
}