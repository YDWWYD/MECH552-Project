#include "CChatterDetection.h"

CChatterDetection::CChatterDetection(double ndMean, double energyThreshold, double energyRatioLimit, double integrationFactor)
{
	ChatterDetected = 0;
	PeriodicEnergy = 0;
	ChatterEnergy = 0;
	EnergyRatio = 0;
	EnergyThreshold = energyThreshold;
	EnergyRatioLimit = energyRatioLimit;

	PrevChatterDetected = 0;
	DetectionDelay = (int)round((ndMean + 1) / 2);
	IntegrationFactor = integrationFactor;
	DelayedPeriodicEnergy = new CMatrix(DetectionDelay, 1);
}

CChatterDetection::~CChatterDetection(void)
{
	delete DelayedPeriodicEnergy;
}

double CChatterDetection::CalculatePeriodicEnergy(CMatrix* amplitude, double spindleSpeed)
{
	double periodicEnergy = 0;

	for (int i = 0; i < amplitude->Row; i++)
	{
		double frequency = (i + 1)*spindleSpeed;
		if (frequency != 0)
			periodicEnergy = periodicEnergy + pow(amplitude->Content[i], 2)*pow(frequency, 2 * IntegrationFactor);
	}

	return periodicEnergy;
}

double CChatterDetection::CalculateChatterEnergy(CMatrix* amplitude, CMatrix* frequency)
{
	double chatterEnergy = 0;

	for (int i = 0; i < amplitude->Row; i++)
	{
		if (frequency->Content[i] != 0)
			chatterEnergy = chatterEnergy + pow(amplitude->Content[i], 2)*pow(frequency->Content[i], 2 * IntegrationFactor);
	}

	return chatterEnergy;
}

void CChatterDetection::RunChatterDetection(CMatrix* periodicAmplitude, double spindleSpeed, CMatrix* chatterAmp, CMatrix* chatterFreq)
{
	//int chatterDetected;

	PeriodicEnergy = CalculatePeriodicEnergy(periodicAmplitude, spindleSpeed);
	ChatterEnergy = CalculateChatterEnergy(chatterAmp, chatterFreq);

	EnergyRatio = ChatterEnergy / (ChatterEnergy + DelayedPeriodicEnergy->Content[0]);
	//fprintf(outputFilePtr, "%.25f\n", energyRatio);

	if (PrevChatterDetected == 0 && EnergyRatio > 0.5) // entry
		ChatterDetected = 0;
	else if (PrevChatterDetected == 1 && EnergyRatio > 0.1) // chatter continues
		ChatterDetected = 1;
	else if (EnergyRatio > EnergyRatioLimit && DelayedPeriodicEnergy->Content[0] > EnergyThreshold) // chatter just detected
		ChatterDetected = 1;
	else if (EnergyRatio < 0.1) //no chatter
		ChatterDetected = 0;
	else
		ChatterDetected = 0;

	// Update if chatter is detected by previous data
	PrevChatterDetected = ChatterDetected;

	// Update delayed periodic energy
	for (int i = 0; i < DetectionDelay - 1; i++)
	{
		DelayedPeriodicEnergy->Content[i] = DelayedPeriodicEnergy->Content[i + 1];
	}
	DelayedPeriodicEnergy->Content[DetectionDelay - 1] = PeriodicEnergy;

	//return chatterDetected;
}