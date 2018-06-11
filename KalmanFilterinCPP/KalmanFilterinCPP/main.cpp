#include <iostream>
#include <iomanip>
#include "CMatrix.h"
#include "CKalmanFilter.h"
#include "CBandpassFilters.h"
#include "CNonlinearEnergyOperator.h"
#include "CMeanFilter.h"
#include "CChatterDetection.h"
#include "CChatterDetectionSystem.h"
#include <fstream>
#include <string>
#include<time.h> // for clock

using namespace std;
double* ReadMeasurements(char* filePath);

int main(void)
{
	int N;
	char* filePath;
	int numOfData;
	double Ts;
	double samplingPeriod; //Ts calculated by DAQ

	int testNum = 1;
	switch (testNum)
	{
	case(1): N = 2;
			 filePath = "..\\CNCData.txt";
			 numOfData = 150016;
			 Ts = (double)1 / 10000;
			 samplingPeriod = (double)1 / 10000; //Ts calculated by DAQ
			 break;
	case(2): N = 2;
			 filePath = "..\\testData.txt";
			 numOfData = 150016;
			 Ts = (double)1 / 10000;
			 samplingPeriod = (double)1 / 10000; //Ts calculated by DAQ
			 break;
	default:
		N = 20;
		filePath = "..\\Acc2.txt";
		samplingPeriod = (double)1 / 25600; //Ts calculated by DAQ
		numOfData = 152788;
		Ts = (double)1 / 25600;
		break;
	}

	double spindleSpeed = 200 * 2 * PI; // [rad/s]
	double lamda = 1E-6;
	double R = 8.539750000056574E-6;

	double* measurements = ReadMeasurements(filePath);
	double startFreq = 200; //Hz
	double stopFreq = 4900;
	double TPE = 400; //Hz Original: 400Hz
	double stopBandMag = 5; //dB Original: 5
	double rippleMag = 2; //dB Original: 2 
	double integrationFactor = -1; // factor is 0 or 1 or -1
	double energyThreshold = 1E-8;
	double energyRatioLimit = 0.2;
	double dt = 0.1; // Original 0.1
	clock_t clockBegin, clockEnd;

	//CKalmanFilter myKalman(N, spindleSpeed, samplingPeriod, lamda, R);
	//CBandpassFilters myBandpassFilters(startFreq, stopFreq, TPE, dt, stopBandMag, rippleMag, Ts);
	//cout << "----------------------Numerator-----------------\n";
	//myBandpassFilters.Numerator->PrintMatrix();
	//cout << "----------------------Denominator------------------\n";
	//myBandpassFilters.Denominator->PrintMatrix();
	//cout << "----------------------Order-------------------------\n";
	//cout << myBandpassFilters.Order << "\n";
	//CNonlinearEnergyOperator myNEO(samplingPeriod, TPE, myBandpassFilters);
	//CMeanFilter myMeanFilter(spindleSpeed, samplingPeriod, myBandpassFilters.NumberOfFilters);
	//CChatterDetection myChatterDetection(myMeanFilter.ndMean, energyThreshold, energyRatioLimit, integrationFactor);

	CChatterDetectionSystem myChatterDetectionSys(N, spindleSpeed, samplingPeriod, lamda, R, 
							startFreq, stopFreq, TPE, dt, stopBandMag, rippleMag, 
							energyThreshold, energyRatioLimit, integrationFactor);

	//FILE* OutFilePointer = fopen("..\\Kalman Output.txt", "w");
	FILE* BandpassPointer = fopen("..\\Bandpass Filter Output.txt", "w");
	FILE* FrequencyPointer = fopen("..\\Frequency Output.txt", "w");
	FILE* AmplitudePointer = fopen("..\\Amplitude Output.txt", "w");
	FILE* FilteredFreqPointer = fopen("..\\Filtered Freq Output.txt", "w");
	FILE* FilteredAmpPointer = fopen("..\\Filtered Amp Output.txt", "w");
	FILE* ChatterDetectionPointer = fopen("..\\Chatter Detection Output.txt", "w");
		
	clockBegin = clock();
	for (int i = 0; i < numOfData; i++)
	//for (int i = 0; i < 10000; i++)
	{
		//double estimation = myKalman.RunKalman(measurements[i]);
		//myBandpassFilters.RunBandpassFilters(myBandpassFilters.CalculateInput(measurements[i], estimation));
		//myNEO.RunNEO(myBandpassFilters.BandpassOutputs);
		//myMeanFilter.RunMeanFilter(myNEO.freq, myNEO.amp);
		//myChatterDetection.RunChatterDetection(myKalman.PeriodicAmp, spindleSpeed, myMeanFilter.FilteredAmp, myMeanFilter.FilteredFreq);

		myChatterDetectionSys.Run(measurements[i], spindleSpeed);

		//for (int j = 0; j < myBandpassFilters.NumberOfFilters; j++)
		//{
		//	fprintf(BandpassPointer, "%.8f	", myBandpassFilters.BandpassOutputs->Content[j]);
		//	fprintf(FrequencyPointer, "%.8f	", myNEO.freq->Content[j] / 2 / PI); //[Hz]
		//	fprintf(AmplitudePointer, "%.8f	", myNEO.amp->Content[j]);
		//	fprintf(FilteredFreqPointer, "%.8f	", myMeanFilter.FilteredFreq->Content[j]);//[rad/s]
		//	fprintf(FilteredAmpPointer, "%.8f	", myMeanFilter.FilteredAmp->Content[j]);
		//}

		//for (int j = 0; j < myChatterDetectionSys.BandpassFilters->NumberOfFilters; j++)
		//{
		//	fprintf(BandpassPointer, "%.8f	", myChatterDetectionSys.BandpassFilters->BandpassOutputs->Content[j]);
		//	fprintf(FrequencyPointer, "%.8f	", myChatterDetectionSys.NEO->freq->Content[j] / 2 / PI); //[Hz]
		//	fprintf(AmplitudePointer, "%.8f	", myChatterDetectionSys.NEO->amp->Content[j]);
		//	fprintf(FilteredFreqPointer, "%.8f	", myChatterDetectionSys.MeanFilter->FilteredFreq->Content[j]);//[rad/s]
		//	fprintf(FilteredAmpPointer, "%.8f	", myChatterDetectionSys.MeanFilter->FilteredAmp->Content[j]);
		//}
		//fprintf(BandpassPointer, "\n");
		//fprintf(FrequencyPointer, "\n");
		//fprintf(AmplitudePointer, "\n");
		//fprintf(FilteredFreqPointer, "\n");
		//fprintf(FilteredAmpPointer, "\n");
		fprintf(ChatterDetectionPointer, "%d\n", myChatterDetectionSys.ChatterDetection->ChatterDetected);
	}		
	clockEnd = clock();
	cout << "done!\n";
	cout << "Each loop takes " << setprecision(10) << ((double)clockEnd - (double)clockBegin) / numOfData << "ms" << endl;
	//cout << setprecision(10) << (double)clockEnd << endl;

	fclose(BandpassPointer);
	fclose(FrequencyPointer);
	fclose(AmplitudePointer);
	fclose(FilteredFreqPointer);
	fclose(FilteredAmpPointer);
	fclose(ChatterDetectionPointer);

	system("pause");
	return 0;
}

//double* ReadMeasurements(char* filePath)
//{
//	double currentData;
//	string line;
//	int dataNumber = 0;
//	int i;
//
//	ifstream myDataFile(filePath);
//
//	if (myDataFile.is_open())
//	{
//		while (getline(myDataFile, line))
//		{
//			dataNumber++;
//		}
//	}
//	else
//		cout << "Unable to open file" << "\n";
//
//	myDataFile.close();
//
//	myDataFile.open(filePath);
//
//	double * measurement = new double[dataNumber];
//
//	while (getline(myDataFile, line))
//	{
//		measurement[i] = atof(line);
//	}
//
//	return measurement;
//}

double* ReadMeasurements(char* filePath)
{
	char data;
	int dataNumber = 0;
	int i;
	FILE* fp;
	fp = fopen(filePath, "r");

	if (fp == NULL)
	{
		printf("Could not open file!\n");
		return NULL;
	}

	while ((data = fgetc(fp)) != EOF)
	{
		if (data == '\n')
			dataNumber++;
	}

	dataNumber++; // ++1 since no '\n'(return line character) for the last data

	fclose(fp);
	//FILE* fp2 = fopen(filePath, "r");
	fp = fopen(filePath, "r");
	double* measurement = (double*)calloc(dataNumber, sizeof(double));

	char* test = (char *)malloc(100 * sizeof(char)); // allocate memory for a string of 100 characters

	for (i = 0; i < dataNumber; i++)
	{
		fscanf(fp, "%s", test);
		measurement[i] = atof(test);
	}

	fclose(fp);
	return measurement;
}