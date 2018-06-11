#include "buttordBandpassOnly.h"
#include "butterBandpassOnly.h"
#include "Matrix.h"
#include "MatrixOperations.h"
#include "Math.h"
#include "Kalman.h"

typedef struct BandpassCoefficients
{
	Matrix Numerator;
	Matrix Denominator;
	int Order;
	int NumOfFilters;
	int StartBand;
}BandpassCoefficients;

typedef struct {
	Matrix Output;
	Matrix EDelayed;
} BandpassOutput;

typedef struct {
	Matrix Frequency;
	Matrix Amplitude;
	Matrix SDelayed;
} NEOOutput;

typedef struct {
	Matrix FilteredFreq;
	Matrix FilteredAmp;
} RecursiveMeanOutput;

BandpassCoefficients BandpassFilterCoefficients(double StartFreq, double StopFreq, double TPE, double dt, double StopBandMag, double RippleMag, double Ts);
Matrix GetBandpassInput(double measurement, double estimation, int numOfFilters);
BandpassOutput BandpassFilterAlgorithm(Matrix input, Matrix A, Matrix B, Matrix EDelayed);
Matrix CalculateLag(double samplingPeriod, double TPE, BandpassCoefficients bandpassFilterCoefficient);
double Max(Matrix matrix);
Matrix CreateHarmonicSpindleSpeed(int N, double spindleSpeed);
NEOOutput NEOAlgorithm(Matrix input, Matrix D, Matrix s_delay, double samplingPeriod, Matrix prevFreq, Matrix prevAmp);
RecursiveMeanOutput RecursiveMeanFilter(Matrix freqIn, Matrix ampIn, double ndMean, double* counter, Matrix prevAveFreq, Matrix prevAveAmp, BandpassCoefficients bandpassCoefficients);
double CalculateVibrationEnergy(Matrix amplitude, Matrix frequency, double integrationFactor);
int ChatterDetection(double periodicEnergy, double chatterEnergy, double energyThreshold, double energyRatioLimit, int* prevChatterDetected, Matrix delayedPeriodicEnergy, FILE* outputFilePtr);

int main(void)
{
	double startFreq = 200; //Hz
	double stopFreq = 4900; 
	double TPE = 400; //Hz Original: 400Hz
	double stopBandMag = 5; //dB Original: 5
	double rippleMag = 2; //dB Original: 2 
	double Ts = (double)1 / 25600;
	//printf("Ts is %f\n", Ts);
	double dt = 0.1; // Original 0.1

	char* filePath = "..\\Acc2.txt";
	int N = 20;
	double spindleSpeed = 200 * 2 * PI;
	//double samplingPeriod = 3.90625E-5;
	double samplingPeriod = (double)1 / 25600;
	double lamda = 1E-6;

	Matrix Phi = CreatePhiMatrix(N, spindleSpeed, samplingPeriod);
	Matrix H = CreateHMatrix(N);

	Matrix measurements = ReadMeasurements(filePath);
	Matrix R = CalculateVariance(measurements, 100, 700);

	Matrix Q = CreateQMatrix(N, lamda, R);
	Matrix q0 = CreateInitial_q(N);
	Matrix P0 = CreateInitialP(N);

	Matrix PPost = P0;
	Matrix qPost = q0;

	Matrix PhiTranspose = Transpose(Phi);
	Matrix identity2N = IdentityMatrix(P0.row); // Identity [2N * 2N]
	Matrix HTranspose = Transpose(H);
	Matrix SpEst = CreateEmptyMatrix(measurements.row, 1);
	Matrix periodicAmp = CreateEmptyMatrix(measurements.row, q0.row / 2); // S_hat_p,n

	//--------------Bandpass Filter Initialization--------------
	BandpassCoefficients bandpassCoefficients = BandpassFilterCoefficients(startFreq, stopFreq, TPE, dt, stopBandMag, rippleMag, Ts);
	printf("----------Denominator----------------\n");
	PrintMatrix(bandpassCoefficients.Denominator);
	printf("----------Numerator------------------\n");
	PrintMatrix(bandpassCoefficients.Numerator);
	printf("----------Order of Filter----------------\n");
	printf("The order of filter is %d\n", bandpassCoefficients.Order);
	printf("----------Num of Filters------------------\n");
	printf("The number of filters is %d\n", bandpassCoefficients.NumOfFilters);

	//Matrix A = ExtractRow(bandpassCoefficients.Denominator, 1); // A: [1 x 2n]
	//Matrix B = ExtractRow(bandpassCoefficients.Numerator, 1); // B: [1 x 2n+1]
	//Matrix EDelayed = CreateEmptyMatrix(1, bandpassCoefficients.Order * 2); // Delayed E: [1 x 2n]
	Matrix A = bandpassCoefficients.Denominator; // A: [m x 2n]
	Matrix B = bandpassCoefficients.Numerator; // B: [m x 2n+1]
	Matrix EDelayed = CreateEmptyMatrix(bandpassCoefficients.NumOfFilters, bandpassCoefficients.Order * 2); // Delayed E: [m x 2n]
	Matrix filteredOutput = CreateEmptyMatrix(measurements.row, bandpassCoefficients.NumOfFilters);
	//Matrix BPInput = ReadMeasurements("..\\BPIn.txt");

	//-------------------Non-linear Energy Operation Initilization---------------------------
	Matrix D = CalculateLag(samplingPeriod, TPE, bandpassCoefficients);
	//PrintMatrix(D);
	int maxDelay = (int)Max(D);
	Matrix SDelayed = CreateEmptyMatrix(bandpassCoefficients.NumOfFilters, 4 * maxDelay);
	Matrix freqOutput = CreateEmptyMatrix(measurements.row, bandpassCoefficients.NumOfFilters);
	Matrix ampOutput = CreateEmptyMatrix(measurements.row, bandpassCoefficients.NumOfFilters);
	Matrix prevFreq = CreateEmptyMatrix(measurements.row, bandpassCoefficients.NumOfFilters);
	Matrix prevAmp = CreateEmptyMatrix(measurements.row, bandpassCoefficients.NumOfFilters);
	//Matrix singleNEOInput = ReadMeasurements("..\\NEOInput.txt");

	//-------------------Recursive Mean Filter Initilization---------------------------
	double* counter  = (double*)calloc(1, sizeof(double));
	double ndMean = round(1/ (spindleSpeed/2/PI)/samplingPeriod);
	Matrix filteredFreqOutput = CreateEmptyMatrix(measurements.row, bandpassCoefficients.NumOfFilters);
	Matrix filteredAmpOutput = CreateEmptyMatrix(measurements.row, bandpassCoefficients.NumOfFilters);
	Matrix prevAveFreq = CreateEmptyMatrix(bandpassCoefficients.NumOfFilters, 1);
	Matrix prevAveAmp = CreateEmptyMatrix(bandpassCoefficients.NumOfFilters, 1);

	// ------------------Energy ratio and chatter dection initilization-----------------
	Matrix harmonicSpindleSpeed = CreateHarmonicSpindleSpeed(N, spindleSpeed);
	double periodicEnergy, chatterEnergy;
	double integrationFactor = -1; // factor is 0 or 1 or -1
	FILE* ChatterEnergyPointer = fopen("..\\Chatter Energy Output.txt", "w");
	FILE* PeriodicEnergyPointer = fopen("..\\Periodic Energy Output.txt", "w");
	int chatterDetected;
	int* prevChatterDetected = (int*) calloc(1, sizeof(int));
	double energyThreshold = 1E-8;
	double energyRatioLimit = 0.2;
	double detectionDelay = round((ndMean + 1) / 2);
	Matrix delayedPeriodicEnergy = CreateEmptyMatrix(1, detectionDelay);
	FILE* ChatterDetectedPointer = fopen("..\\Chatter Detection Output.txt", "w");
	FILE* EnergyRatioPointer = fopen("..\\Energy Ratio Output.txt", "w");

	for (int i = 0; i < measurements.row; i++)
	//for (int i = 0; i < 20000; i++)
	{
		//----------------Kalman Filter-----------------
		KalmanOutput KalmanOut = KalmanAlgorithm(measurements.content[i][0], Phi, PhiTranspose, H, HTranspose, identity2N, R, Q, PPost, qPost);
		FreeMatrixMemory(PPost);
		FreeMatrixMemory(qPost);
		PPost = KalmanOut.PPost;
		qPost = KalmanOut.qPost;
		SpEst.content[i][0] = KalmanOut.SqEst;

		for (int j = 0; j < periodicAmp.column; j++)
			periodicAmp.content[i][j] = KalmanOut.periodicAmp.content[0][j];

		//----------------Bandpass Filter-----------------
		//Matrix bandpassIn = ExtractRow(BPInput, i+1);
		Matrix bandpassIn = GetBandpassInput(measurements.content[i][0], KalmanOut.SqEst, bandpassCoefficients.NumOfFilters);
		
		BandpassOutput bandpassOutput = BandpassFilterAlgorithm(bandpassIn, A, B, EDelayed);
		EDelayed = bandpassOutput.EDelayed;	

		for (int j = 0; j < bandpassOutput.Output.row; j++)
		{
			filteredOutput.content[i][j] = bandpassOutput.Output.content[j][0];
		}

		FreeMatrixMemory(bandpassOutput.Output);
		FreeMatrixMemory(bandpassIn);

		//----------------NEO-----------------
		Matrix NEOInput = ExtractRowToColumn(filteredOutput, i + 1);
		//Matrix NEOInput = ExtractRowToColumn(singleNEOInput, i + 1);
		//printf("%.25f", NEOInput.content[0][0]);
		NEOOutput NEOOuput = NEOAlgorithm(NEOInput, D, SDelayed, samplingPeriod, prevFreq, prevAmp);
		SDelayed = NEOOuput.SDelayed;
		for (int j = 0; j < NEOOuput.Frequency.row; j++)
		{
			freqOutput.content[i][j] = NEOOuput.Frequency.content[j][0];
			ampOutput.content[i][j] = NEOOuput.Amplitude.content[j][0];
		}
		FreeMatrixMemory(NEOInput);

		//------------------------Recursive mean filter-----------------------------
		RecursiveMeanOutput recursiveMeanOutput = RecursiveMeanFilter(NEOOuput.Frequency, NEOOuput.Amplitude, ndMean, counter, prevAveFreq, prevAveAmp, bandpassCoefficients);
		for (int j = 0; j < recursiveMeanOutput.FilteredFreq.row; j++)
		{
			filteredFreqOutput.content[i][j] = recursiveMeanOutput.FilteredFreq.content[j][0];
			filteredAmpOutput.content[i][j] = recursiveMeanOutput.FilteredAmp.content[j][0];
		}
		FreeMatrixMemory(NEOOuput.Amplitude);
		FreeMatrixMemory(NEOOuput.Frequency);
		//FreeMatrixMemory(recursiveMeanOutput.FilteredFreq);
		//FreeMatrixMemory(recursiveMeanOutput.FilteredAmp);

		//------------------------Energy Ratio and Detection-----------------------------
		periodicEnergy = CalculateVibrationEnergy(KalmanOut.periodicAmp, harmonicSpindleSpeed, integrationFactor);

		Matrix filteredAmpTrans = Transpose(recursiveMeanOutput.FilteredAmp);
		Matrix filteredFreqTrans = Transpose(recursiveMeanOutput.FilteredFreq);
		chatterEnergy = CalculateVibrationEnergy(filteredAmpTrans, filteredFreqTrans, integrationFactor);
		fprintf(PeriodicEnergyPointer, "%.25f\n", periodicEnergy);
		fprintf(ChatterEnergyPointer, "%.25f\n", chatterEnergy);

		chatterDetected = ChatterDetection(periodicEnergy, chatterEnergy, energyThreshold, energyRatioLimit, prevChatterDetected, delayedPeriodicEnergy, EnergyRatioPointer);
		fprintf(ChatterDetectedPointer, "%d\n", chatterDetected);

		FreeMatrixMemory(KalmanOut.periodicAmp);
		FreeMatrixMemory(recursiveMeanOutput.FilteredFreq);
		FreeMatrixMemory(recursiveMeanOutput.FilteredAmp);
		FreeMatrixMemory(filteredAmpTrans);
		FreeMatrixMemory(filteredFreqTrans);

		printf("Iteration # %d done\n", i);
	}

	fclose(PeriodicEnergyPointer);
	fclose(ChatterEnergyPointer);
	fclose(ChatterDetectedPointer);
	fclose(EnergyRatioPointer);

	//FILE* KalmanFilePointer = fopen("..\\Kalman Output.txt", "w");
	////fprintf(OutFilePointer, "Output from the Kalman filter is:\n");
	////for (int i = 0; i < 10000; i++)
	////for (int i = 0; i < SpEst.row; i++)
	////{
	////	fprintf(KalmanFilePointer, "%.10f\n", SpEst.content[i][0]);
	////}
	////fclose(KalmanFilePointer);

	//FILE* NEOFreqPointer = fopen("..\\Frequency Output.txt", "w");
	//FILE* NEOAmpPointer = fopen("..\\Amplitude Output.txt", "w");
	//FILE* FilteredFreqPointer = fopen("..\\Filtered Freq Output.txt", "w");
	//FILE* FilteredAmpPointer = fopen("..\\Filtered Amp Output.txt", "w");
	////fprintf(OutFilePointer, "Output from the bandpass filter is:\n");
	////for (int i = 0; i < 10000; i++)
	//for (int i = 0; i < measurements.row; i++)
	//{
	//	fprintf(KalmanFilePointer, "%.10f\n", SpEst.content[i][0]);
	//	for (int j = 0; j < bandpassCoefficients.NumOfFilters; j++)
	//	{
	//		fprintf(NEOFreqPointer, "%.8f	", freqOutput.content[i][j]);
	//		fprintf(NEOAmpPointer, "%.8f	", ampOutput.content[i][j]);
	//		fprintf(FilteredFreqPointer, "%.8f	", filteredFreqOutput.content[i][j]);
	//		fprintf(FilteredAmpPointer, "%.8f	", filteredAmpOutput.content[i][j]);
	//	}
	//	fprintf(NEOFreqPointer, "\n");
	//	fprintf(NEOAmpPointer, "\n");
	//	fprintf(FilteredFreqPointer, "\n");
	//	fprintf(FilteredAmpPointer, "\n");
	//}
	//fclose(KalmanFilePointer);
	//fclose(NEOFreqPointer);
	//fclose(NEOAmpPointer);
	//fclose(FilteredFreqPointer);
	//fclose(FilteredAmpPointer);

	//FILE* NEOAmpPointer = fopen("..\\Amplitude Output.txt", "w");
	////fprintf(OutFilePointer, "Output from the bandpass filter is:\n");
	////for (int i = 0; i < 10000; i++)
	//for (int i = 0; i < filteredOutput.row; i++)
	//{
	//	for (int j = 0; j < bandpassCoefficients.NumOfFilters; j++)
	//	{
	//		fprintf(NEOAmpPointer, "%.8f	", ampOutput.content[i][j]);
	//	}
	//	fprintf(NEOAmpPointer, "\n");
	//}
	//fclose(NEOAmpPointer);

	system("pause");

	return 0;
}

BandpassCoefficients BandpassFilterCoefficients(double StartFreq, double StopFreq, double TPE, double dt, double StopBandMag, double RippleMag, double Ts)
{
	BandpassCoefficients bandpassCoefficients;
	double PassBand[2];
	double StopBand[2];
	double SamplingFreq = (double)1 / Ts;
	double m1 = ceil(StartFreq / TPE);
	double m = (ceil(StopFreq / TPE))*TPE / TPE - m1;
	double wn;
	double Wn[2];
	double* order = (double *)malloc(1 * sizeof(double));
	double* order1 = (double *)malloc(1 * sizeof(double));
	emxArray_real_T *num;
	emxArray_creal_T *den;
	emxInit_real_T1(&num, 2);
	emxInit_creal_T1(&den, 2);

	for (int i = 0; i < m; i++)
	{
		wn = m1 + i;
		PassBand[0] = (wn + dt)*TPE/(SamplingFreq/2);
		PassBand[1] = (wn + 1 - dt)*TPE/(SamplingFreq/2);

		StopBand[0] = wn*TPE / (SamplingFreq / 2);
		StopBand[1] = (wn + 1)*TPE / (SamplingFreq / 2);

		if (i == 0)
		{
			buttordBandpassOnly(PassBand, StopBand, RippleMag, StopBandMag, order, Wn);
			printf("Order is %f\n", *order);
			bandpassCoefficients.Numerator = CreateEmptyMatrix((int)m, 2 * (int)*order + 1);
			bandpassCoefficients.Denominator = CreateEmptyMatrix((int)m, 2 * (int)*order);
		}

		buttordBandpassOnly(PassBand, StopBand, RippleMag, StopBandMag, order1, Wn);
		butterBandpassOnly(*order, Wn, num, den);

		for (int j = 0; j < bandpassCoefficients.Numerator.column; j++)
			bandpassCoefficients.Numerator.content[i][j] = num->data[j];

		for (int k = 0; k < bandpassCoefficients.Denominator.column; k++)
			bandpassCoefficients.Denominator.content[i][k] = den->data[k+1].re;
	}

	bandpassCoefficients.Order = (int)*order;
	bandpassCoefficients.NumOfFilters = (int)m;
	bandpassCoefficients.StartBand = (int)m1;

	free(order);
	order = NULL;
	free(order1);
	order1 = NULL;

	return bandpassCoefficients;
}

Matrix GetBandpassInput(double measurement, double estimation, int numOfFilters)
{
	Matrix BandpassInput = CreateEmptyMatrix(numOfFilters, 1);

	for (int i = 0; i < numOfFilters; i++)
	{
		BandpassInput.content[i][0] = measurement - estimation;
	}
	return BandpassInput;
}

BandpassOutput BandpassFilterAlgorithm(Matrix input, Matrix A, Matrix B, Matrix EDelayed) // A:Denominator, B: Numerator
{
	BandpassOutput bandpassOutput;
	Matrix E = CreateEmptyMatrix(input.row, EDelayed.column + 1); // E: [m x 2n+1]
	Matrix AETrans = CreateEmptyMatrix(input.row, 1); // A*(Ed)' : [m x 1]

	// ----- For loop to calculate A*(E_delayed)'
	for (int i = 0; i < input.row; i++) // # of filters
	{
		Matrix EDelayedTrans = ExtractRowToColumn(EDelayed, i + 1); // Ed': [2n x 1]
		Matrix ASingle = ExtractRow(A, i + 1); // Extract A of #(i+1) filter: [1 x 2n]
		Matrix AETransSingle = MatrixMultiplication(ASingle, EDelayedTrans); // (A*(Ed)') of #(i+1) filter [1 x 1]
		AETrans.content[i][0] = AETransSingle.content[0][0];

		FreeMatrixMemory(EDelayedTrans);
		FreeMatrixMemory(ASingle);
		FreeMatrixMemory(AETransSingle);
	}

	Matrix ECurrent = MatrixSubtraction(input, AETrans); // [m x 1] - [m x 1] = [m x 1]

	// Update E signal
	for (int i = 0; i < input.row; i++) // # of filters
	{
		E.content[i][0] = ECurrent.content[i][0];
		for (int j = 0; j < EDelayed.column; j++)
		{
			E.content[i][j + 1] = EDelayed.content[i][j];
		}
	}

	// Calculate output signal -- y
	Matrix signalOutput = CreateEmptyMatrix(input.row, 1);

	// ---- for loop to calculate each filter's output---------
	for (int i = 0; i < input.row; i++)
	{
		Matrix ESingleTrans = ExtractRowToColumn(E, i + 1); // E': [2n+1 x 1]
		Matrix BSingle = ExtractRow(B, i + 1); // B of #(i+1) filter: [1 x 2n+1]
		Matrix signalOutputSingle = MatrixMultiplication(BSingle, ESingleTrans); // Output of #(i+1) filter: [1 x 1]
		signalOutput.content[i][0] = signalOutputSingle.content[0][0];

		FreeMatrixMemory(ESingleTrans);
		FreeMatrixMemory(BSingle);
		FreeMatrixMemory(signalOutputSingle);
	}

	// Update delayed E
	for (int i = 0; i < input.row; i++)
	{
		for (int j = 0; j < EDelayed.column - 1; j++)
		{
			EDelayed.content[i][EDelayed.column - 1 - j] = EDelayed.content[i][EDelayed.column - 2 - j];
		}
		EDelayed.content[i][0] = ECurrent.content[i][0];
	}

	bandpassOutput.Output = signalOutput;
	bandpassOutput.EDelayed = EDelayed;

	FreeMatrixMemory(E);
	FreeMatrixMemory(AETrans);
	FreeMatrixMemory(ECurrent);

	return bandpassOutput;
}

Matrix CalculateLag(double samplingPeriod, double TPE, BandpassCoefficients bandpassFilterCoefficient)
{
	Matrix D = CreateEmptyMatrix(1, bandpassFilterCoefficient.NumOfFilters); // Initialize lag matrix 

	for (int i = 0; i < bandpassFilterCoefficient.NumOfFilters; i++)
	{
		D.content[0][i] = round((1 / samplingPeriod) / (4 * (bandpassFilterCoefficient.StartBand + i + 1)*TPE + TPE / 2) + 0.5);
	}

	return D;
}

double Max(Matrix matrix)
{
	double max = matrix.content[0][0];

	for (int i = 0; i < matrix.row; i++)
	{
		for (int j = 0; j < matrix.column; j++)
		{
			if (matrix.content[i][j] > max)
				max = matrix.content[i][j];
		}
	}

	return max;
}

NEOOutput NEOAlgorithm(Matrix input, Matrix D, Matrix s_delay, double T, Matrix prevFreq, Matrix prevAmp)
{
	Matrix phi_y_kD = CreateEmptyMatrix(input.row, 1);
	Matrix phi_y_k2D = CreateEmptyMatrix(input.row, 1);
	Matrix phi_s_k2D = CreateEmptyMatrix(input.row, 1);

	Matrix s_k = CreateEmptyMatrix(input.row, 1);
	Matrix s_kD = CreateEmptyMatrix(input.row, 1);
	Matrix s_k2D = CreateEmptyMatrix(input.row, 1);
	Matrix s_k3D = CreateEmptyMatrix(input.row, 1);
	Matrix s_k4D = CreateEmptyMatrix(input.row, 1);

	Matrix freq = CreateEmptyMatrix(input.row, 1);
	Matrix amp = CreateEmptyMatrix(input.row, 1);

	// Update input from bandpass filter
	for (int i = 0; i < input.row; i++)
	{
		s_k.content[i][0] = input.content[i][0];
		s_kD.content[i][0] = s_delay.content[i][s_delay.column - (int)D.content[0][i]];
		s_k2D.content[i][0] = s_delay.content[i][s_delay.column - 2 * (int)D.content[0][i]];
		s_k3D.content[i][0] = s_delay.content[i][s_delay.column - 3 * (int)D.content[0][i]];
		s_k4D.content[i][0] = s_delay.content[i][s_delay.column - 4 * (int)D.content[0][i]];
	}

	// Calculate y_k, y_k-D, y_k-2D and y_k-3D // Must initialize here to avoid memory issue
	Matrix y_k = MatrixSubtraction(s_k, s_kD);
	Matrix y_kD = MatrixSubtraction(s_kD, s_k2D);
	Matrix y_k2D = MatrixSubtraction(s_k2D, s_k3D);
	Matrix y_k3D = MatrixSubtraction(s_k3D, s_k4D);

	// Calculate Phi[y_(k-D)], Phi[y_(k-2D)], Phi[s_(k-2D)]
	for (int i = 0; i < input.row; i++)
	{
		//phi_y_kD.content[i][0] = 1 / pow(T, 2)*(pow(y_kD.content[i][0], 2) - y_k2D.content[i][0] * y_k.content[i][0]);
		//phi_y_k2D.content[i][0] = 1 / pow(T, 2)*(pow(y_k2D.content[i][0], 2) - y_k3D.content[i][0] * y_kD.content[i][0]);
		//phi_s_k2D.content[i][0] = 1 / pow(T, 2)*(pow(s_k2D.content[i][0], 2) - s_k3D.content[i][0] * s_kD.content[i][0]);
		phi_y_kD.content[i][0] = pow(y_kD.content[i][0], 2) - y_k2D.content[i][0] * y_k.content[i][0];
		phi_y_k2D.content[i][0] = pow(y_k2D.content[i][0], 2) - y_k3D.content[i][0] * y_kD.content[i][0];
		phi_s_k2D.content[i][0] = pow(s_k2D.content[i][0], 2) - s_k3D.content[i][0] * s_kD.content[i][0];
	}

	// Calculate frequency and amplitude
	for (int i = 0; i < input.row; i++)
	{
		double phi_xd = phi_y_k2D.content[i][0] + phi_y_kD.content[i][0];
		if (phi_s_k2D.content[i][0] > 0 && phi_xd > 0 && phi_xd < 8 * phi_s_k2D.content[i][0])
		{
			freq.content[i][0] = acos(1 - phi_xd / (4 * phi_s_k2D.content[i][0])) / 2.0 / PI / T / D.content[0][i]; // in Hz
			amp.content[i][0] = pow(phi_s_k2D.content[i][0]/(1 - pow(1- phi_xd/(4*phi_s_k2D.content[i][0]), 2)), 0.5);
		}
		else
		{
			freq.content[i][0] = prevFreq.content[i][0]; // in Hz
			amp.content[i][0] = prevAmp.content[i][0];
		}
	}

	// Update delayed input from bandpass filter
	for (int i = 0; i < input.row; i++)
	{
		for (int j = 0; j < s_delay.column - 1; j++)
		{
			s_delay.content[i][j] = s_delay.content[i][j+1];
		}
		s_delay.content[i][s_delay.column - 1] = input.content[i][0];
	}

	// Update previous frequency and amplitude
	for (int i = 0; i < input.row; i++)
	{
		prevFreq.content[i][0] = freq.content[i][0];
		prevAmp.content[i][0] = amp.content[i][0];
	}

	NEOOutput output = {freq, amp, s_delay};
	//NEOOutput output = { freq, phi_s_k2D, s_delay };
	FreeMatrixMemory(y_k);
	FreeMatrixMemory(y_kD);
	FreeMatrixMemory(y_k2D);
	FreeMatrixMemory(y_k3D);
	FreeMatrixMemory(phi_y_kD);
	FreeMatrixMemory(phi_y_k2D);
	FreeMatrixMemory(phi_s_k2D);

	FreeMatrixMemory(s_k);
	FreeMatrixMemory(s_kD);
	FreeMatrixMemory(s_k2D);
	FreeMatrixMemory(s_k3D);
	FreeMatrixMemory(s_k4D);

	return output;
}

RecursiveMeanOutput RecursiveMeanFilter(Matrix freqIn, Matrix ampIn, double ndMean, double* counter, Matrix prevAveFreq, Matrix prevAveAmp, BandpassCoefficients bandpassCoefficients)
{
	Matrix filteredFreq = CreateEmptyMatrix(bandpassCoefficients.NumOfFilters, 1);
	Matrix filteredAmp = CreateEmptyMatrix(bandpassCoefficients.NumOfFilters, 1);

	if (*counter < ndMean)
	{
		if (*counter != 0)
		{
			for (int i = 0; i < bandpassCoefficients.NumOfFilters; i++)
			{
				filteredFreq.content[i][0] = (*counter - 1) / (*counter)*prevAveFreq.content[i][0] + 1 / (*counter)*freqIn.content[i][0];
				filteredAmp.content[i][0] = (*counter - 1) / (*counter)*prevAveAmp.content[i][0] + 1 / (*counter)*ampIn.content[i][0];
			}
		}
	}

	else
	{
		for (int i = 0; i < bandpassCoefficients.NumOfFilters; i++)
		{
			filteredFreq.content[i][0] = ((ndMean - 1) / ndMean*prevAveFreq.content[i][0] + 1 / ndMean*freqIn.content[i][0]) * 2* PI; // In [rad/s]
			filteredAmp.content[i][0] = (ndMean - 1) / ndMean*prevAveAmp.content[i][0] + 1 / ndMean*ampIn.content[i][0];
		}
	}

	// Update previous average frequency and amplitude
	for (int i = 0; i < bandpassCoefficients.NumOfFilters; i++)
	{
		prevAveFreq.content[i][0] = filteredFreq.content[i][0] / 2 / PI; // in Hz
		prevAveAmp.content[i][0] = filteredAmp.content[i][0];
	}

	//Update counter;
	(*counter)++;

	RecursiveMeanOutput recursiveMeanOutput = { filteredFreq, filteredAmp };
	return recursiveMeanOutput;
}

Matrix CreateHarmonicSpindleSpeed(int N, double spindleSpeed) // spindleSpeed in [rad/s]
{
	Matrix HarmonicSpindleSpeed = CreateEmptyMatrix(1, N);

	for (int i = 0; i < HarmonicSpindleSpeed.column; i++)
	{
		HarmonicSpindleSpeed.content[0][i] = (i + 1)*spindleSpeed; // [rad/s]
	}

	return HarmonicSpindleSpeed;
}

double CalculateVibrationEnergy(Matrix amplitude, Matrix frequency, double integrationFactor) //amplitude: [1xN] frequency: [1xN]
{
	double sum = 0;

	for (int i = 0; i < amplitude.column; i++)
	{
		if (frequency.content[0][i] != 0)
			sum = sum + pow(amplitude.content[0][i], 2)*pow(frequency.content[0][i], 2 * integrationFactor);
	}

	return sum;
}

int ChatterDetection(double periodicEnergy, double chatterEnergy, double energyThreshold, double energyRatioLimit, int* prevChatterDetected, Matrix delayedPeriodicEnergy, FILE* outputFilePtr)
{
	int chatterDetected;

	double energyRatio = chatterEnergy / (chatterEnergy + delayedPeriodicEnergy.content[0][0]);
	fprintf(outputFilePtr, "%.25f\n", energyRatio);

	if (*prevChatterDetected == 0 && energyRatio > 0.5) // entry
		chatterDetected = 0;
	else if (*prevChatterDetected == 1 && energyRatio > 0.1) // chatter continues
		chatterDetected = 1;
	else if (energyRatio > energyRatioLimit && delayedPeriodicEnergy.content[0][0] > energyThreshold) // chatter just detected
		chatterDetected = 1;
	else if (energyRatio < 0.1) //no chatter
		chatterDetected = 0;
	else
		chatterDetected = 0;

	// Update if chatter is detected by previous data
	*prevChatterDetected = chatterDetected;

	// Update delayed periodic energy
	for (int i = 0; i < delayedPeriodicEnergy.column - 1; i++)
	{
		delayedPeriodicEnergy.content[0][i] = delayedPeriodicEnergy.content[0][i+1];
	}
	delayedPeriodicEnergy.content[0][delayedPeriodicEnergy.column - 1] = periodicEnergy;

	return chatterDetected;
}