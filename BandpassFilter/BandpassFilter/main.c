//#include "rt_nonfinite.h"
#include "buttordBandpassOnly.h"
//#include "sortButtord.h"
#include "stdio.h"
//#include "rt_nonfinite.h"
#include "butterBandpassOnly.h"
//#include "butterBandpassOnly_emxutil.h"
//#include "exp.h"
//#include "isequal.h"
//#include "sortButter.h"
//#include "relop.h"
//#include "mypoly.h"
//#include "mybilinear.h"
//#include "mylp2bs.h"
//#include "myzp2ss.h"
//#include "butterBandpassOnly_rtwutil.h"
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
}BandpassCoefficients;

typedef struct {
	Matrix Output;
	Matrix EDelayed;
} BandpassOutput;

BandpassCoefficients BandpassFilterCofficients(double StartFreq, double StopFreq, double TPE, double dt, double StopBandMag, double RippleMag, double Ts);
Matrix GetBandpassInput(double measurement, double estimation, int numOfFilters);
//Matrix BandpassFilterAlgorithm(Matrix input, BandpassCoefficients bandpassCoefficients);
BandpassOutput BandpassFilterAlgorithm(Matrix input, Matrix A, Matrix B, Matrix EDelayed);

int main(void)
{
	/*double wp[] = { 0.3781, 0.4031 };
	double ws[] = { 0.3750, 0.4063 };
	double rp = 2, rs = 5;
	double n = 1;
	double* order = (double *)malloc(1 * sizeof(double));
	double wn[2];
	buttordBandpassOnly(wp, ws, rp, rs, order, wn);

	emxArray_real_T *num;
	emxArray_creal_T *den;
	emxInit_real_T1(&num, 2);
	emxInit_creal_T1(&den, 2);

	butterBandpassOnly(*order, wn, num, den);

	for (int i = 0; i < num->size[1]; i++)
	{
		printf("num[%d] is %.10f\n", i, num->data[i]);
	}
		
	for (int i = 0; i < den->size[1]; i++)
	{
		printf("Re(den[%d]) is %.8f		Im(den[%d]) is %f\n", i, den->data[i].re, i, den->data[i].im);
	}

	free(order);
	order = NULL;*/

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
	double samplingPeriod = 3.90625E-5;
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

	BandpassCoefficients bandpassCoefficients = BandpassFilterCofficients(startFreq, stopFreq, TPE, dt, stopBandMag, rippleMag, Ts);
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
	//Matrix filteredOutput = CreateEmptyMatrix(measurements.row, 1);
	Matrix filteredOutput = CreateEmptyMatrix(measurements.row, bandpassCoefficients.NumOfFilters);

	for (int i = 0; i < measurements.row; i++)
	//for (int i = 0; i < 10000; i++)
	{
		//----------------Kalman Filter-----------------
		KalmanOutput KalmanOut = KalmanAlgorithm(measurements.content[i][0], Phi, PhiTranspose, H, HTranspose, identity2N, R, Q, PPost, qPost);
		FreeMatrixMemory(PPost);
		FreeMatrixMemory(qPost);
		PPost = KalmanOut.PPost;
		qPost = KalmanOut.qPost;
		SpEst.content[i][0] = KalmanOut.SqEst;

		//----------------Bandpass Filter-----------------
		Matrix bandpassIn = GetBandpassInput(measurements.content[i][0], KalmanOut.SqEst, bandpassCoefficients.NumOfFilters);
		//Matrix singleFilterIn = ExtractRow(bandpassIn, 1);

		BandpassOutput bandpassOutput = BandpassFilterAlgorithm(bandpassIn, A, B, EDelayed);
		EDelayed = bandpassOutput.EDelayed;	

		for (int j = 0; j < bandpassOutput.Output.row; j++)
		{
			filteredOutput.content[i][j] = bandpassOutput.Output.content[j][0];
		}

		FreeMatrixMemory(bandpassOutput.Output);
		FreeMatrixMemory(bandpassIn);

		printf("Iteration # %d done\n", i);
	}

	FILE* OutFilePointer = fopen("..\\Bandpass Filter Output.txt", "w");
	//fprintf(OutFilePointer, "Output from the bandpass filter is:\n");
	//for (int i = 0; i < 10000; i++)
	for (int i = 0; i < filteredOutput.row; i++)
	{
		for (int j = 0; j < bandpassCoefficients.NumOfFilters; j++)
		{
			fprintf(OutFilePointer, "%.8f	", filteredOutput.content[i][j]);
		}
		fprintf(OutFilePointer, "\n");
	}

	fclose(OutFilePointer);

	system("pause");

	return 0;
}

BandpassCoefficients BandpassFilterCofficients(double StartFreq, double StopFreq, double TPE, double dt, double StopBandMag, double RippleMag, double Ts)
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

//Matrix BandpassFilterAlgorithm(Matrix input, BandpassCoefficients bandpassCoefficients)
//{
//	//Matrix filterOutput = CreateEmptyMatrix(input.row, 1); //Out: [m x 1]
//	//Matrix E = CreateEmptyMatrix(input.row, bandpassCoefficients.Order * 2 + 1); // E: [m x 2n+1]
//	//Matrix EDelayed = CreateEmptyMatrix(input.row, bandpassCoefficients.Order * 2); // Delayed E: [m x 2n]
//	//Matrix ECurrent = CreateEmptyMatrix(input.row, 1); //Ek: [m x 1]
//
//	Matrix filterOutput = CreateEmptyMatrix(1, 1); //y: [1 x 1]
//	Matrix E = CreateEmptyMatrix(1, bandpassCoefficients.Order * 2 + 1); // E: [1 x 2n+1]
//	//Matrix ETrans = Transpose(E); // E': [2n+1 x 1]
//	Matrix EDelayed = CreateEmptyMatrix(1, bandpassCoefficients.Order * 2); // Delayed E: [1 x 2n]
//	//Matrix EDelayedTrans = Transpose(EDelayed); // (DelayedE)': [2n x 1]
//	Matrix ECurrent = CreateEmptyMatrix(1, 1); //Ek: [1 x 1]
//	Matrix A = ExtractRow(bandpassCoefficients.Denominator, 1); // A: [1 x 2n]
//	Matrix B = ExtractRow(bandpassCoefficients.Numerator, 1); // B: [1 x 2n+1]
//
//
//	Matrix EDelayedTrans = Transpose(EDelayed); // (DelayedE)': [2n x 1]
//	Matrix AETrans = MatrixMultiplication(A, EDelayedTrans);
//	ECurrent = MatrixSubtraction(input, AETrans);
//
//	// Update E signal
//	E.content[0][0] = ECurrent.content[0][0];
//	for (int i = 0; i < EDelayed.column; i++)
//	{
//		E.content[0][i + 1] = EDelayed.content[0][i];
//	}
//
//	// Calculate output signal - y
//	Matrix ETrans = Transpose(E); // E': [2n+1 x 1]
//	filterOutput = MatrixMultiplication(B, ETrans);
//
//	// Update delayed E
//	for (int i = 0; i < EDelayed.column - 1; i++)
//	{
//		EDelayed.content[0][EDelayed.column - 1 - i] = EDelayed.content[0][EDelayed.column - 2 - i];
//	}
//	EDelayed.content[0][0] = ECurrent.content[0][0];
//
//	return filterOutput;
//}	

//BandpassOutput BandpassFilterAlgorithm(Matrix input, Matrix A, Matrix B, Matrix EDelayed)
//{
//	BandpassOutput bandpassOutput;
//	Matrix E = CreateEmptyMatrix(1, EDelayed.column + 1); // E: [1 x 2n+1]
//
//	Matrix EDelayedTrans = Transpose(EDelayed); // (DelayedE)': [2n x 1]
//	Matrix AETrans = MatrixMultiplication(A, EDelayedTrans);
//	Matrix ECurrent = MatrixSubtraction(input, AETrans);
//
//	// Update E signal
//	E.content[0][0] = ECurrent.content[0][0];
//	for (int i = 0; i < EDelayed.column; i++)
//	{
//		E.content[0][i + 1] = EDelayed.content[0][i];
//	}
//
//	// Calculate output signal - y
//	Matrix ETrans = Transpose(E); // E': [2n+1 x 1]
//	Matrix signalOutput = MatrixMultiplication(B, ETrans);
//
//	// Update delayed E
//	for (int i = 0; i < EDelayed.column - 1; i++)
//	{
//		EDelayed.content[0][EDelayed.column - 1 - i] = EDelayed.content[0][EDelayed.column - 2 - i];
//	}
//	EDelayed.content[0][0] = ECurrent.content[0][0];
//
//	bandpassOutput.Output = signalOutput;
//	bandpassOutput.EDelayed = EDelayed;
//
//	FreeMatrixMemory(E);
//	FreeMatrixMemory(EDelayedTrans);
//	FreeMatrixMemory(AETrans); 
//	FreeMatrixMemory(ECurrent);
//	FreeMatrixMemory(ETrans);
//
//	return bandpassOutput;
//}

BandpassOutput BandpassFilterAlgorithm(Matrix input, Matrix A, Matrix B, Matrix EDelayed)
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