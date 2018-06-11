#include "Matrix.h"
#include "MatrixOperations.h"
#include "Kalman.h"

int main(void)
{
	int N = 20;
	double spindleSpeed = 200 * 2 * PI;
	double samplingPeriod = (double)1 / 25600; //Ts calculated by DAQ
	double lamda = 1E-6;
	char* filePath = "..\\Acc2.txt";
	//char* filePath = "..\\KalmanAcc2.txt";

	Matrix Phi = CreatePhiMatrix(N, spindleSpeed, samplingPeriod);
	Matrix H = CreateHMatrix(N);

	Matrix measurements = ReadMeasurements(filePath);
	printf("%.15f\n", measurements.content[0][0]);
	Matrix R = CalculateVariance(measurements, 100, 700);
	printf("R = %.25f\n", R.content[0][0]);

	Matrix Q = CreateQMatrix(N, lamda, R);
	printf("Q = %.28f\n", Q.content[0][0]);
	Matrix q0 = CreateInitial_q(N);
	Matrix P0 = CreateInitialP(N);

	Matrix PPost = P0;
	Matrix qPost = q0;

	Matrix PhiTranspose = Transpose(Phi);
	Matrix identity2N = IdentityMatrix(P0.row); // Identity [2N * 2N]
	Matrix HTranspose = Transpose(H);
	Matrix SpEst = CreateEmptyMatrix(measurements.row, 1);
	Matrix periodicAmp = CreateEmptyMatrix(measurements.row, q0.row/2);
	//Matrix SpEst = KalmanFilter(filePath, N, spindleSpeed, samplingPeriod, lamda);

	for (int i = 0; i < measurements.row; i++)
	//for (int i = 0; i < 1000; i++)
	{	
		//printf("Input %d is %.10f\n", i+1, measurements.content[i][0]);
		KalmanOutput KalmanOut = KalmanAlgorithm(measurements.content[i][0], Phi, PhiTranspose, H, HTranspose, identity2N, R, Q, PPost, qPost);

		FreeMatrixMemory(PPost);
		FreeMatrixMemory(qPost);

		PPost = KalmanOut.PPost;
		qPost = KalmanOut.qPost;

		SpEst.content[i][0] = KalmanOut.SqEst;

		for (int j = 0; j < periodicAmp.column; j++)
		{
			periodicAmp.content[i][j] = KalmanOut.periodicAmp.content[0][j];
		}
		FreeMatrixMemory(KalmanOut.periodicAmp);

		//printf("Iteration number %d...\n", i);
	}

	FILE* OutFilePointer = fopen("..\\Kalman Output.txt", "w");
	FILE* PeriodicAmplitudePointer = fopen("..\\Kalman Periodic Amplitude.txt", "w");
	//fprintf(OutFilePointer, "Output from the Kalman filter is:\n");
	//for (int i = 0; i < 1000; i++)
	for (int i = 0; i < SpEst.row; i++)
	{
		fprintf(OutFilePointer, "%.10f\n", SpEst.content[i][0]);
		for (int j = 0; j < N; j++)
		{
			fprintf(PeriodicAmplitudePointer, "%.8f	", periodicAmp.content[i][j]);
		}
		fprintf(PeriodicAmplitudePointer, "\n");
	}

	fclose(OutFilePointer);
	fclose(PeriodicAmplitudePointer);

	system("pause");
	return 0;
}
