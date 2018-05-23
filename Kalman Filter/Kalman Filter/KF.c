#include "Matrix.h"
#include "MatrixOperations.h"
#include "Kalman.h"

int main(void)
{
	int N = 20;
	double spindleSpeed = 200 * 2 * PI;
	double samplingPeriod = 3.90625E-5;
	double lamda = 1E-6;
	char* filePath = "..\\Acc2.txt";

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
	//Matrix SpEst = KalmanFilter(filePath, N, spindleSpeed, samplingPeriod, lamda);

	//for (int i = 0; i < measurements.row; i++)
	for (int i = 0; i < 100; i++)
	{		
		KalmanOutput KalmanOut = KalmanAlgorithm(measurements.content[i][0], Phi, PhiTranspose, H, HTranspose, identity2N, R, Q, PPost, qPost);

		FreeMatrixMemory(PPost);
		FreeMatrixMemory(qPost);

		PPost = KalmanOut.PPost;
		qPost = KalmanOut.qPost;

		SpEst.content[i][0] = KalmanOut.SqEst;
		printf("Iteration number %d...\n", i);
	}

	FILE* OutFilePointer = fopen("..\\Kalman Output.txt", "w");
	fprintf(OutFilePointer, "Output from the Kalman filter is:\n");
	for (int i = 0; i < 100; i++)
	//for (int i = 0; i < SpEst.row; i++)
	{
		fprintf(OutFilePointer, "%.6f\n", SpEst.content[i][0]);
	}

	fclose(OutFilePointer);

	system("pause");
	return 0;
}
