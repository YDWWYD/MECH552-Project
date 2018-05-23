#include "Kalman.h"

Matrix CreatePhiMatrix(int N, double spindleSpeed, double samplingPeriod)
{
	Matrix phiMatrix = CreateEmptyMatrix(2 * N, 2 * N);

	for (int n = 1; n <= N; n++)
	{
		phiMatrix.content[2 * (n - 1)][2 * (n - 1)] = cos(n*spindleSpeed*samplingPeriod);
		phiMatrix.content[2 * (n - 1)][2 * (n - 1) + 1] = -sin(n*spindleSpeed*samplingPeriod);
		phiMatrix.content[2 * (n - 1) + 1][2 * (n - 1)] = sin(n*spindleSpeed*samplingPeriod);
		phiMatrix.content[2 * (n - 1) + 1][2 * (n - 1) + 1] = cos(n*spindleSpeed*samplingPeriod);
	}

	//printf("------ Phi Matrix------\n");
	//PrintMatrix(phiMatrix);

	return phiMatrix;
}