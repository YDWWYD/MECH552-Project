#include "Kalman.h"

Matrix CreateQMatrix(int N, double lamda, Matrix R)
{
	Matrix Q = CreateEmptyMatrix(2 * N, 2 * N);

	for (int i = 0; i < 2 * N; i++)
		Q.content[i][i] = lamda * R.content[0][0];

	//printf("------QMatrix------\n");
	//for (int i = 0; i < 2 * N; i++)
	//{
	//	for (int j = 0; j < 2 * N; j++)
	//		printf("%.15f  ", Q.content[i][j]);
	//	printf("\n");
	//}

	return Q;
}