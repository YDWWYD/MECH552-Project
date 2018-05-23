#include "Kalman.h"

Matrix CreateHMatrix(int N)
{
	Matrix H = CreateEmptyMatrix(1, 2 * N);

	for (int i = 0; i < N; i++)
		H.content[0][2 * i] = 1;

	//printf("------ H Matrix------\n");
	//PrintMatrix(H);

	return H;
}