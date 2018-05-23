#include "Kalman.h"
Matrix CreateInitialP(int N)
{
	Matrix initialP = CreateEmptyMatrix(2 * N, 2 * N);

	for (int i = 0; i < initialP.row; i++)
		initialP.content[i][i] = 10;

	//printf("------ P0 Matrix------\n");
	//PrintMatrix(initialP);

	return initialP;
}