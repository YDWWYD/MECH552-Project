#include "Kalman.h"
Matrix CreateInitial_q(int N)
{
	Matrix initial_q = CreateEmptyMatrix(2 * N, 1);

	//printf("------q0 Matrix-----\n");
	//PrintMatrix(initial_q);

	return initial_q;
}