#include "MatrixOperations.h"

Matrix IdentityMatrix(int N)
{
	Matrix identity = CreateEmptyMatrix(N, N);

	for (int i = 0; i < N; i++)
		identity.content[i][i] = 1;

	return identity;
}