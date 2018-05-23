#include "MatrixOperations.h"

void FreeMatrixMemory(Matrix matrix)
{
	for (int i = 0; i < matrix.row; i++)
	{
		free(matrix.content[i]);
		matrix.content[i] = NULL;
	}

	free(matrix.content);
	matrix.content = NULL;
}
