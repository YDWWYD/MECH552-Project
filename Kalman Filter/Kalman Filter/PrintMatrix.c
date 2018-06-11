#include "MatrixOperations.h"

void PrintMatrix(Matrix matrix)
{
	for (int i = 0; i < matrix.row; i++)
	{
		for (int j = 0; j < matrix.column; j++)
			printf("%.15f  ", matrix.content[i][j]);
		printf("\n");
	}
}
