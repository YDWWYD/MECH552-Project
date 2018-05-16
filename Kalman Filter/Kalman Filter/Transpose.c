#include "MatrixOperations.h"

Matrix Transpose(Matrix matrix)
{
	Matrix matrixTransposed = CreateEmptyMatrix(matrix.column, matrix.row);

	for (int i = 0; i < matrixTransposed.row; i++)
	{
		for (int j = 0; j < matrixTransposed.column; j++)
			matrixTransposed.content[i][j] = matrix.content[j][i];
	}

	//printf("------ Transposed Matrix------\n");
	//PrintMatrix(matrixTransposed);

	return matrixTransposed;
}