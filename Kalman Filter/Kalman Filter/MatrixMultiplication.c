#include "MatrixOperations.h"

Matrix MatrixMultiplication(Matrix matrix1, Matrix matrix2)
{
	if (matrix1.column != matrix2.row)
	{
		printf("Dimension Invalid!\n");
		return;
	}

	Matrix result = CreateEmptyMatrix(matrix1.row, matrix2.column);

	for (int i = 0; i < result.row; i++)
	{
		for (int j = 0; j < result.column; j++)
		{
			//result[i, j] = 0; //Default values are 0
			for (int k = 0; k < matrix1.column; k++)
				result.content[i][j] = result.content[i][j] + matrix1.content[i][k] * matrix2.content[k][j];
		}
	}

	//printf("------ Multiplication Result------\n");
	//PrintMatrix(result);

	return result;
}