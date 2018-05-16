#include "MatrixOperations.h"

Matrix MatrixSubtraction(Matrix matrix1, Matrix matrix2)
{
	if (matrix1.row != matrix2.row || matrix1.column != matrix2.column)
	{
		printf("Dimension Invalid!\n");
		return;
	}

	Matrix result = CreateEmptyMatrix(matrix1.row, matrix1.column);

	for (int i = 0; i < matrix1.row; i++)
	{
		for (int j = 0; j < matrix1.column; j++)
			result.content[i][j] = matrix1.content[i][j] - matrix2.content[i][j];
	}

	//printf("------ Subtraction Result------\n");
	//PrintMatrix(result);
	return result;
}