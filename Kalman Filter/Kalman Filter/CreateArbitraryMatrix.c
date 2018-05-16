#include "MatrixOperations.h"

Matrix CreateArbitraryMatrix(int row, int column)
{
	Matrix matrix = CreateEmptyMatrix(row, column);

	printf("Please enter %d rows and %d columns of numbers: ", row, column);

	for (int i = 0; i < matrix.row; i++)
	{
		for (int j = 0; j < matrix.column; j++)
			scanf("%lf", &(matrix.content[i][j]));
	}

	printf("------------ Matrix is£º---------------\n");
	PrintMatrix(matrix);

	return matrix;
}