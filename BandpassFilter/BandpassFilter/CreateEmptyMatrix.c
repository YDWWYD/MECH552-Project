#include "MatrixOperations.h"

Matrix CreateEmptyMatrix(int row, int column)
{
	Matrix matrix;

	matrix.row = row;
	matrix.column = column;
	matrix.content = (double**)calloc(matrix.row, sizeof(double*));
	for (int i = 0; i < matrix.row; i++)
		*(matrix.content + i) = (double *)calloc(matrix.column, sizeof(double));

	return matrix;
}