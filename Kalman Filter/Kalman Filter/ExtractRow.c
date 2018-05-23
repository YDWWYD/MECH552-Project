#include "MatrixOperations.h"

Matrix ExtractRow(Matrix matrix, int row)
{
	if (row > matrix.row)
	{
		printf("Row exceeds the max row.\n");
		return;
	}

	Matrix rowMatrix = CreateEmptyMatrix(1, matrix.column);

	for (int i = 0; i < rowMatrix.column; i++)
	{
		rowMatrix.content[0][i] = matrix.content[row-1][i];
	}
	return rowMatrix;
}