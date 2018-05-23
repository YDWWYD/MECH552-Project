#include "MatrixOperations.h"

Matrix ExtractRowToColumn(Matrix inputMatrix, int row)
{
	if (row > inputMatrix.row)
	{
		printf("Row exceeds the max row.\n");
		return;
	}

	Matrix columnMatrix = CreateEmptyMatrix(inputMatrix.column, 1);

	for (int i = 0; i < columnMatrix.row; i++)
	{
		columnMatrix.content[i][0] = inputMatrix.content[row - 1][i];
	}
	return columnMatrix;
}