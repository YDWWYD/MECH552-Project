#include "MatrixOperations.h"

double Determinant(Matrix matrix)
{
	double det = 0, coefficient = 1;

	if (matrix.row != matrix.column)
	{
		printf("Invalid Dimension. Determinant does not exist!\n");
		return;
	}

	if (matrix.row == 1)
		det = matrix.content[0][0];
	else
	{
		Matrix minor = CreateEmptyMatrix(matrix.row - 1, matrix.column - 1);
		for (int c = 0; c < matrix.row; c++) // Loop through the 1st row
		{
			int m = 0;
			int n = 0;

			for (int i = 0; i < matrix.row; i++)
			{
				for (int j = 0; j < matrix.row; j++)
				{
					if (i != 0 && j != c) // Exclude the current row and column of elements
					{
						minor.content[m][n] = matrix.content[i][j];
						if (n < matrix.column - 2)
							n++;
						else
						{
							n = 0;
							m++;
						}
					}
				}
			}

			det = det + coefficient*(matrix.content[0][c] * Determinant(minor));
			coefficient = -1 * coefficient;
		}
	}

	return det;
}