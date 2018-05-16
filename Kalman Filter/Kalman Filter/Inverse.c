#include "MatrixOperations.h"

Matrix Inverse(Matrix matrix)
{
	Matrix inverse = CreateEmptyMatrix(matrix.row, matrix.column);
	Matrix cofactor = CreateEmptyMatrix(matrix.row, matrix.column);
	Matrix minor = CreateEmptyMatrix(matrix.row - 1, matrix.column - 1);

	double det = Determinant(matrix);

	if (det == 0)
	{
		printf("Singular matrix, no inverse!\n");
		return;
	}

	if (matrix.row == 1) // matrix is a number -> inverse = 1/det
	{
		inverse.content[0][0] = 1 / det;
	}

	else
	{
		for (int p = 0; p < matrix.row; p++)
		{
			for (int q = 0; q < matrix.column; q++)
			{
				int m = 0;
				int n = 0;

				for (int i = 0; i < matrix.row; i++)
				{
					for (int j = 0; j < matrix.column; j++)
					{
						if (i != p && j != q)
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
				cofactor.content[p][q] = pow(-1, p + q) * Determinant(minor);
			}
		}

		Matrix cofactorTransposed = Transpose(cofactor);

		for (int i = 0; i < matrix.row; i++)
		{
			for (int j = 0; j < matrix.column; j++)
				inverse.content[i][j] = cofactorTransposed.content[i][j] / det;
		}
	}

	//printf("------------Inversed Maitrx:----------\n");
	//PrintMatrix(inverse);

	return inverse;
}