#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define PI acos(-1.0)

typedef struct {
	int row;
	int column;
	double** content;
} Matrix;

Matrix CreateArbitraryMatrix(int row, int column);
Matrix CreateMatrix(int row, int column);
Matrix CreateEmptyMatrix(int row, int column);
Matrix CreatePhiMatrix(int N, double spindleSpeed, double samplingPeriod);
Matrix CreateHMatrix(int N);
Matrix CreateQMatrix(int N, double lamda, double R);
double* ReadMeasurements(char* filePath);
double CalculateVariance(double* measurement, int rangeBottom, int rangeTop);
Matrix CreateInitialQ(int N);
Matrix CreateInitialP(int N);
Matrix MatrixTranspose(Matrix matrix);
Matrix MatirxMultiplication(Matrix matrix1, Matrix matrix2);
Matrix MatirxAddition(Matrix matrix1, Matrix matrix2);
Matrix MatirxSubtraction(Matrix matrix1, Matrix matrix2);
double Determinant(Matrix matrix);
Matrix InverseMatrix(Matrix matrix);
void PrintMatrix(Matrix matrix);

int main(void)
{
	int N = 10;
	double spindleSpeed = 200 * 2 * PI;
	double samplingPeriod = 3.90625E-5;
	double lamda = 1E-6;

	Matrix matrix = CreateMatrix(2 * N, 1 * N );
	Matrix phiMatrix = CreatePhiMatrix(N, spindleSpeed, samplingPeriod);
	Matrix HMatrix = CreateHMatrix(N);

	char* filePath = "G:\\MECH552\\yudi\\data\\Acc2.txt";

	double* meausrements = ReadMeasurements(filePath);
	double R = CalculateVariance(meausrements, 100, 700);

	Matrix QMatrix = CreateQMatrix(N, lamda, R);
	//Matrix initialQ = CreateInitialQ(N);
	//Matrix initialP = CreateInitialP(N);

	Matrix PhiTransposed = MatrixTranspose(matrix);

	Matrix matrix1 = CreateArbitraryMatrix(3, 3);
	Matrix matrix2 = CreateArbitraryMatrix(3, 3);
//	Matrix result = MatirxMultiplication(matrix1, matrix2);
	Matrix addResult = MatirxAddition(matrix1, matrix2);
	Matrix minusResult = MatirxSubtraction(matrix1, matrix2);

	//Matrix matrix3 = CreateArbitraryMatrix(4, 4);
	//double det = Determinant(matrix3);
	//printf("----Determinant is: %.3f.\n", det);
	//Matrix inverse = InverseMatrix(matrix3);

	system("pause");
	return 0;
}

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

void PrintMatrix(Matrix matrix)
{
	for (int i = 0; i < matrix.row; i++)
	{
		for (int j = 0; j < matrix.column; j++)
			printf("%.3f  ", matrix.content[i][j]);
		printf("\n");
	}
}

Matrix CreateArbitraryMatrix(int row, int column)
{
	Matrix matrix;

	matrix.row = row;
	matrix.column = column;
	matrix.content = (double**)calloc(matrix.row, sizeof(double*));
	for (int i = 0; i < matrix.row; i++)
		*(matrix.content + i) = (double *)calloc(matrix.column, sizeof(double));

	printf("Please enter %d rows and %d columns of numbers: ", row, column);
	
	for (int i = 0; i < matrix.row; i++)
	{
		for (int j = 0; j < matrix.column; j++)
			scanf("%lf", &(matrix.content[i][j]));
	}

	printf("------------ Matrix is£º---------------\n");
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < column; j++)
			printf("%.5f  ", matrix.content[i][j]);
		printf("\n");
	}

	return matrix;
}

Matrix CreateMatrix(int row, int column)
{
	Matrix matrix;
	int count = 0;
	
	matrix.row = row;
	matrix.column = column;

	matrix.content = (double**)calloc(row, sizeof(double*));
	for (int i = 0; i < row; i++)
		*(matrix.content + i) = (double *)calloc(column, sizeof(double));

	for ( int i = 0; i < row; i++)
		for ( int j = 0; j < column; j++)
		{	
			count++;
			matrix.content[i][j] = count;
		}
	
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < column; j++)
			printf("%.2f  ", matrix.content[i][j]);
		printf("\n");
	}
	
	return matrix;
}

Matrix CreatePhiMatrix(int N, double spindleSpeed, double samplingPeriod)
{
	Matrix phiMatrix;
	phiMatrix.row = 2 * N;
	phiMatrix.column = 2 * N;

	phiMatrix.content = (double**)calloc(phiMatrix.row, sizeof(double*));
	for (int i = 0; i < phiMatrix.row; i++)
		*(phiMatrix.content + i) = (double *)calloc(phiMatrix.column, sizeof(double));

	for (int n = 1; n <= N; n++)
	{
		phiMatrix.content[2 * (n - 1)][2 * (n - 1)] = cos(n*spindleSpeed*samplingPeriod);
		phiMatrix.content[2 * (n - 1)][2 * (n - 1) + 1] = -sin(n*spindleSpeed*samplingPeriod);
		phiMatrix.content[2 * (n - 1) + 1][2 * (n - 1)] = sin(n*spindleSpeed*samplingPeriod);
		phiMatrix.content[2 * (n - 1) + 1][2 * (n - 1) + 1] = cos(n*spindleSpeed*samplingPeriod);
	}

	printf("------ Phi Matrix------\n");
	for (int i = 0; i < phiMatrix.row; i++)
	{
		for (int j = 0; j < phiMatrix.column; j++)
			printf("%.5f  ", phiMatrix.content[i][j]);
		printf("\n");
	}

	return phiMatrix;
}

Matrix CreateHMatrix(int N)
{
	Matrix H;
	H.row = 1;
	H.column = 2 * N;

	H.content = (double**)calloc(H.row, sizeof(double*));
	for (int i = 0; i < H.row; i++)
		*(H.content + i) = (double *)calloc(H.column, sizeof(double));

	for (int i = 0; i < N; i++)
		H.content[0][2 * i] = 1;

	printf("------ H Matrix------\n");
	for (int i = 0; i < H.row; i++)
	{
		for (int j = 0; j < H.column; j++)
			printf("%.2f  ", H.content[i][j]);
		printf("\n");
	}

	printf("\n");

	return H;
}

double* ReadMeasurements(char* filePath)
{
	char data;
	int dataNumber = 0;
	int i;
	FILE* fp = fopen(filePath, "r");

	if (fp == NULL)
	{
		printf("Could not open file!\n");
		return 0;
	}

	while ( (data = fgetc(fp)) != EOF)
	{
		if (data == '\n')
			dataNumber++;
	}

	dataNumber++; // ++1 since no '\n'(return line character) for the last data

	fclose(fp);
	FILE* fp2 = fopen(filePath, "r");
	double* measurement = (double *)malloc(dataNumber * sizeof(double));

	char* test = (char *)malloc(100 * sizeof(char)); // allocate memory for a string of 100 characters
	//char test[100];

	for ( i = 0; i < dataNumber; i++)
	{
		fscanf(fp2, "%s", test);
		//printf("Test\n");
		measurement[i] = atof(test);
		//printf("%.8f\n", measurement[i]);
	}

	fclose(fp2);
	return measurement;
}

double CalculateVariance(double* measurement, int rangeBottom, int rangeTop)
{
	int numOfData = rangeTop - rangeBottom + 1;
	double sum = 0, sum1 = 0;
	double average = 0;
	double variance = 0;

	for (int i = 0; i < numOfData; i++)
		sum = sum + measurement[rangeBottom - 1 + i] * PI/30;

	average = sum / numOfData;

	for (int i = 0; i < numOfData; i++)
		sum1 = sum1 + pow((measurement[rangeBottom - 1 + i] * PI / 30 - average), 2);

	variance = sum1 / (numOfData - 1);

	printf("Variance (R) is %.10f\n", variance);
	return variance;
}

Matrix CreateQMatrix(int N, double lamda, double R)
{
	Matrix Q;

	Q.row = 2 * N;
	Q.column = 2 * N;

	Q.content = (double**)calloc(Q.row, sizeof(double*));
	for (int i = 0; i < Q.row; i++)
		*(Q.content + i) = (double *)calloc(Q.column, sizeof(double));

	for (int i = 0; i < 2 * N; i++)
		Q.content[i][i] = lamda * R;

	printf("------QMatrix------\n");
	for (int i = 0; i < 2 * N; i++)
	{
		for (int j = 0; j < 2 * N; j++)
			printf("%.15f  ", Q.content[i][j]);
		printf("\n");
	}

	return Q;
}

Matrix CreateInitialQ(int N)
{
	Matrix initialQ;

	initialQ.row = 2 * N;
	initialQ.column = 1;

	initialQ.content = (double**)calloc(initialQ.row, sizeof(double*));
	for (int i = 0; i < initialQ.row; i++)
		*(initialQ.content + i) = (double *)calloc(initialQ.column, sizeof(double));

	printf("------Q0 Matrix-----\n");
	for (int i = 0; i < initialQ.row; i++)
	{
		for (int j = 0; j < initialQ.column; j++)
			printf("%.5f  ", initialQ.content[i][j]);
		printf("\n");
	}

	printf("\n");

	return initialQ;
}

Matrix CreateInitialP(int N)
{
	Matrix initialP;

	initialP.row = 2 * N;
	initialP.column = 2 * N;
	
	initialP.content = (double**)calloc(initialP.row, sizeof(double*));
	for (int i = 0; i <  initialP.row; i++)
		*(initialP.content + i) = (double *)calloc(initialP.column, sizeof(double));

	for (int i = 0; i < initialP.row; i++)
	{
		for (int j = 0; j < initialP.column; j++)
			initialP.content[i][j] = 1;
	}

	printf("------ P0 Matrix------\n");
	for (int i = 0; i < initialP.row; i++)
	{
		for (int j = 0; j < initialP.column; j++)
			printf("%.5f  ", initialP.content[i][j]);
		printf("\n");
	}

	return initialP;
}

Matrix MatrixTranspose(Matrix matrix)
{
	Matrix matrixTransposed;
	matrixTransposed.row = matrix.column;
	matrixTransposed.column = matrix.row;

	matrixTransposed.content = (double**)calloc( matrixTransposed.row , sizeof(double*));
	for (int i = 0; i < matrixTransposed.row; i++)
		*(matrixTransposed.content + i) = (double *)calloc(matrixTransposed.column, sizeof(double));

	for (int i = 0; i < matrixTransposed.row; i++)
	{
		for (int j = 0; j < matrixTransposed.column; j++)
			matrixTransposed.content[i][j] = matrix.content[j][i];			
	}

	printf("------ Transposed Matrix------\n");
	for (int i = 0; i < matrixTransposed.row; i++)
	{
		for (int j = 0; j < matrixTransposed.column; j++)
			printf("%.5f  ", matrixTransposed.content[i][j]);
		printf("\n");
	}

	return matrixTransposed;
}

Matrix MatirxMultiplication(Matrix matrix1, Matrix matrix2) 
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

	printf("------ Multiplication Result------\n");
	PrintMatrix(result);

	return result;
}

Matrix MatirxAddition(Matrix matrix1, Matrix matrix2)
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
			result.content[i][j] = matrix1.content[i][j] + matrix2.content[i][j];
	}

	printf("------ Addition Result------\n");
	PrintMatrix(result);
	return result;
}

Matrix MatirxSubtraction(Matrix matrix1, Matrix matrix2)
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

	printf("------ Subtraction Result------\n");
	PrintMatrix(result);
	return result;
}

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

Matrix InverseMatrix(Matrix matrix)
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

	Matrix cofactorTransposed = MatrixTranspose(cofactor);

	for (int i = 0; i < matrix.row; i++)
	{
		for (int j = 0; j < matrix.column; j++)
			inverse.content[i][j] = cofactorTransposed.content[i][j] / det;
	}

	printf("------------Inversed Maitrx:----------\n");

	PrintMatrix(inverse);

	return inverse;
}