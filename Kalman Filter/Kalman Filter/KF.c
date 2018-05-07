#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define PI acos(-1.0)

typedef struct {
	int row;
	int column;
	double** content;
} Matrix;

//double** CreateMatrix(int row, int column);
Matrix CreateMatrix(int row, int column);
double** CreatePhiMatrix(int N, double spindleSpeed, double samplingPeriod);
double* CreateHMatrix(int N);
double** CreateQMatrix(int N, double lamda, double R);
double* ReadMeasurements(char* filePath);
double CalculateVariance(double* measurement, int rangeBottom, int rangeTop);
double* CreateInitialQ(int N);
double** CreateInitialP(int N);
Matrix MatrixTranspose(Matrix matrix);

int main(void)
{
	int N = 10;
	//int i, j, count;
	double spindleSpeed = 200 * 2 * PI;
	double samplingPeriod = 3.90625E-5;
	double lamda = 1E-6;

	//double** matrix = CreateMatrix(2 * N, 1 * N);
	Matrix matrix = CreateMatrix(2 * N, 1 );
	double** phiMatrix = CreatePhiMatrix(N, spindleSpeed, samplingPeriod);
	double* HMatrix = CreateHMatrix(N);

	char* filePath = "G:\\MECH552\\yudi\\data\\Acc2.txt";

	double* meausrements = ReadMeasurements(filePath);
	double R = CalculateVariance(meausrements, 100, 700);

	double** QMatrix = CreateQMatrix(N, lamda, R);
	//double* initialQ = CreateInitialQ(N);
	//double** initialP = CreateInitialP(N);

	Matrix Phi = { 2 * N, 2 * N, phiMatrix };
	Matrix PhiTransposed = MatrixTranspose(matrix);

	system("pause");
	return 0;
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

//double** CreateMatrix(int row, int column)
//{
//	int i, j, count;
//
//	double** matrix = (double**)calloc(row, sizeof(double*));
//	for (i = 0; i < row; i++)
//		*(matrix + i) = (double *)calloc(column, sizeof(double));
//
//	count = 0;
//	for ( i = 0; i < row; i++)
//		for ( j = 0; j < column; j++)
//		{	
//			count++;
//			matrix[i][j] = count;
//		}
//
//	for (i = 0; i < row; i++)
//	{
//		for (j = 0; j < column; j++)
//			printf("%.2f  ", matrix[i][j]);
//		printf("\n");
//	}
//
//	return matrix;
//}

double** CreatePhiMatrix(int N, double spindleSpeed, double samplingPeriod)
{
	int i,j,n;

	double** phiMatrix = (double**)calloc(2*N, sizeof(double*));
	for (i = 0; i < 2*N; i++)
		*(phiMatrix + i) = (double *)calloc(2*N, sizeof(double));

	for ( n = 1; n <= N; n++)
	{
		phiMatrix[2 * (n - 1)][2 * (n - 1)] = cos(n*spindleSpeed*samplingPeriod);
		phiMatrix[2 * (n - 1)][2 * (n - 1) + 1] = -sin(n*spindleSpeed*samplingPeriod);
		phiMatrix[2 * (n - 1) + 1][2 * (n - 1)] = sin(n*spindleSpeed*samplingPeriod);
		phiMatrix[2 * (n - 1) + 1][2 * (n - 1) + 1] = cos(n*spindleSpeed*samplingPeriod);
	}

	for (i = 0; i < 2 * N; i++)
	{
		for (j = 0; j < 2 * N; j++)
			printf("%.5f  ", phiMatrix[i][j]);
		printf("\n");
	}

	return phiMatrix;
}

double* CreateHMatrix(int N)
{
	double * HMatrix = (double *)calloc(2 * N, sizeof(double));

	for (int i = 0; i < N; i++)
		HMatrix[2 * i] = 1;

	for (int i = 0; i < 2 * N; i++)
		printf("%.1f  ", HMatrix[i]);

	printf("\n");

	return HMatrix;
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

double** CreateQMatrix(int N, double lamda, double R)
{
	double** QMatrix = (double**)calloc(2 * N, sizeof(double*));
	for (int i = 0; i < 2 * N; i++)
		*(QMatrix + i) = (double *)calloc(2 * N, sizeof(double));

	for (int i = 0; i < 2 * N; i++)
		QMatrix[i][i] = lamda * R;

	printf("------QMatrix------\n");
	for (int i = 0; i < 2 * N; i++)
	{
		for (int j = 0; j < 2 * N; j++)
			printf("%.15f  ", QMatrix[i][j]);
		printf("\n");
	}

	return QMatrix;
}

double* CreateInitialQ(int N)
{
	double* initialQ = (double *)calloc(2 * N, sizeof(double));
	printf("------Q0 Matrix-----\n");
	for (int i = 0; i < 2 * N; i++)		
			printf("%.5f  ", initialQ[i]);

	printf("\n");

	return initialQ;
}

double** CreateInitialP(int N)
{
	double** initialP = (double**)calloc(2 * N, sizeof(double*));
	for (int i = 0; i < 2 * N; i++)
		*(initialP + i) = (double *)calloc(2 * N, sizeof(double));

	for (int i = 0; i < 2 * N; i++)
	{
		for (int j = 0; j < 2*N; j++)
			initialP[i][j] = 1;
	}

	printf("------ P0 Matrix------\n");
	for (int i = 0; i < 2 * N; i++)
	{
		for (int j = 0; j < 2 * N; j++)
			printf("%.5f  ", initialP[i][j]);
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

	//double** matrixTransposed = (double**)calloc( column , sizeof(double*));
	//for (int i = 0; i < row; i++)
	//	*(matrixTransposed + i) = (double *)calloc(row, sizeof(double));

	//for (int i = 0; i < column; i++)
	//{
	//	for (int j = 0; j < row; j++)
	//		matrixTransposed[i][j] = matrix[j][i];			
	//}

	printf("------ Transposed Matrix------\n");
	for (int i = 0; i < matrixTransposed.row; i++)
	{
		for (int j = 0; j < matrixTransposed.column; j++)
			printf("%.5f  ", matrixTransposed.content[i][j]);
		printf("\n");
	}

	return matrixTransposed;
}

