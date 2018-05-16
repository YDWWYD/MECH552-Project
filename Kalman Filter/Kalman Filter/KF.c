#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "Matrix.h"
#include "MatrixOperations.h"

#define PI acos(-1.0)

//typedef struct {
//	int row;
//	int column;
//	double** content;
//} Matrix;

//Matrix CreateArbitraryMatrix(int row, int column);
//Matrix CreateEmptyMatrix(int row, int column);
Matrix CreatePhiMatrix(int N, double spindleSpeed, double samplingPeriod);
Matrix CreateHMatrix(int N);
Matrix CreateQMatrix(int N, double lamda, Matrix R);
Matrix ReadMeasurements(char* filePath);
Matrix CalculateVariance(Matrix measurement, int rangeBottom, int rangeTop);
Matrix CreateInitial_q(int N);
Matrix CreateInitialP(int N);
//Matrix Transpose(Matrix matrix);
//Matrix MatrixMultiplication(Matrix matrix1, Matrix matrix2);
//Matrix MatrixAddition(Matrix matrix1, Matrix matrix2);
//Matrix MatrixSubtraction(Matrix matrix1, Matrix matrix2);
//double Determinant(Matrix matrix);
//Matrix Inverse(Matrix matrix);
//void PrintMatrix(Matrix matrix);
Matrix IdentityMatrix(int N);
Matrix Kalman(Matrix measurement, Matrix q0, Matrix P0, Matrix Q, Matrix R, Matrix Phi, Matrix H);
//void FreeMatrixMemory(Matrix matrix);

int main(void)
{
	int N = 20;
	double spindleSpeed = 200 * 2 * PI;
	double samplingPeriod = 3.90625E-5;
	double lamda = 1E-6;

	Matrix phiMatrix = CreatePhiMatrix(N, spindleSpeed, samplingPeriod);
	Matrix HMatrix = CreateHMatrix(N);

	char* filePath = "..\\Acc2.txt";

	Matrix measurements = ReadMeasurements(filePath);
	Matrix R = CalculateVariance(measurements, 100, 700);

	Matrix QMatrix = CreateQMatrix(N, lamda, R);
	Matrix initial_q = CreateInitial_q(N);
	Matrix initialP = CreateInitialP(N);

	//Matrix matrix1 = CreateArbitraryMatrix(3, 3);
	//Matrix matrix2 = CreateArbitraryMatrix(3, 3);
	//	Matrix result = MatrixMultiplication(matrix1, matrix2);

	Matrix SpEst = Kalman(measurements, initial_q, initialP, QMatrix, R, phiMatrix, HMatrix);

	FreeMatrixMemory(phiMatrix);
	FreeMatrixMemory(HMatrix);
	FreeMatrixMemory(measurements);
	FreeMatrixMemory(R);
	FreeMatrixMemory(QMatrix);
	FreeMatrixMemory(initial_q);
	FreeMatrixMemory(initialP);

	system("pause");
	return 0;
}

//Matrix CreateEmptyMatrix(int row, int column)
//{
//	Matrix matrix;
//
//	matrix.row = row;
//	matrix.column = column;
//	matrix.content = (double**)calloc(matrix.row, sizeof(double*));
//	for (int i = 0; i < matrix.row; i++)
//		*(matrix.content + i) = (double *)calloc(matrix.column, sizeof(double));
//
//	return matrix;
//}

//void PrintMatrix(Matrix matrix)
//{
//	for (int i = 0; i < matrix.row; i++)
//	{
//		for (int j = 0; j < matrix.column; j++)
//			printf("%.5f  ", matrix.content[i][j]);
//		printf("\n");
//	}
//}

Matrix IdentityMatrix(int N)
{
	Matrix identity = CreateEmptyMatrix(N, N);

	for (int i = 0; i < N; i++)
		identity.content[i][i] = 1;

	return identity;
}

//Matrix CreateArbitraryMatrix(int row, int column)
//{
//	Matrix matrix =  CreateEmptyMatrix(row, column);
//
//	printf("Please enter %d rows and %d columns of numbers: ", row, column);
//	
//	for (int i = 0; i < matrix.row; i++)
//	{
//		for (int j = 0; j < matrix.column; j++)
//			scanf("%lf", &(matrix.content[i][j]));
//	}
//
//	printf("------------ Matrix is£º---------------\n");
//	PrintMatrix(matrix);
//
//	return matrix;
//}

Matrix CreatePhiMatrix(int N, double spindleSpeed, double samplingPeriod)
{
	Matrix phiMatrix = CreateEmptyMatrix(2*N, 2*N);

	for (int n = 1; n <= N; n++)
	{
		phiMatrix.content[2 * (n - 1)][2 * (n - 1)] = cos(n*spindleSpeed*samplingPeriod);
		phiMatrix.content[2 * (n - 1)][2 * (n - 1) + 1] = -sin(n*spindleSpeed*samplingPeriod);
		phiMatrix.content[2 * (n - 1) + 1][2 * (n - 1)] = sin(n*spindleSpeed*samplingPeriod);
		phiMatrix.content[2 * (n - 1) + 1][2 * (n - 1) + 1] = cos(n*spindleSpeed*samplingPeriod);
	}

	printf("------ Phi Matrix------\n");
	PrintMatrix(phiMatrix);

	return phiMatrix;
}

Matrix CreateHMatrix(int N)
{
	Matrix H = CreateEmptyMatrix(1, 2*N);

	for (int i = 0; i < N; i++)
		H.content[0][2 * i] = 1;

	printf("------ H Matrix------\n");
	PrintMatrix(H);

	return H;
}

Matrix ReadMeasurements(char* filePath)
{
	char data;
	int dataNumber = 0;
	int i;
	FILE* fp = fopen(filePath, "r");

	if (fp == NULL)
	{
		printf("Could not open file!\n");
		return;
	}

	while ( (data = fgetc(fp)) != EOF)
	{
		if (data == '\n')
			dataNumber++;
	}

	dataNumber++; // ++1 since no '\n'(return line character) for the last data

	fclose(fp);
	FILE* fp2 = fopen(filePath, "r");
	Matrix measurement = CreateEmptyMatrix(dataNumber,1);

	char* test = (char *)malloc(100 * sizeof(char)); // allocate memory for a string of 100 characters

	for ( i = 0; i < dataNumber; i++)
	{
		fscanf(fp2, "%s", test);
		measurement.content[i][0] = atof(test);
	}

	fclose(fp2);
	return measurement;
}

Matrix CalculateVariance(Matrix measurement, int rangeBottom, int rangeTop)
{
	int numOfData = rangeTop - rangeBottom + 1;
	double sum = 0, sum1 = 0;
	double average = 0;

	Matrix variance = CreateEmptyMatrix(1, 1);

	for (int i = 0; i < numOfData; i++)
		sum = sum + measurement.content[rangeBottom - 1 + i][0] * PI/30;

	average = sum / numOfData;

	for (int i = 0; i < numOfData; i++)
		sum1 = sum1 + pow((measurement.content[rangeBottom - 1 + i][0] * PI / 30 - average), 2);

	variance.content[0][0] = sum1 / (numOfData - 1);

	printf("Variance (R) is %.10f\n", variance.content[0][0]);
	return variance;
}

Matrix CreateQMatrix(int N, double lamda, Matrix R)
{
	Matrix Q = CreateEmptyMatrix(2*N, 2*N);

	for (int i = 0; i < 2 * N; i++)
		Q.content[i][i] = lamda * R.content[0][0];

	//printf("------QMatrix------\n");
	//for (int i = 0; i < 2 * N; i++)
	//{
	//	for (int j = 0; j < 2 * N; j++)
	//		printf("%.15f  ", Q.content[i][j]);
	//	printf("\n");
	//}

	return Q;
}

Matrix CreateInitial_q(int N)
{
	Matrix initial_q =  CreateEmptyMatrix(2*N, 1);

	printf("------q0 Matrix-----\n");
	PrintMatrix(initial_q);

	return initial_q;
}

Matrix CreateInitialP(int N)
{
	Matrix initialP = CreateEmptyMatrix(2*N, 2*N);

	for (int i = 0; i < initialP.row; i++)
			initialP.content[i][i] = 10;

	printf("------ P0 Matrix------\n");
	PrintMatrix(initialP);

	return initialP;
}

//Matrix Transpose(Matrix matrix)
//{
//	Matrix matrixTransposed = CreateEmptyMatrix(matrix.column, matrix.row);
//
//	for (int i = 0; i < matrixTransposed.row; i++)
//	{
//		for (int j = 0; j < matrixTransposed.column; j++)
//			matrixTransposed.content[i][j] = matrix.content[j][i];			
//	}
//
//	//printf("------ Transposed Matrix------\n");
//	//PrintMatrix(matrixTransposed);
//
//	return matrixTransposed;
//}

//Matrix MatrixMultiplication(Matrix matrix1, Matrix matrix2) 
//{
//	if (matrix1.column != matrix2.row)
//	{
//		printf("Dimension Invalid!\n");
//		return;
//	}
//
//	Matrix result = CreateEmptyMatrix(matrix1.row, matrix2.column);
//
//	for (int i = 0; i < result.row; i++)
//	{
//		for (int j = 0; j < result.column; j++)
//		{
//			//result[i, j] = 0; //Default values are 0
//			for (int k = 0; k < matrix1.column; k++)
//				result.content[i][j] = result.content[i][j] + matrix1.content[i][k] * matrix2.content[k][j];
//		}
//	}
//
//	//printf("------ Multiplication Result------\n");
//	//PrintMatrix(result);
//
//	return result;
//}

//Matrix MatrixAddition(Matrix matrix1, Matrix matrix2)
//{
//	if (matrix1.row != matrix2.row || matrix1.column != matrix2.column)
//	{
//		printf("Dimension Invalid!\n");
//		return;
//	}
//
//	Matrix result = CreateEmptyMatrix(matrix1.row, matrix1.column);
//
//	for (int i = 0; i < matrix1.row; i++)
//	{
//		for (int j = 0; j < matrix1.column; j++)
//			result.content[i][j] = matrix1.content[i][j] + matrix2.content[i][j];
//	}
//
//	//printf("------ Addition Result------\n");
//	//PrintMatrix(result);
//	return result;
//}

//Matrix MatrixSubtraction(Matrix matrix1, Matrix matrix2)
//{
//	if (matrix1.row != matrix2.row || matrix1.column != matrix2.column)
//	{
//		printf("Dimension Invalid!\n");
//		return;
//	}
//
//	Matrix result = CreateEmptyMatrix(matrix1.row, matrix1.column);
//
//	for (int i = 0; i < matrix1.row; i++)
//	{
//		for (int j = 0; j < matrix1.column; j++)
//			result.content[i][j] = matrix1.content[i][j] - matrix2.content[i][j];
//	}
//
//	//printf("------ Subtraction Result------\n");
//	//PrintMatrix(result);
//	return result;
//}

//double Determinant(Matrix matrix)
//{
//	double det = 0, coefficient = 1;
//
//	if (matrix.row != matrix.column)
//	{
//		printf("Invalid Dimension. Determinant does not exist!\n");
//		return;
//	}
//
//	if (matrix.row == 1)
//		det = matrix.content[0][0];
//	else
//	{
//		Matrix minor = CreateEmptyMatrix(matrix.row - 1, matrix.column - 1);
//		for (int c = 0; c < matrix.row; c++) // Loop through the 1st row
//		{
//			int m = 0;
//			int n = 0;
//
//			for (int i = 0; i < matrix.row; i++)
//			{
//				for (int j = 0; j < matrix.row; j++)
//				{
//					if (i != 0 && j != c) // Exclude the current row and column of elements
//					{
//						minor.content[m][n] = matrix.content[i][j];
//						if (n < matrix.column - 2)
//							n++;
//						else
//						{
//							n = 0;
//							m++;
//						}
//					}
//				}
//			}
//
//			det = det + coefficient*(matrix.content[0][c] * Determinant(minor));
//			coefficient = -1 * coefficient;
//		}
//	}
//	
//	return det;
//}

//Matrix Inverse(Matrix matrix)
//{
//	Matrix inverse = CreateEmptyMatrix(matrix.row, matrix.column);
//	Matrix cofactor = CreateEmptyMatrix(matrix.row, matrix.column);
//	Matrix minor = CreateEmptyMatrix(matrix.row - 1, matrix.column - 1);
//
//	double det = Determinant(matrix);
//
//	if (det == 0)
//	{
//		printf("Singular matrix, no inverse!\n");
//		return;
//	}
//
//	if (matrix.row == 1) // matrix is a number -> inverse = 1/det
//	{
//		inverse.content[0][0] = 1 / det;
//	}
//	
//	else
//	{
//		for (int p = 0; p < matrix.row; p++)
//		{
//			for (int q = 0; q < matrix.column; q++)
//			{
//				int m = 0;
//				int n = 0;
//
//				for (int i = 0; i < matrix.row; i++)
//				{
//					for (int j = 0; j < matrix.column; j++)
//					{
//						if (i != p && j != q)
//						{
//							minor.content[m][n] = matrix.content[i][j];
//							if (n < matrix.column - 2)
//								n++;
//							else
//							{
//								n = 0;
//								m++;
//							}
//						}
//					}
//				}
//				cofactor.content[p][q] = pow(-1, p + q) * Determinant(minor);
//			}
//		}
//
//		Matrix cofactorTransposed = Transpose(cofactor);
//
//		for (int i = 0; i < matrix.row; i++)
//		{
//			for (int j = 0; j < matrix.column; j++)
//				inverse.content[i][j] = cofactorTransposed.content[i][j] / det;
//		}
//	}
//
//	//printf("------------Inversed Maitrx:----------\n");
//	//PrintMatrix(inverse);
//
//	return inverse;
//}

Matrix Kalman(Matrix measurement, Matrix q0, Matrix P0, Matrix Q, Matrix R, Matrix Phi, Matrix H) 
{
	Matrix currentMeasurement = CreateEmptyMatrix(1, 1);

	Matrix K ; // Kalman Filter Gain
	Matrix qHatPrior; // q prior estimate [2N * 1]
	Matrix PHatPrior; // P prior estimate [2N * 2N]
	Matrix qPost; // q posteriori estimate [2N * 1]
	Matrix PPost; // P posteriori estimate [2N * 2N]

	Matrix PhiTranspose = Transpose(Phi);
	Matrix identity2N = IdentityMatrix(P0.row); // Identity [2N * 2N]
	Matrix HTranspose = Transpose(H);

	Matrix PPhiTrans; // PPhiTrans: P_k-1*Phi' [2N * 2N]
	Matrix PhiPPhiTrans; //PhiPPhiTrans: Phi*P_k-1*Phi' [2N * 2N]

	Matrix PHTrans; // PHTrans: P_k*H' [2N * 1]
	Matrix HPHTrans; //HPHTrans: H*P_k*H' [1 * 1]
	Matrix HPHTransPlusR; //HPHTrans: H*P_k*H' + R [1 * 1]
	Matrix invHPHTransPlusR; //invHPHTrans: inv(H*P_k*H' + R) [1 * 1]
	Matrix HTransInv; //HTransInv: H'*inv(H*P_k*H' + R) [2N*1]

	Matrix Hq; // Hq: H*q_k [1 * 1]
	Matrix sMinusHq; // sMinusHq: sk- H*q_k  [1 * 1]
	Matrix KsMinusHq; // KsMinusHq: K_k*(s_k- H*q_k) [2N * 1]

	Matrix KH; // KH: K_k*H [2N*2N]
	Matrix IminusKH; // IminusKH: I-K_k*H [2N*2N]

	Matrix Hq_post; // Hq: H*q_k [1 * 1]

	Matrix SpEst = CreateEmptyMatrix(measurement.row, 1);

	int count = 0;

	for (int i = 0; i < 100; i++)
	//for (int i = 0; i < measurement.row; i++)
	{
		if (i == 0)
		{
			qHatPrior = MatrixMultiplication(Phi, q0);
			//printf("----------------qHatPrior:----------------\n");
			//PrintMatrix(qHatPrior);

			PPhiTrans = MatrixMultiplication(P0, PhiTranspose); // PPhiTrans: P_k-1*Phi'
			PhiPPhiTrans = MatrixMultiplication(Phi, PPhiTrans);
			PHatPrior = MatrixAddition(PhiPPhiTrans, Q);
			//printf("----------------PHatPrior:----------------\n");
			//for (int i = 0; i < PHatPrior.row; i++)
			//{
			//	for (int j = 0; j < PHatPrior.column; j++)
			//		printf("%.25f  ", PHatPrior.content[i][j]);
			//	printf("\n");
			//}
		}

		else
		{
			qHatPrior = MatrixMultiplication(Phi, qPost);
			
			PPhiTrans = MatrixMultiplication(PPost, PhiTranspose); 
			PhiPPhiTrans = MatrixMultiplication(Phi, PPhiTrans);
			PHatPrior = MatrixAddition(PhiPPhiTrans, Q);

			FreeMatrixMemory(qPost);
			FreeMatrixMemory(PPost);
		}

		PHTrans = MatrixMultiplication(PHatPrior, HTranspose);
		HPHTrans = MatrixMultiplication(H, PHTrans);
		HPHTransPlusR = MatrixAddition(HPHTrans, R); //HPHTrans: H*P_k*H' + R [1 * 1]
		invHPHTransPlusR = Inverse(HPHTransPlusR); //invHPHTrans: inv(H*P_k*H' + R) [1 * 1]
		HTransInv = MatrixMultiplication(HTranspose, invHPHTransPlusR);

		K = MatrixMultiplication(PHatPrior, HTransInv);
		//printf("----------------Before Inverse:----------------\n");
		//PrintMatrix(MatrixMultiplication(H, MatrixMultiplication(PHatPrior, Transpose(H))));
		//printf("----------------K:----------------\n");
		//PrintMatrix(K);
		currentMeasurement.content[0][0] = measurement.content[i][0];

		Hq = MatrixMultiplication(H, qHatPrior); // Hq: H*q_k [1 * 1]
		sMinusHq = MatrixSubtraction(currentMeasurement, Hq); // sMinusHq: sk- H*q_k  [1 * 1]
		KsMinusHq = MatrixMultiplication(K, sMinusHq); // KsMinusHq: Kk*(sk- H*q_k) [2N * 1]

		qPost = MatrixAddition(qHatPrior, KsMinusHq);

		KH = MatrixMultiplication(K, H); // KH: K_k*H [2N*2N]
		IminusKH = MatrixSubtraction(identity2N, KH); // IminusKH: I-K_k*H [2N*2N]

		PPost = MatrixMultiplication(IminusKH, PHatPrior);

		Hq_post = MatrixMultiplication(H, qPost);
		SpEst.content[i][0] = Hq_post.content[0][0];

		count++;
		printf("counter is %d...\n",count);

		FreeMatrixMemory(PPhiTrans);
		FreeMatrixMemory(PhiPPhiTrans);

		FreeMatrixMemory(PHTrans);
		FreeMatrixMemory(HPHTrans);
		FreeMatrixMemory(HPHTransPlusR);
		FreeMatrixMemory(invHPHTransPlusR);
		FreeMatrixMemory(HTransInv);

		FreeMatrixMemory(Hq);
		FreeMatrixMemory(sMinusHq);
		FreeMatrixMemory(KsMinusHq);
		FreeMatrixMemory(KH);
		FreeMatrixMemory(IminusKH);
		FreeMatrixMemory(Hq_post);

		FreeMatrixMemory(K);
		FreeMatrixMemory(qHatPrior);
		FreeMatrixMemory(PHatPrior);
	}
	//qHatPrior = MatrixMultiplication(Phi, q0);
	//PHatPrior = MatrixAddition(MatrixMultiplication(Phi, MatrixMultiplication(P0, Transpose(Phi))), Q);
	//K = MatrixMultiplication(PHatPrior, MatrixMultiplication(Transpose(H), Inverse(MatrixAddition(MatrixMultiplication(H, MatrixMultiplication(PHatPrior, Transpose(H))), R))));
	//currentMeasurement.content[0][0] = measurement.content[0][0];
	//qPost = MatrixAddition(qHatPrior, MatrixMultiplication(K, MatrixSubtraction(currentMeasurement, MatrixMultiplication(H, qHatPrior))));
	//PPost = MatrixMultiplication(MatrixSubtraction(IdentityMatrix(P0.row), MatrixMultiplication(K, H)), PHatPrior);

	printf("----------Sp_estimation--------------\n");
	PrintMatrix(SpEst);

	FILE* OutFilePointer = fopen("..\\Kalman Output.txt", "w");

	fprintf(OutFilePointer, "Output from the Kalman filter is:\n");
	for (int i = 0; i < SpEst.row; i++)
	{
		fprintf(OutFilePointer, "%.6f\n", SpEst.content[i][0]);
	}

	fclose(OutFilePointer);

	FreeMatrixMemory(PhiTranspose);
	FreeMatrixMemory(identity2N);
	FreeMatrixMemory(HTranspose);
	FreeMatrixMemory(SpEst);
	FreeMatrixMemory(currentMeasurement);

	return SpEst;
}

//void FreeMatrixMemory(Matrix matrix)
//{
//	for (int i = 0; i < matrix.row; i++)
//	{
//		free(matrix.content[i]);
//		matrix.content[i] = NULL;
//	}
//
//	free(matrix.content);
//	matrix.content = NULL;
//}