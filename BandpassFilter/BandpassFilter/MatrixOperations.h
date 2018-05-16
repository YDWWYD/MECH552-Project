#pragma once
#include "Matrix.h"

Matrix CreateArbitraryMatrix(int row, int column);
Matrix CreateEmptyMatrix(int row, int column);
//Matrix CreatePhiMatrix(int N, double spindleSpeed, double samplingPeriod);
//Matrix CreateHMatrix(int N);
//Matrix CreateQMatrix(int N, double lamda, Matrix R);
//Matrix ReadMeasurements(char* filePath);
//Matrix CalculateVariance(Matrix measurement, int rangeBottom, int rangeTop);
//Matrix CreateInitial_q(int N);
//Matrix CreateInitialP(int N);
Matrix Transpose(Matrix matrix);
Matrix MatrixMultiplication(Matrix matrix1, Matrix matrix2);
Matrix MatrixAddition(Matrix matrix1, Matrix matrix2);
Matrix MatrixSubtraction(Matrix matrix1, Matrix matrix2);
double Determinant(Matrix matrix);
Matrix Inverse(Matrix matrix);
void PrintMatrix(Matrix matrix);
Matrix IdentityMatrix(int N);
//Matrix Kalman(Matrix measurement, Matrix q0, Matrix P0, Matrix Q, Matrix R, Matrix Phi, Matrix H);
void FreeMatrixMemory(Matrix matrix);
