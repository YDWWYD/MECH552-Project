#pragma once
#include "Matrix.h"
#include "MatrixOperations.h"
#include <stdlib.h>
#include <stdio.h>

typedef struct {
	double SqEst;
	Matrix PPost;
	Matrix qPost;
	Matrix periodicAmp;
} KalmanOutput;

Matrix CreatePhiMatrix(int N, double spindleSpeed, double samplingPeriod);
Matrix CreateHMatrix(int N);
Matrix CreateQMatrix(int N, double lamda, Matrix R);
Matrix ReadMeasurements(char* filePath);
Matrix CalculateVariance(Matrix measurement, int rangeBottom, int rangeTop);
Matrix CreateInitial_q(int N);
Matrix CreateInitialP(int N);
Matrix KalmanFilter(char* filePath, int N, double spindleSpeed, double samplingPeriod, double lamda);

KalmanOutput KalmanAlgorithm(double measurement, Matrix Phi, Matrix PhiTranspose, Matrix H, Matrix HTranspose, Matrix identity2N, Matrix R, Matrix Q, Matrix PPost, Matrix qPost);
double KalmanFilterSingle(double measurement, Matrix Phi, Matrix PhiTranspose, Matrix H, Matrix HTranspose, Matrix identity2N, Matrix R, Matrix Q, Matrix PPost, Matrix qPost);