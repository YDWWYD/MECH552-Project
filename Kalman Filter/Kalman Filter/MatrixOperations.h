#pragma once
#include "Matrix.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define PI acos(-1.0)

Matrix CreateArbitraryMatrix(int row, int column);
Matrix CreateEmptyMatrix(int row, int column);
Matrix Transpose(Matrix matrix);
Matrix MatrixMultiplication(Matrix matrix1, Matrix matrix2);
Matrix MatrixAddition(Matrix matrix1, Matrix matrix2);
Matrix MatrixSubtraction(Matrix matrix1, Matrix matrix2);
double Determinant(Matrix matrix);
Matrix Inverse(Matrix matrix);
void PrintMatrix(Matrix matrix);
Matrix IdentityMatrix(int N);
void FreeMatrixMemory(Matrix matrix);
Matrix ExtractRow(Matrix matrix, int row);
Matrix ExtractRowToColumn(Matrix inputMatrix, int tow);