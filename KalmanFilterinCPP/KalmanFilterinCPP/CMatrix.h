#ifndef __CMatrix_
#define __CMatrix_

#include <iostream>
#include <iomanip>

using namespace std;

class CMatrix
{
public:
	int Row, Column;
	double* Content;

public:
	CMatrix();
	CMatrix(int row, int column);
	CMatrix(CMatrix& sourceMatrix);
	~CMatrix();
	void PrintMatrix(void);
	int Index(int row, int column);
	CMatrix* ExtractRow(CMatrix* sourceMatrix, int row);
	CMatrix* ExtractRowToColumn(CMatrix* sourceMatrix, int row);
	CMatrix operator*(const CMatrix& M1);
	CMatrix operator-(const CMatrix& M1);
	CMatrix operator+(const CMatrix& M1);
	CMatrix& operator=(const CMatrix& M1);
	//void operator=(const CMatrix& M1);
	double max(void);
};

#endif
