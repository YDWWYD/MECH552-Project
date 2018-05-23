#include "Kalman.h"

Matrix CalculateVariance(Matrix measurement, int rangeBottom, int rangeTop)
{
	int numOfData = rangeTop - rangeBottom + 1;
	double sum = 0, sum1 = 0;
	double average = 0;

	Matrix variance = CreateEmptyMatrix(1, 1);

	for (int i = 0; i < numOfData; i++)
		sum = sum + measurement.content[rangeBottom - 1 + i][0] * PI / 30;

	average = sum / numOfData;

	for (int i = 0; i < numOfData; i++)
		sum1 = sum1 + pow((measurement.content[rangeBottom - 1 + i][0] * PI / 30 - average), 2);

	variance.content[0][0] = sum1 / (numOfData - 1);

	//printf("Variance (R) is %.10f\n", variance.content[0][0]);
	return variance;
}