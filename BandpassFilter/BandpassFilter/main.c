//#include "rt_nonfinite.h"
#include "buttordBandpassOnly.h"
//#include "sortButtord.h"
#include "stdio.h"
//#include "rt_nonfinite.h"
#include "butterBandpassOnly.h"
//#include "butterBandpassOnly_emxutil.h"
//#include "exp.h"
//#include "isequal.h"
//#include "sortButter.h"
//#include "relop.h"
//#include "mypoly.h"
//#include "mybilinear.h"
//#include "mylp2bs.h"
//#include "myzp2ss.h"
//#include "butterBandpassOnly_rtwutil.h"
#include "Matrix.h"
#include "MatrixOperations.h"
#include "Math.h"

typedef struct BandpassCoefficients
{
	Matrix Numerator;
	Matrix Denominator;
}BandpassCoefficients;

BandpassCoefficients BandpassFilterCofficients(double StartFreq, double StopFreq, double TPE, double dt, double StopBandMag, double RippleMag, double Ts);

int main(void)
{
	/*double wp[] = { 0.3781, 0.4031 };
	double ws[] = { 0.3750, 0.4063 };
	double rp = 2, rs = 5;
	double n = 1;
	double* order = (double *)malloc(1 * sizeof(double));
	double wn[2];
	buttordBandpassOnly(wp, ws, rp, rs, order, wn);

	emxArray_real_T *num;
	emxArray_creal_T *den;
	emxInit_real_T1(&num, 2);
	emxInit_creal_T1(&den, 2);

	butterBandpassOnly(*order, wn, num, den);

	for (int i = 0; i < num->size[1]; i++)
	{
		printf("num[%d] is %.10f\n", i, num->data[i]);
	}
		
	for (int i = 0; i < den->size[1]; i++)
	{
		printf("Re(den[%d]) is %.8f		Im(den[%d]) is %f\n", i, den->data[i].re, i, den->data[i].im);
	}

	free(order);
	order = NULL;*/
	//BandpassCoefficients bandpassCofficients = 
	//Matrix denominatorCoeffi = CreateEmptyMatrix(1, 1);
	//Matrix numeratorCoeffi = CreateEmptyMatrix(1, 1);

	double startFreq = 1000; //Hz
	double stopFreq = 6900; 
	double TPE = 400;
	double stopBandMag = 5; //dB
	double rippleMag = 2;
	double Ts = (double)1 / 25600;
	//printf("Ts is %f\n", Ts);
	double dt = 0.1;

	BandpassCoefficients bandpassCofficients = BandpassFilterCofficients(startFreq, stopFreq, TPE, dt, stopBandMag, rippleMag, Ts);
	printf("----------Denominator----------------\n");
	PrintMatrix(bandpassCofficients.Denominator);
	printf("----------Numerator----------------\n");
	PrintMatrix(bandpassCofficients.Numerator);

	system("pause");

	return 0;
}

BandpassCoefficients BandpassFilterCofficients(double StartFreq, double StopFreq, double TPE, double dt, double StopBandMag, double RippleMag, double Ts)
{
	BandpassCoefficients bandpassCoefficients;
	double PassBand[2];
	double StopBand[2];
	double SamplingFreq = (double)1 / Ts;
	double m1 = ceil(StartFreq / TPE);
	double m = (ceil(StopFreq / TPE))*TPE / TPE - m1;
	double wn;
	double Wn[2];
	double* order = (double *)malloc(1 * sizeof(double));
	double* order1 = (double *)malloc(1 * sizeof(double));
	emxArray_real_T *num;
	emxArray_creal_T *den;
	emxInit_real_T1(&num, 2);
	emxInit_creal_T1(&den, 2);

	for (int i = 0; i < m; i++)
	{
		wn = m1 + i;
		PassBand[0] = (wn + dt)*TPE/(SamplingFreq/2);
		PassBand[1] = (wn + 1 - dt)*TPE/(SamplingFreq/2);

		StopBand[0] = wn*TPE / (SamplingFreq / 2);
		StopBand[1] = (wn + 1)*TPE / (SamplingFreq / 2);

		if (i == 0)
		{
			buttordBandpassOnly(PassBand, StopBand, RippleMag, StopBandMag, order, Wn);
			printf("Order is %f\n", *order);
			bandpassCoefficients.Numerator = CreateEmptyMatrix((int)m, 2 * (int)*order + 1);
			bandpassCoefficients.Denominator = CreateEmptyMatrix((int)m, 2 * (int)*order);
		}

		buttordBandpassOnly(PassBand, StopBand, RippleMag, StopBandMag, order1, Wn);
		butterBandpassOnly(*order, Wn, num, den);

		for (int j = 0; j < bandpassCoefficients.Numerator.column; j++)
			bandpassCoefficients.Numerator.content[i][j] = num->data[j];

		for (int k = 0; k < bandpassCoefficients.Denominator.column; k++)
			bandpassCoefficients.Denominator.content[i][k] = den->data[k+1].re;
	}

	return bandpassCoefficients;
}


