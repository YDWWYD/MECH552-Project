#include "CBandpassFilters.h"
#include "buttordBandpassOnly.h"
#include "butterBandpassOnly.h"

CBandpassFilters::CBandpassFilters(double StartFreq, double StopFreq, double TPE, double dt, double StopBandMag, double RippleMag, double Ts)
{
	double PassBand[2];
	double StopBand[2];
	double SamplingFreq = (double)1 / Ts;	
	double wn;
	double Wn[2];

	double* order = new double[1];
	double* order1 = new double[1];
	emxArray_real_T *num;
	emxArray_creal_T *den;
	emxInit_real_T1(&num, 2);
	emxInit_creal_T1(&den, 2);

	StartBand = (int)ceil(StartFreq / TPE);
	NumberOfFilters = (int)((ceil(StopFreq / TPE))*TPE / TPE - StartBand);
	//NumberOfFilters = (int)((ceil(StopFreq / TPE))*TPE / TPE - StartBand) - 1;

	for (int i = 0; i < NumberOfFilters; i++)
	{
		wn = StartBand + i;
		PassBand[0] = (wn + dt)*TPE / (SamplingFreq / 2);
		PassBand[1] = (wn + 1 - dt)*TPE / (SamplingFreq / 2);

		StopBand[0] = wn*TPE / (SamplingFreq / 2);
		StopBand[1] = (wn + 1)*TPE / (SamplingFreq / 2);

		if (i == 0)
		{
			buttordBandpassOnly(PassBand, StopBand, RippleMag, StopBandMag, order, Wn);
			Order = (int)*order;
			cout << "Order is " << *order << "\n";
			Numerator = new CMatrix(NumberOfFilters, 2 * Order + 1);
			Denominator = new CMatrix(NumberOfFilters, 2 * Order);
		}

		buttordBandpassOnly(PassBand, StopBand, RippleMag, StopBandMag, order1, Wn);
		butterBandpassOnly(*order, Wn, num, den);

		for (int j = 0; j < Numerator->Column; j++)
			Numerator->Content[Numerator->Index(i, j)] = num->data[j];

		for (int k = 0; k < Denominator->Column; k++)
			Denominator->Content[Denominator->Index(i, k)] = den->data[k + 1].re;
	}

	E = new CMatrix(NumberOfFilters, 2 * Order + 1);
	EDelayed = new CMatrix(NumberOfFilters, 2 * Order);
	EDelayedSingle = new CMatrix(2 * Order, 1);
	ASingle = new CMatrix(1, 2 * Order);
	BSingle = new CMatrix(1, 2 * Order + 1);
	ESingleTrans = new CMatrix(2 * Order + 1, 1);
	BandpassOutputs = new CMatrix(NumberOfFilters, 1);

	delete order;
	delete order1;
}

CBandpassFilters:: ~CBandpassFilters()
{
	delete Numerator;
	delete Denominator;
	delete E;
	delete EDelayed;
	delete EDelayedSingle;
	delete ASingle;
	delete BSingle;
	delete ESingleTrans;
	delete BandpassOutputs;
}

double CBandpassFilters::CalculateInput(double measurement, double estimation)
{
	return measurement - estimation;
}

void CBandpassFilters::RunBandpassFilters(double input)
{
	// ----- For loop to calculate A*(E_delayed)' of each filter and update 1st column of E signal
	for (int i = 0; i < NumberOfFilters; i++)
	{
		EDelayedSingle->ExtractRowToColumn(EDelayed, i); // Ed': [2n x 1]
		ASingle->ExtractRow(Denominator, i); // Extract A of #(i+1) filter: [1 x 2n]
		E->Content[i*E->Column] = input - (*ASingle * *EDelayedSingle).Content[0];// (INPUT - A*Ed') of #i filter [1x1]
		//AETrans->Content[i] = (*ASingle * *EDelayed).Content[0]; // A*Ed' of #i filter [1x1]
	}

	// Update E signal of the remaining columns: copy from EDelayed signal
	for (int i = 0; i < NumberOfFilters; i++) // # of filters
	{
		for (int j = 0; j < EDelayed->Column; j++)
		{
			E->Content[i*E->Column + j + 1] = EDelayed->Content[i*EDelayed->Column + j];
		}
	}

	// ---- for loop to calculate each filter's output---------
	for (int i = 0; i < NumberOfFilters; i++)
	{
		ESingleTrans->ExtractRowToColumn(E, i); // E': [2n+1 x 1]
		BSingle->ExtractRow(Numerator, i); // B of #(i+1) filter: [1 x 2n+1]
		//Matrix signalOutputSingle = MatrixMultiplication(BSingle, ESingleTrans); // Output of #(i+1) filter: [1 x 1]
		BandpassOutputs->Content[i] = (*BSingle**ESingleTrans).Content[0];// Output of #(i+1) filter: [1 x 1]
	}

	// Update delayed E
	for (int i = 0; i < NumberOfFilters; i++)
	{
		for (int j = 0; j < EDelayed->Column - 1; j++)
		{
			int index1 = i*EDelayed->Column + EDelayed->Column - 1 - j;
			EDelayed->Content[index1] = EDelayed->Content[index1 - 1];
		}
		EDelayed->Content[i*EDelayed->Column] = E->Content[i*E->Column];
	}
}