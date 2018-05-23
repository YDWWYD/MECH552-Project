#include "Kalman.h"
Matrix KalmanFilter(char* filePath, int N, double spindleSpeed, double samplingPeriod, double lamda)
{
	Matrix Phi = CreatePhiMatrix(N, spindleSpeed, samplingPeriod);
	Matrix H = CreateHMatrix(N);

	Matrix measurements = ReadMeasurements(filePath);
	Matrix R = CalculateVariance(measurements, 100, 700);

	Matrix Q = CreateQMatrix(N, lamda, R);
	Matrix q0 = CreateInitial_q(N);
	Matrix P0 = CreateInitialP(N);
	Matrix currentMeasurement = CreateEmptyMatrix(1, 1);

	Matrix K; // Kalman Filter Gain
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

	Matrix SpEst = CreateEmptyMatrix(measurements.row, 1);

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
		currentMeasurement.content[0][0] = measurements.content[i][0];

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
		printf("counter is %d...\n", count);

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

	//printf("----------Sp_estimation--------------\n");
	//PrintMatrix(SpEst);

	FILE* OutFilePointer = fopen("..\\Kalman Output.txt", "w");

	fprintf(OutFilePointer, "Output from the Kalman filter is:\n");
	//for (int i = 0; i < 100; i++)
	for (int i = 0; i < SpEst.row; i++)
	{
		fprintf(OutFilePointer, "%.6f\n", SpEst.content[i][0]);
	}

	fclose(OutFilePointer);

	FreeMatrixMemory(PhiTranspose);
	FreeMatrixMemory(identity2N);
	FreeMatrixMemory(HTranspose);
	//FreeMatrixMemory(SpEst);
	FreeMatrixMemory(currentMeasurement);

	FreeMatrixMemory(Phi);
	FreeMatrixMemory(H);
	FreeMatrixMemory(measurements);
	FreeMatrixMemory(R);
	FreeMatrixMemory(Q);
	FreeMatrixMemory(q0);
	FreeMatrixMemory(P0);

	return SpEst;
}
