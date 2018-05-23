#include "Kalman.h"

KalmanOutput KalmanAlgorithm(double measurement, Matrix Phi, Matrix PhiTranspose, Matrix H, Matrix HTranspose, Matrix identity2N, Matrix R, Matrix Q, Matrix PPost, Matrix qPost)
{
	Matrix currentMeasurement = CreateEmptyMatrix(1, 1);
	currentMeasurement.content[0][0] = measurement;

	KalmanOutput KalmanOut;
	Matrix K; // Kalman Filter Gain
	Matrix qHatPrior; // q prior estimate [2N * 1]
	Matrix PHatPrior; // P prior estimate [2N * 2N]
	//Matrix qPost; // q posteriori estimate [2N * 1]
	//Matrix PPost; // P posteriori estimate [2N * 2N]

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

	qHatPrior = MatrixMultiplication(Phi, qPost);

	PPhiTrans = MatrixMultiplication(PPost, PhiTranspose);
	PhiPPhiTrans = MatrixMultiplication(Phi, PPhiTrans);
	PHatPrior = MatrixAddition(PhiPPhiTrans, Q);

	PHTrans = MatrixMultiplication(PHatPrior, HTranspose);
	HPHTrans = MatrixMultiplication(H, PHTrans);
	HPHTransPlusR = MatrixAddition(HPHTrans, R); //HPHTrans: H*P_k*H' + R [1 * 1]
	invHPHTransPlusR = Inverse(HPHTransPlusR); //invHPHTrans: inv(H*P_k*H' + R) [1 * 1]
	HTransInv = MatrixMultiplication(HTranspose, invHPHTransPlusR);

	K = MatrixMultiplication(PHatPrior, HTransInv);

	Hq = MatrixMultiplication(H, qHatPrior); // Hq: H*q_k [1 * 1]
	sMinusHq = MatrixSubtraction(currentMeasurement, Hq); // sMinusHq: sk- H*q_k  [1 * 1]
	KsMinusHq = MatrixMultiplication(K, sMinusHq); // KsMinusHq: Kk*(sk- H*q_k) [2N * 1]

	qPost = MatrixAddition(qHatPrior, KsMinusHq);

	KH = MatrixMultiplication(K, H); // KH: K_k*H [2N*2N]
	IminusKH = MatrixSubtraction(identity2N, KH); // IminusKH: I-K_k*H [2N*2N]

	PPost = MatrixMultiplication(IminusKH, PHatPrior);

	Hq_post = MatrixMultiplication(H, qPost);

	KalmanOut.SqEst = Hq_post.content[0][0];
	KalmanOut.qPost = qPost;
	KalmanOut.PPost = PPost;

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
	//FreeMatrixMemory(Hq_post);

	FreeMatrixMemory(K);
	FreeMatrixMemory(qHatPrior);
	FreeMatrixMemory(PHatPrior);

	FreeMatrixMemory(currentMeasurement);

	return KalmanOut;
	//SpEst.content[i][0] = Hq_post.content[0][0];
}