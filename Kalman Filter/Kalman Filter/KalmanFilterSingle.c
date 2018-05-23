#include "Kalman.h"

double KalmanFilterSingle(double measurement, Matrix Phi, Matrix PhiTranspose, Matrix H, Matrix HTranspose, Matrix identity2N, Matrix R, Matrix Q, Matrix PPost, Matrix qPost)
{
	KalmanOutput KalmanOut = KalmanAlgorithm(measurement, Phi, PhiTranspose, H, HTranspose, identity2N, R, Q, PPost, qPost);

	FreeMatrixMemory(PPost);
	FreeMatrixMemory(qPost);

	PPost = KalmanOut.PPost;
	qPost = KalmanOut.qPost;

	return KalmanOut.SqEst;
}