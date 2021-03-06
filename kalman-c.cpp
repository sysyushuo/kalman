﻿#include "stdafx.h"
#include <stdio.h>
#include <math.h>

#define  SIZE 4 
#define	 INTERVAL 0.0005
#define	 TIME 0.5
#define  NOISE 0.01
typedef double MAT[SIZE][SIZE], VEC[SIZE];

//void choldc1(int, MAT, VEC);

struct kalman_state
{
	VEC x;		// state vector
	VEC c;

	MAT P;  // prediction error covariance
	MAT Q;  // process noise covariance 
	MAT R;  // measurement error covariance

	MAT G;  // Kalman gain; a.k.a. K

	MAT F;  // Matrix of process model
	MAT H;  // Matrix of measurement model

	MAT Ht; // transpose of measurement Matrix
	MAT Ft; // transpose of process Matrix
	MAT Pp; // P, post-prediction, pre-update

	MAT tmp0;
	MAT tmp1;
	MAT tmp2;
	MAT tmp3;
	MAT tmp4;
	VEC tmp5;
	VEC tmp6;
	VEC tmp7;
};

void Matinv(int n,MAT  src, MAT dst)
{
	double det=0.0;

	/* Compute adjoint: */

	dst[0][0] =
		+src[1][1] * src[2][2] * src[3][3]
		- src[1][1] * src[2][3] * src[3][2]
		- src[2][1] * src[1][2] * src[3][3]
		+ src[2][1] * src[1][3] * src[3][2]
		+ src[3][1] * src[1][2] * src[2][3]
		- src[3][1] * src[1][3] * src[2][2];

	dst[0][1] =
		-src[0][1] * src[2][2] * src[3][3]
		+ src[0][1] * src[2][3] * src[3][2]
		+ src[2][1] * src[0][2] * src[3][3]
		- src[2][1] * src[0][3] * src[3][2]
		- src[3][1] * src[0][2] * src[2][3]
		+ src[3][1] * src[0][3] * src[2][2];

	dst[0][2] =
		+src[0][1] * src[1][2] * src[3][3]
		- src[0][1] * src[1][3] * src[3][2]
		- src[1][1] * src[0][2] * src[3][3]
		+ src[1][1] * src[0][3] * src[3][2]
		+ src[3][1] * src[0][2] * src[1][3]
		- src[3][1] * src[0][3] * src[1][2];

	dst[0][3] =
		-src[0][1] * src[1][2] * src[2][3]
		+ src[0][1] * src[1][3] * src[2][2]
		+ src[1][1] * src[0][2] * src[2][3]
		- src[1][1] * src[0][3] * src[2][2]
		- src[2][1] * src[0][2] * src[1][3]
		+ src[2][1] * src[0][3] * src[1][2];

	dst[1][0] =
		-src[1][0] * src[2][2] * src[3][3]
		+ src[1][0] * src[2][3] * src[3][2]
		+ src[2][0] * src[1][2] * src[3][3]
		- src[2][0] * src[1][3] * src[3][2]
		- src[3][0] * src[1][2] * src[2][3]
		+ src[3][0] * src[1][3] * src[2][2];

	dst[1][1] =
		+src[0][0] * src[2][2] * src[3][3]
		- src[0][0] * src[2][3] * src[3][2]
		- src[2][0] * src[0][2] * src[3][3]
		+ src[2][0] * src[0][3] * src[3][2]
		+ src[3][0] * src[0][2] * src[2][3]
		- src[3][0] * src[0][3] * src[2][2];

	dst[1][2] =
		-src[0][0] * src[1][2] * src[3][3]
		+ src[0][0] * src[1][3] * src[3][2]
		+ src[1][0] * src[0][2] * src[3][3]
		- src[1][0] * src[0][3] * src[3][2]
		- src[3][0] * src[0][2] * src[1][3]
		+ src[3][0] * src[0][3] * src[1][2];

	dst[1][3] =
		+src[0][0] * src[1][2] * src[2][3]
		- src[0][0] * src[1][3] * src[2][2]
		- src[1][0] * src[0][2] * src[2][3]
		+ src[1][0] * src[0][3] * src[2][2]
		+ src[2][0] * src[0][2] * src[1][3]
		- src[2][0] * src[0][3] * src[1][2];

	dst[2][0] =
		+src[1][0] * src[2][1] * src[3][3]
		- src[1][0] * src[2][3] * src[3][1]
		- src[2][0] * src[1][1] * src[3][3]
		+ src[2][0] * src[1][3] * src[3][1]
		+ src[3][0] * src[1][1] * src[2][3]
		- src[3][0] * src[1][3] * src[2][1];

	dst[2][1] =
		-src[0][0] * src[2][1] * src[3][3]
		+ src[0][0] * src[2][3] * src[3][1]
		+ src[2][0] * src[0][1] * src[3][3]
		- src[2][0] * src[0][3] * src[3][1]
		- src[3][0] * src[0][1] * src[2][3]
		+ src[3][0] * src[0][3] * src[2][1];

	dst[2][2] =
		+src[0][0] * src[1][1] * src[3][3]
		- src[0][0] * src[1][3] * src[3][1]
		- src[1][0] * src[0][1] * src[3][3]
		+ src[1][0] * src[0][3] * src[3][1]
		+ src[3][0] * src[0][1] * src[1][3]
		- src[3][0] * src[0][3] * src[1][1];

	dst[2][3] =
		-src[0][0] * src[1][1] * src[2][3]
		+ src[0][0] * src[1][3] * src[2][1]
		+ src[1][0] * src[0][1] * src[2][3]
		- src[1][0] * src[0][3] * src[2][1]
		- src[2][0] * src[0][1] * src[1][3]
		+ src[2][0] * src[0][3] * src[1][1];

	dst[3][0] =
		-src[1][0] * src[2][1] * src[3][2]
		+ src[1][0] * src[2][2] * src[3][1]
		+ src[2][0] * src[1][1] * src[3][2]
		- src[2][0] * src[1][2] * src[3][1]
		- src[3][0] * src[1][1] * src[2][2]
		+ src[3][0] * src[1][2] * src[2][1];

	dst[3][1] =
		+src[0][0] * src[2][1] * src[3][2]
		- src[0][0] * src[2][2] * src[3][1]
		- src[2][0] * src[0][1] * src[3][2]
		+ src[2][0] * src[0][2] * src[3][1]
		+ src[3][0] * src[0][1] * src[2][2]
		- src[3][0] * src[0][2] * src[2][1];

	dst[3][2] =
		-src[0][0] * src[1][1] * src[3][2]
		+ src[0][0] * src[1][2] * src[3][1]
		+ src[1][0] * src[0][1] * src[3][2]
		- src[1][0] * src[0][2] * src[3][1]
		- src[3][0] * src[0][1] * src[1][2]
		+ src[3][0] * src[0][2] * src[1][1];

	dst[3][3] =
		+src[0][0] * src[1][1] * src[2][2]
		- src[0][0] * src[1][2] * src[2][1]
		- src[1][0] * src[0][1] * src[2][2]
		+ src[1][0] * src[0][2] * src[2][1]
		+ src[2][0] * src[0][1] * src[1][2]
		- src[2][0] * src[0][2] * src[1][1];

	/* Compute determinant: */

	det = +src[0][0] * dst[0][0] + src[0][1] * dst[1][0] + src[0][2] * dst[2][0] + src[0][3] * dst[3][0];

	/* Multiply adjoint with reciprocal of determinant: */

	det = 1.0 / det;

	dst[0][0] *= det;
	dst[0][1] *= det;
	dst[0][2] *= det;
	dst[0][3] *= det;
	dst[1][0] *= det;
	dst[1][1] *= det;
	dst[1][2] *= det;
	dst[1][3] *= det;
	dst[2][0] *= det;
	dst[2][1] *= det;
	dst[2][2] *= det;
	dst[2][3] *= det;
	dst[3][0] *= det;
	dst[3][1] *= det;
	dst[3][2] *= det;
	dst[3][3] *= det;
}



//print a square real matrix A of size n with caption s
//(n items per line).
void MatPrint(const char *s, int n, MAT A) {
	int i, j;
	printf("\n %s\n", s);
	for (i = 0; i<n; i++) {
		for (j = 0; j<n; j++)
			printf(" %10.6f", A[i][j]);
		printf("\n");
	}
}

void VecPrint(const char *s, int n, VEC A) {
	int  j;
	printf("\n %s\n", s);
	for (j = 0; j<n; j++)
		printf(" %10.6f", A[j]);
	printf("\n");
}


void MatMult(int n, MAT A, MAT B, MAT C) {
	double SUM;
	int I, J, K;
	for (I = 0; I<n; I++)
		for (J = 0; J<n; J++) {
			SUM = 0.0;
			for (K = 0; K<n; K++)
				SUM += A[I][K] * B[K][J];
			C[I][J] = SUM;
		}
}

//copy MAT A in MAT A1
void MatCopy(int n, MAT A, MAT A1) {
	int i, j;
	for (i = 0; i<n; i++)
		for (j = 0; j<n; j++)
			A1[i][j] = A[i][j];
}

void Mattranspose(int n, MAT A, MAT At)
{
	int i, j;
	for (i = 0; i < n; ++i) {
		for (j = 0; j<n; ++j) {
			At[j][i] = A[i][j];
		}
	}
}

void Mataccum(int n,MAT a, MAT b )
{
	int i, j;

	for (i = 0; i < n; ++i) {
		for (j = 0; j < n; ++j) {
			a[i][j] += b[i][j];
		}
	}
}
void Matnega(int n, MAT a)
{
	int i, j;
	for (i = 0; i < n; ++i) {
		for (j = 0; j<n; ++j)
			a[i][j] = -a[i][j];
	}
}

void Mataddeye(int n, MAT a)
{
	int i,j;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			if (i == j)
			{
				a[i][j] += 1.0;
			}
		}
	}	
}

void MatVec(int n, MAT a, VEC b,VEC c)
{
	int i, j;
	double SUM=0.0;
	for (i = 0; i < n; ++i) {
		for (j = 0; j<n; ++j)
			SUM += a[i][j] * b[j];
		c[i] = SUM;
	}
}

void MatNum(int n, MAT a, double b)
{
	int i, j;
	for (i = 0; i < n; ++i) {
		for (j = 0; j<n; ++j)
			a[i][j] *= b;
	}
}

void Vecsub( int n,VEC a, VEC b,VEC c)
{
	int j;
	for (j = 0; j < n; ++j) {
		c[j] = a[j] - b[j];
	}	 
}
void Vecadd(int n, VEC  a, VEC b,VEC c)
{
	int j;
	for (j = 0; j < n; ++j) {
		c[j] = a[j] + b[j];
	}
		
}

void VecNumMult(int n, VEC a,double b,VEC c)
{
	int j;
	double SUM = 0.0;
	for (j = 0; j < n; ++j) {
		c[j] = a[j] * b;
	}
}

int kalman_init(int n, kalman_state *kalman_state_s)
{
	Mataddeye(n, kalman_state_s->R);
	Mataddeye(n, kalman_state_s->F);
	Mataddeye(n, kalman_state_s->H);
	
	kalman_state_s->c[0] = TIME*TIME/2.0;
	kalman_state_s->c[1] = TIME*TIME/2.0;
	kalman_state_s->c[2] = TIME;
	kalman_state_s->c[3] = TIME;

	kalman_state_s->F[0][2] = TIME;
	kalman_state_s->F[1][3] = TIME;

	kalman_state_s->Q[0][0] = TIME*TIME*TIME*TIME*NOISE*NOISE/4.0;
	kalman_state_s->Q[0][2] = TIME*TIME*TIME*NOISE*NOISE/2.0;
	kalman_state_s->Q[1][1] = TIME*TIME*TIME*TIME*NOISE*NOISE / 4.0;
	kalman_state_s->Q[1][3] = TIME*TIME*TIME*NOISE*NOISE / 2.0;
	kalman_state_s->Q[2][0] = TIME*TIME*TIME*NOISE*NOISE / 2.0;
	kalman_state_s->Q[2][2] = TIME*TIME*NOISE*NOISE;
	kalman_state_s->Q[3][1] = TIME*TIME*TIME*NOISE*NOISE / 2.0;
	kalman_state_s->Q[3][3] = TIME*TIME*NOISE*NOISE;

	MatCopy(n, kalman_state_s->Q, kalman_state_s->P);

	return 0;
}

int kalman_step(int n, kalman_state *kalman_state_s,VEC z)
{
	/* P_k = F_{k-1} P_{k-1} F^T_{k-1} + Q_{k-1} */
	MatMult(n, kalman_state_s->F, kalman_state_s->P, kalman_state_s->tmp0);
	Mattranspose(n, kalman_state_s->F, kalman_state_s->Ft);
	MatMult(n, kalman_state_s->tmp0, kalman_state_s->Ft, kalman_state_s->Pp);
	Mataccum(n, kalman_state_s->Pp, kalman_state_s->Q);

	/* G_k = P_k H^T_k (H_k P_k H^T_k + R)^{-1} */
	Mattranspose(n, kalman_state_s->H, kalman_state_s->Ht );
	MatMult(n, kalman_state_s->Pp, kalman_state_s->Ht, kalman_state_s->tmp1);
	MatMult(n, kalman_state_s->H, kalman_state_s->Pp, kalman_state_s->tmp2);
	MatMult(n, kalman_state_s->tmp2, kalman_state_s->Ht, kalman_state_s->tmp3 );
	Mataccum(n, kalman_state_s->tmp3, kalman_state_s->R);
	Matinv(n, kalman_state_s->tmp3, kalman_state_s->tmp4);

	//MAT tmp;
	//MatMult(n, kalman_state_s->tmp3, kalman_state_s->tmp4, tmp);
	//MatPrint("tmp3*tmp4=", n, tmp);

	MatMult(n, kalman_state_s->tmp1, kalman_state_s->tmp4, kalman_state_s->G );

	/* \hat{x}_k = \hat{x_k} + G_k(z_k - H_k*hat{x}_k)) */
	MatVec(n, kalman_state_s->F, kalman_state_s->x, kalman_state_s->tmp5);
	VecNumMult(n, kalman_state_s->c, INTERVAL, kalman_state_s->tmp7);
	Vecadd(n, kalman_state_s->tmp5, kalman_state_s->tmp7, kalman_state_s->tmp5);
	
	MatVec(n, kalman_state_s->H, kalman_state_s->tmp5, kalman_state_s->tmp6);
	Vecsub(n, z, kalman_state_s->tmp6, kalman_state_s->tmp6);
	MatVec(n, kalman_state_s->G,  kalman_state_s->tmp6, kalman_state_s->tmp7);
	Vecadd(n, kalman_state_s->tmp5, kalman_state_s->tmp7, kalman_state_s->x);
	

	/* P_k = (I - G_k H_k) P_k */
	MatMult(n, kalman_state_s->G, kalman_state_s->H, kalman_state_s->tmp0);
	Matnega(n, kalman_state_s->tmp0);
	Mataddeye(n, kalman_state_s->tmp0);
	MatMult(n, kalman_state_s->tmp0, kalman_state_s->Pp, kalman_state_s->P);

	

	/* success */
	return 0;
}

// main program to demonstrate the use of function cholsl()
int main() {
	int i;
	kalman_state kalman_state_s;
	memset(&kalman_state_s, 0, sizeof(kalman_state_s));
	VEC z = {0.0,0.0,0.1,0.1};

	kalman_init(SIZE, &kalman_state_s);

	for (i = 0; i < 10; i++)
	{
		VecPrint("x hat is", SIZE, kalman_state_s.x);
		VEC tmp = { 0.1*(1 + i),0.1*(1 + i),0.1*i,0.1*i };
		kalman_step(SIZE, &kalman_state_s, z);
		Vecadd(SIZE, z, tmp,z);
		VecPrint("z is", SIZE, z);
	}
}
