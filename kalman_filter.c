
#include "kalman_filter.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

typedef struct {
	
    double * x;    /* state vector */
	
    double * P;  /* prediction error covariance */ 
    double * Q;  /* process noise covariance */
    double * R;  /* measurement error covariance */
	
    double * G;  /* Kalman gain; a.k.a. K */
	
    double * F;  /* Jacobian of process model */
    double * H;  /* Jacobian of measurement model */
	
   double * Ht; /* transpose of measurement Jacobian */
    double * Ft; /* transpose of process Jacobian */
    double * Pp; /* P, post-prediction, pre-update */

	/* temporary storage */
    double * tmp0;
    double * tmp1;
    double * tmp2;
    double * tmp3;
    double * tmp4;
    double * tmp5; 

} kalman_state_t;

static int choldc1(double * a, double * p, int n) {
    int i,j,k;
    double sum;

    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            sum = a[i*n+j];
            for (k = i - 1; k >= 0; k--) {
                sum -= a[i*n+k] * a[j*n+k];
            }
            if (i == j) {
                if (sum <= 0) {
                    return 1; /* error */
                }
                p[i] = sqrt(sum);
            }
            else {
                a[j*n+i] = sum / p[i];
            }
        }
    }

    return 0; /* success */
}

static int choldcsl(double * A, double * a, double * p, int n) 
{
    int i,j,k; double sum;
    for (i = 0; i < n; i++) 
        for (j = 0; j < n; j++) 
            a[i*n+j] = A[i*n+j];
    if (choldc1(a, p, n)) return 1;
    for (i = 0; i < n; i++) {
        a[i*n+i] = 1 / p[i];
        for (j = i + 1; j < n; j++) {
            sum = 0;
            for (k = i; k < j; k++) {
                sum -= a[j*n+k] * a[k*n+i];
            }
            a[j*n+i] = sum / p[j];
        }
    }

    return 0; /* success */
}


static int cholsl(double * A, double * a, double * p, int n) 
{
    int i,j,k;
    if (choldcsl(A,a,p,n)) return 1;
    for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) {
            a[i*n+j] = 0.0;
        }
    }
    for (i = 0; i < n; i++) {
        a[i*n+i] *= a[i*n+i];
        for (k = i + 1; k < n; k++) {
            a[i*n+i] += a[k*n+i] * a[k*n+i];
        }
        for (j = i + 1; j < n; j++) {
            for (k = j; k < n; k++) {
                a[i*n+j] += a[k*n+i] * a[k*n+j];
            }
        }
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < i; j++) {
            a[i*n+j] = a[j*n+i];
        }
    }

    return 0; /* success */
}

static void zeros(double * a, int m, int n)
{
    int j;
    for (j=0; j<m*n; ++j)
        a[j] = 0;
}

#ifdef DEBUG
static void dump(double * a, int m, int n, const char * fmt)
{
    int i,j;

    char f[100];
    sprintf(f, "%s ", fmt);
    for(i=0; i<m; ++i) {
        for(j=0; j<n; ++j)
            printf(f, a[i*n+j]);
        printf("\n");
    }
}
#endif

/* C <- A * B */
static void mulmat(double * a, double * b, double * c, int arows, int acols, int bcols)
{
    int i, j,l;

    for(i=0; i<arows; ++i)
        for(j=0; j<bcols; ++j) {
            c[i*bcols+j] = 0;
            for(l=0; l<acols; ++l)
                c[i*bcols+j] += a[i*acols+l] * b[l*bcols+j];
        }
}

static void mulvec(double * a, double * x, double * y, int m, int n)
{
    int i, j;

    for(i=0; i<m; ++i) {
        y[i] = 0;
        for(j=0; j<n; ++j)
            y[i] += x[j] * a[i*n+j];
    }
}

static void transpose(double * a, double * at, int m, int n)
{
    int i,j;

    for(i=0; i<m; ++i)
        for(j=0; j<n; ++j) {
            at[j*m+i] = a[i*n+j];
        }
}

/* A <- A + B */
static void accum(double * a, double * b, int m, int n)
{        
    int i,j;

    for(i=0; i<m; ++i)
        for(j=0; j<n; ++j)
            a[i*n+j] += b[i*n+j];
}

/* C <- A + B */
static void add(double * a, double * b, double * c, int n)
{
    int j;

    for(j=0; j<n; ++j)
        c[j] = a[j] + b[j];
}


/* C <- A - B */
static void sub(double * a, double * b, double * c, int n)
{
    int j;

    for(j=0; j<n; ++j)
        c[j] = a[j] - b[j];
}

static void negate(double * a, int m, int n)
{        
    int i, j;

    for(i=0; i<m; ++i)
        for(j=0; j<n; ++j)
            a[i*n+j] = -a[i*n+j];
}

static void mat_addeye(double * a, int n)
{
    int i;
    for (i=0; i<n; ++i)
        a[i*n+i] += 1;
}

void  eyes_matrix(double *a,int n)
{
	int i=0,j=0;
	for (i=0;i<n;i++)
	{
		for (j=0;j<n;j++)
		{
			if (i==j)
			{
				a[i*n+j]=1;
			} 
			else
			{
				a[i*n+j]=0;
			}
			
		}
	}
}

void meas_noise_set(double *a,double noise,int n)
{
	int i=0,j=0;
	for (i=0;i<n;i++)
	{
		for (j=0;j<n;j++)
		{
			a[i*n+j]*=noise;
		}
	}
}

static void unpack(void * v, kalman_state_t * kalman_state, int n,int m)
{
	 /* skip over n, m in data structure */
	char * cptr = (char *)v;
	cptr += 2*sizeof(int);

	double * dptr = (double *)cptr;
	kalman_state->x = dptr;
	dptr += n;
	kalman_state->P = dptr;
	dptr += n*n;
	kalman_state->Q = dptr;
	dptr += n*n;
	kalman_state->R = dptr;
	dptr += m*m;
	kalman_state->G = dptr;
	dptr += n*m;
	kalman_state->F = dptr;
	dptr += n*n;
	kalman_state->H = dptr;
	dptr += m*n;
	kalman_state->Ht = dptr;
	dptr += n*m;
	kalman_state->Ft = dptr;
	dptr += n*n;
	kalman_state->Pp = dptr;
	dptr += n*n;
	kalman_state->tmp0 = dptr;
	dptr += n*n;
	kalman_state->tmp1 = dptr;
	dptr += n*m;
	kalman_state->tmp2 = dptr;
	dptr += m*n;
	kalman_state->tmp3 = dptr;
	dptr += m*m;
	kalman_state->tmp4 = dptr;
	dptr += m*m;
	kalman_state->tmp5 = dptr;
 }

void kalman_init(void * v, int n,int m)
{
	/* retrieve n and set them in incoming data structure */
	int * ptr = (int *)v;
	*ptr = n;
	ptr++;
	*ptr = n;

	/* unpack rest of incoming structure for initlization */
	kalman_state_t kalman_state;
	unpack(v, &kalman_state, n, m);

	/* zero-out matrices */
	zeros(kalman_state.P, n, n);
	zeros(kalman_state.Q, n, n);
	zeros(kalman_state.R, m, m);
	zeros(kalman_state.G, n, m);
	zeros(kalman_state.F, n, n);
	mat_addeye(kalman_state.F, n);
	zeros(kalman_state.H, m, n);
	mat_addeye(kalman_state.H, n);
}

void update_system_func(kalman_state_t *state, double  time,int n)
{
	state->F[0*n+3] = time;
	state->F[1*n+4] = time;
}

void update_system_noise(kalman_state_t *state, double time,double acc_noise,int n)
{
	state->Q[0*n+0] = acc_noise*time*time*time*time/n;
	state->Q[0*n+2] = 2.0*acc_noise*time*time*time/n;
	state->Q[1*n+1] = acc_noise*time*time*time*time/n;
	state->Q[1*n+3] = acc_noise*time*time*time/n;
	state->Q[2*n+0] = acc_noise*time*time*time/n;
	state->Q[2*n+2] = acc_noise*time*time;
	state->Q[3*n+1] = acc_noise*time*time*time/n;
	state->Q[3*n+3] = acc_noise*time*time;
}

int kalman_step(void * v, double * z)
{        
    /* unpack incoming structure */
	
    int * ptr = (int *)v;
    int n = *ptr;
    ptr++;
    int m = *ptr;
	
    kalman_state_t kalman_state;
    unpack(v, &kalman_state, n, m); 
    update_system_func(&kalman_state,1.0,n);
    update_system_noise(&kalman_state , 1.0,0.1, n);
    
    /* P_k = F_{k-1} P_{k-1} F^T_{k-1} + Q_{k-1} */
    mulmat(kalman_state.F, kalman_state.P, kalman_state.tmp0, n, n, n);
    transpose(kalman_state.F, kalman_state.Ft, n, n);
    mulmat(kalman_state.tmp0, kalman_state.Ft, kalman_state.Pp, n, n, n);
    accum(kalman_state.Pp, kalman_state.Q, n, n);
	
    /* G_k = P_k H^T_k (H_k P_k H^T_k + R)^{-1} */
    transpose(kalman_state.H, kalman_state.Ht, m, n);
    mulmat(kalman_state.Pp, kalman_state.Ht, kalman_state.tmp1, n, n, m);
    mulmat(kalman_state.H, kalman_state.Pp, kalman_state.tmp2, m, n, n);
    mulmat(kalman_state.tmp2, kalman_state.Ht, kalman_state.tmp3, m, n, m);
    accum(kalman_state.tmp3, kalman_state.R, m, m);
    if (cholsl(kalman_state.tmp3, kalman_state.tmp4, kalman_state.tmp5, m)) return 1;
    mulmat(kalman_state.tmp1, kalman_state.tmp4, kalman_state.G, n, m, m);
	
    /* \hat{x}_k = \hat{x_k} + G_k(z_k - h(\hat{x}_k)) */
    sub(z, kalman_state.H, kalman_state.tmp5, m);
    mulvec(kalman_state.G, kalman_state.tmp5, kalman_state.tmp2, n, m);
    add(kalman_state.F, kalman_state.tmp2, kalman_state.x, n);
	
    /* P_k = (I - G_k H_k) P_k */
    mulmat(kalman_state.G, kalman_state.H, kalman_state.tmp0, n, m, n);
    negate(kalman_state.tmp0, n, n);
    mat_addeye(kalman_state.tmp0, n);
    mulmat(kalman_state.tmp0, kalman_state.Pp, kalman_state.P, n, n, n);
	
    /* success */
    return 0;
}
