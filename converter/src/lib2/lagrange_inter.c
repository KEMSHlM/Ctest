#include <stdio.h>
#include <stdlib.h>
#include "lagrange_inter.h"
#include "alm.h"

DONOR donor_make_2dim (double x0, double x1, double val) {
    DONOR self;

    self.pos[0] = x0;
    self.pos[1] = x1;
    self.val = val;

    return self;
}

DONOR *donor_make (int n) {
    DONOR *self;

    self = malloc(sizeof(DONOR)*n);
    if (self == NULL) {
        fprintf(stderr, "error in donor make!\n");
    } 

    return self;
}

DONOR *donor_free (DONOR *self) {
    free(self); 
    return NULL;
}

double conb_inter(int i, int j, double *x, DONOR *donor)
{
    int k;
    double ans = 1;
    double sum = 1;

    for(k=0; k<DONOR_NUM; k++) {
        if (i==k) continue;
        ans *= (x[j] - donor[i].pos[j]);
        sum *= (donor[k].pos[j] - donor[i].pos[j]);        
    } 
    return ans / sum;
}

double Lagrange_inter (double *x, DONOR *donor) {
    int i, j;
    double ans = 0.0;
    double a = 1.0;

    for (i=0; i<DONOR_NUM; i++) {
        for (j=0; j<LG_DIM; j++) {
            a *= conb_inter(i, j, x, donor); 
        }
        ans += a * donor[i].val;
        a = 1.0;
    }

    return ans;
}

void basic_vector (double *A0, DONOR donor) {
    A0[0] = 1;
    A0[1] = donor.pos[0];
    A0[2] = donor.pos[1];
    A0[3] = donor.pos[0]*donor.pos[0];
    A0[4] = donor.pos[0]*donor.pos[1];
    A0[5] = donor.pos[1]*donor.pos[1];
}

double assignment (double *c, double *x) {
    return c[0] + c[1]*x[0] + c[2]*x[1] + c[3]*x[0]*x[0] + c[4]*x[0]*x[1] + c[5]*x[1]*x[1]; 
}

void normal_matrix (double **A, double *y, DONOR *donor) {
    int n, i, j;
    double *A0 = dvec(LSA_BASIC_VECTOR_NUM);

    for (n=0;n<DONOR_NUM;n++) {
        basic_vector(A0, donor[n]);
        for (i=0;i<LSA_BASIC_VECTOR_NUM;i++) {
            for (j=0;j<LSA_BASIC_VECTOR_NUM;j++) {
                A[i][j] += A0[i]*A0[j];
            }
            y[i] += A0[i]*donor[n].val;
        }
    }

    freedvec(A0, LSA_BASIC_VECTOR_NUM);
}

void LUdecomp (double **A) {
    int i, j, k;
    double val;

    for (i=0;i<LSA_BASIC_VECTOR_NUM;i++) {
        for (j=0;j<=i;j++) {
           val = A[i][j];
            for (k=0;k<j;k++) {
                val -= A[i][k]*A[k][j];
            }
            A[i][j] = val;
        }
        for(j=i+1;j<LSA_BASIC_VECTOR_NUM;j++){
			val = A[i][j];
			for(k=0;k<i;++k){
				val -= A[i][k]*A[k][j];
			}
			A[i][j] = val/A[i][i];
		}
    }
}

void LUsolve (double **A, double *x, double *b) {
    int i, j;
    double val;

	// 前進代入(back substitution)
	//   LY=bからYを計算
    for(i=0; i<LSA_BASIC_VECTOR_NUM; i++){
		val = b[i];
		for(j=0; j<i; j++){
			val -= A[i][j]*x[j];
		}
		x[i] = val/A[i][i];
	}

	// 後退代入(back substitution)
	//  UX=YからXを計算
	for(i=LSA_BASIC_VECTOR_NUM-1; i>=0; --i){
		val = x[i];
		for(j=i+1; j<LSA_BASIC_VECTOR_NUM; j++){
			val -= A[i][j]*x[j];
		}
		x[i] = val;
	}
}

double least_squares_approximation (double *x, DONOR *donor) {
    double ans;
    double **A = dmat(LSA_BASIC_VECTOR_NUM, LSA_BASIC_VECTOR_NUM);
    double *c  = dvec(LSA_BASIC_VECTOR_NUM);
    double *y  = dvec(LSA_BASIC_VECTOR_NUM);

    normal_matrix(A, y, donor);
    LUdecomp(A);
    LUsolve(A, c, y);
    ans = assignment(c, x);

    freedmat(A, LSA_BASIC_VECTOR_NUM, LSA_BASIC_VECTOR_NUM);
    freedvec(c, LSA_BASIC_VECTOR_NUM);
    freedvec(y, LSA_BASIC_VECTOR_NUM);

    return ans;
}