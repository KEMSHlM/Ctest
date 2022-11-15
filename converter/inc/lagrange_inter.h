#ifndef LG_H
#define LG_H

/* setting */ 
#define DONOR_NUM 6

#define LG_DEBUG 0
#define LG_DIM 2

/* 
*   最小二乗近似
*   基底ベクトルは，1, x, y, x^2, xy, y^2
*/
#define LSA_BASIC_VECTOR_NUM 6
#define LSA_DIM 2

typedef struct {
    double pos[LG_DIM];
    double val;
} DONOR;

DONOR donor_make_2dim (double x0, double x1, double val);
DONOR *donor_make(int n);
DONOR *donor_free(DONOR *self);
double Lagrange_inter (double *x, DONOR *donor);
double least_squares_approximation (double *x, DONOR *donor);

#endif