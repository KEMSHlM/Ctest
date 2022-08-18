#include <stdio.h>
#include <stdlib.h>

#define NUM 5
#define X 2

typedef struct {
    double x;
    double y;
}point;

point s[NUM] = {{1, 2},
                {2, 5},
                {3, 7},
                {5, 8},
                {7, 10}};

double xconb(int, double, point *);
double lagrange_inter(double, point *);
double lagrange_diff(double, point *);

int main (int argc, char *argv[]) {
    int i;
    double m, n;
    double x;

    if(argc != 2) {
        printf("Usage:\n");
        printf("\texample: %s x[-]\n", argv[0]);
        exit(1);
    }
    x = atol(argv[1]);

    m = lagrange_inter(x, s);
    n = lagrange_diff(x, s);
    printf("x=%2.1lf:\n",x);
    printf(" f(x)    = %lf\n",m);
    printf(" diff(f) = %lf\n",n);

    return 0;
}

double xconb_inter(int i, double x, point *s)
{
    double ans = 1;
    double sum = 1;

    for(int j=0; j<NUM; j++) {
        if (i==j) continue;
        ans *= (x - s[j].x);
        sum *= (s[i].x - s[j].x);        
    } 
    return ans / sum;
}

double xconb_diff(int i, double x, point *s)
{
    double ans = 0;
    double sum1 = 1;
    double sum2 = 1;

    for (int j=0; j<NUM; j++) {
        if (j==i) continue;
        sum1 *= (s[i].x - s[j].x);
        for (int k=0; k<NUM; k++) {
            if (k==i || k==j) continue;
            sum2 *= (x - s[k].x);
        }
        ans += sum2 ;
        sum2 = 1;
    }
    return ans / sum1;
}

/* complete */
double lagrange_inter(double x, point *s)
{
    double ans = 0;
    double ai;
    for (int i=0; i<NUM; i++) {
        ai = xconb_inter(i, x, s);
        ans += s[i].y * ai;
    }
    return ans;
}

double lagrange_diff(double x, point *s)
{
    double ans = 0;
    double ai;
    for (int i=0; i<NUM; i++) {
        ai = xconb_diff(i, x, s);
        ans += s[i].y * ai;
    }
    return ans;
}

