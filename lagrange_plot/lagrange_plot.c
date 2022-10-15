#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

// define は，小数点をつけないとint型，つけるとdouble型になる
#define Xmin -5.0
#define Xmax 10.0
#define Ymin -5.0
#define Ymax 15.0
#define N 100

typedef struct {
    double x;
    double y;
}point;

#define NUM 2
point s[NUM] =  {
                    {0.0, 3},
                    {0.5, 1},
                };

double lagrange_inter(double, point *, int num);
double lagrange_diff(double, point *, int num);
point *point_diffCal(point *s);
/* for debug */
void plot(point *s, double (*funcptr)(double, point *, int num), int);

int main (void) {
    int i;
    double ans;

    ans = lagrange_diff(1.0, s, NUM);
    // plot(s, lagrange_inter, NUM);

    printf("ans = %lf\n", ans);

    return 0;
}

point* point_diffCal(point *s)
{   
    point *r = malloc(sizeof (point) * NUM-1);

    for (int i=0; i<NUM-1; i++) {
        r[i].x = (s[i].x + s[i+1].x) / 2;
        r[i].y = (s[i+1].y - s[i].y) / (s[i+1].x - s[i].x);       
    }

    return r;
}

double xconb_inter(int i, double x, point *s, int num)
{
    double ans = 1;
    double sum = 1;

    for(int j=0; j<num; j++) {
        if (i==j) continue;
        ans *= (x - s[j].x);
        sum *= (s[i].x - s[j].x);        
    } 
    return ans / sum;
}

double xconb_diff(int i, double x, point *s, int num)
{
    double ans = 0;
    double sum1 = 1;
    double sum2 = 1;

    for (int j=0; j<num; j++) {
        if (j==i) continue;
        sum1 *= (s[i].x - s[j].x);
        for (int k=0; k<num; k++) {
            if (k==i || k==j) continue;
            sum2 *= (x - s[k].x);
        }
        ans += sum2 ;
        sum2 = 1;
    }
    return ans / sum1;
}

/* complete */
double lagrange_inter(double x, point *s, int num)
{
    double ans = 0;
    double ai;
    for (int i=0; i<num; i++) {
        ai = xconb_inter(i, x, s, num);
        ans += s[i].y * ai;
        printf("a[%d]=%lf\n", i, ai);
    }
    return ans;
}

double lagrange_diff(double x, point *s, int num)
{
    double ans = 0;
    double ai;
    for (int i=0; i<num; i++) {
        ai = xconb_diff(i, x, s, num);
        ans += s[i].y * ai;
        printf("a[%d]=%lf\n", i, ai);
    }
    return ans;
}

void plot_init(FILE *gp)
{
    fprintf(gp, "set terminal png\n");
    fprintf(gp, "set output 'output.png'\n");
    fprintf(gp, "set title 'lagrange interpolation'\n");
    fprintf(gp, "set pointsize 2\n");
    fprintf(gp, "set xlabel 'x-axis'\n");
    fprintf(gp, "set ylabel 'f(x)'\n");
	fprintf(gp, "set xrange [%2.1lf:%2.1lf]\n", Xmin, Xmax);
    fprintf(gp, "set yrange [%2.1lf:%2.1lf]\n", Ymin, Ymax);
}

void plot(point *s, double (*funcptr)(double, point *, int), int num)
{
    FILE *fp;
    FILE *fp2;
    FILE *gp;

    gp = popen("gnuplot -persist","w");
    plot_init(gp);

	if ((fp = fopen("plotdata.dat", "w")) == NULL) {
		printf("Can't open datafile");	exit(1);
	}
    if ((fp2 = fopen("plotdata2.dat", "w")) == NULL) {
		printf("Can't open datafile");	exit(1);
	} 

    double dx = (Xmax - Xmin) / (double) N;

    for (int i=0; i<=N; i++) {
		fprintf(fp, "%lf %f\n", Xmin + (int)i * dx, funcptr(Xmin + (int)i * dx, s, num)); 
	}
    fprintf(fp, "\n\n");
    for (int i=0; i<num; i++) {
        fprintf(fp2, "%lf %lf\n", s[i].x, s[i].y);
    }
    fprintf(fp2, "\n\n");
    fclose(fp2);
    fclose(fp);

    fprintf(gp, "plot 'plotdata.dat' using 1:2 axes x1y1 with lines linewidth 1.5 lc 'black' title 'f(x)',");
    fprintf(gp, "'plotdata2.dat' with points pt 2 lc 'black' title 'sample data'\n");
    usleep(200000);
    fflush(gp);

    fprintf(gp,"set output\n");
	fprintf(gp, "exit\n");

    pclose(gp);
}