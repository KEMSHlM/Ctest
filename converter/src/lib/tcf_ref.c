#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/select.h>
#include "common.h"
#include "lagrange_inter.h"
#include "alm.h"
#include "tcf_ref.h"
#define SLEEP_TIME 1.0

#if VERSION_SELECT
	#define MIN_T 290
	#define MAX_T 790
	#define MIN_PHI 0.0
	#define MAX_PHI 15.0
	#define TNUM 501
	#define PHINUM 501 
#else 
	#define MIN_T 290
	#define MAX_T 760
	#define MIN_PHI 0.0
	#define MAX_PHI 15.0
	#define TNUM 48
	#define PHINUM 34 
#endif

double **tcfref, *x0;
DONOR *donor;
struct timeval tv1;
#if VERSION_SELECT
	double *Tlist;
	double *philist;
#else
	const double Tlist[TNUM] = {
		290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 
			400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 
			500, 510, 520, 530, 540, 550, 560, 570, 580, 590, 
			600, 610, 620, 630, 640, 650, 660, 670, 680, 690, 
			700, 710, 720, 730, 740, 750, 760
	};
	const double philist[PHINUM] = {
		0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 
			0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9, 
			1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0, 
			10.0, 11.0, 12.0, 13.0, 14.0, 15.0
	};
#endif

void tcf_plot() {
	int i, j;
    FILE *fp;
	FILE *gp;
	char *png = "5atm.png";
	char *filename = "./output/5atm_plot.dat";
    tv1.tv_sec = SLEEP_TIME;
    tv1.tv_usec = 0;

    fp = fopen(filename,"w");
	gp = popen("gnuplot -persist","w");
    if (!gp)  { printf("error! failed to pipe gnuplot\n"); exit(1); }

    fprintf(gp, "set title 'tcf'\n");
    fprintf(gp, "set terminal png\n");
    fprintf(gp, "set output '%s'\n", png);
    fprintf(gp, "unset key\n");
	fprintf(gp, "set pm3d map\n");
	fprintf(gp, "set logscale zcb\n");
	fprintf(gp, "set xrange[290:760]\n");
	fprintf(gp, "set yrange[0:15]\n");
	fprintf(gp, "set xlabel 'T'\n");
	fprintf(gp, "set ylabel 'phi'\n");

	for (i=0; i<TNUM; i++) {
		for (j=0; j<PHINUM; j++) {
			fprintf(fp, "%lf %lf %lf\n", Tlist[i], philist[j],  tcfref[i][j]);
		}
		fprintf(fp,"\n");
	}
    fprintf(gp, "set zrange[0.00001:100000]\n");
    fprintf(gp, "set cbrange[0.00001:100000]\n");
    fprintf(gp, "splot '%s' with pm3d\n", filename);

    select(0, NULL, NULL, NULL, &tv1);
    fprintf(gp, "exit\n");
    fclose(fp);
    pclose(gp);
}

void tcf_read (FILE *fp) {
	int i, j;

	fscanf(fp, "TITLE = ""HEAT MASS TRANSFER""\n");
	fscanf(fp, "VARIABLES = ""phi"", ""T"", ""tcf"",""Tcoolf"",""tig""\n");
	fscanf(fp, "ZONE T= ONLY ZONE, I=48, J=33,F=POINT\n");
	for (i=0; i<TNUM; i++) {
		for (j=0; j<PHINUM; j++) {
			if (j==0) {
				tcfref[i][j] = 1.0e5;
			} else {
				fscanf(fp, "%*lf\t%*lf\t%lf\t%*lf\t%*lf\n", &tcfref[i][j]);
			}
		}
	}
}

int tcf_init (double P) {
	FILE *fp;

	tcfref = dmat(TNUM, PHINUM);
	x0 = dvec(DIM);
	donor = donor_make(DONOR_NUM);

	#if VERSION_SELECT
		int i;
		Tlist = dvec(TNUM);
		philist = dvec(PHINUM);
		for (i=0; i<TNUM; i++) {
			Tlist[i] = MIN_T + (double)i;
		}
		for (i=0; i<PHINUM; i++) {
			philist[i] = MIN_PHI + (double)i*0.03;
		}
	#endif

	if (P==1.0e5) {
		fp = fopen("./tcfref/1atm.dat", "r");
	} else if (P==5.0e5) {
		fp = fopen("./tcfref/5atm.dat", "r");
	} else if (P==10.0e5) {
		fp = fopen("./tcfref/10atm.dat", "r");
	} else {
		fprintf(stderr, "Tcf_reference is not available\n");
		return 1;
	}
	if (!fp) {
		fprintf(stderr, "can't open atm.dat\n");
		return 1;
	}

	tcf_read(fp);
	fclose(fp);
	return 0;
}

/* ok */
int binarySearch_list (double p, const double *list, int N) {
	int a1 = 0, a2 = N-1;
	if (list[a1]>p) return a1; 
	if (list[a2]<p) return a2; 

	do {
		if (list[(a1+a2)/2] > p) a2 = (a1+a2)/2;
		else a1 = (a1+a2)/2;
		if (a2-a1 == 1) return a2;
	} while (1);
}

int search_list (double p, const double *list, int N) {
	int i;
	for (i=0; i<N; i++) {
		if (p<list[i]) return i;
	}
	return N-1;
}

double tcf_referance (double T, double phi) {
	if (T<Tlist[0]||phi<=0) return 1.0e5; 

	int ii=0, jj=0;
	double ans;

	#if VERSION_SELECT
		ii = binarySearch_list(T, Tlist, TNUM);
		jj = binarySearch_list(phi, philist, PHINUM);
	#else
		ii = search_list(T, Tlist, TNUM);
		jj = search_list(phi, philist, PHINUM);
	#endif

	if (ii < TNUM/2) {
		if (jj < PHINUM/2) {
			donor[0] = donor_make_2dim(Tlist[ii-1], philist[jj-1], 	tcfref[ii-1][jj-1]);
			donor[1] = donor_make_2dim(Tlist[ii-1], philist[jj], 	tcfref[ii-1][jj]);
			donor[2] = donor_make_2dim(Tlist[ii-1], philist[jj+1], 	tcfref[ii-1][jj+1]);
			donor[3] = donor_make_2dim(Tlist[ii], 	philist[jj-1], 	tcfref[ii][jj-1]);
			donor[4] = donor_make_2dim(Tlist[ii], 	philist[jj], 	tcfref[ii][jj]);
			donor[5] = donor_make_2dim(Tlist[ii+1], philist[jj-1], 	tcfref[ii+1][jj-1]);
		} else {
			donor[0] = donor_make_2dim(Tlist[ii-1], philist[jj-1], 	tcfref[ii-1][jj-1]);
			donor[1] = donor_make_2dim(Tlist[ii-1], philist[jj], 	tcfref[ii-1][jj]);
			donor[2] = donor_make_2dim(Tlist[ii-1], philist[jj+1], 	tcfref[ii-1][jj+1]);
			donor[3] = donor_make_2dim(Tlist[ii], 	philist[jj], 	tcfref[ii][jj]);
			donor[4] = donor_make_2dim(Tlist[ii], 	philist[jj+1], 	tcfref[ii][jj+1]);
			donor[5] = donor_make_2dim(Tlist[ii+1], philist[jj+1], 	tcfref[ii+1][jj+1]);
		}
	} else {
		if (jj < PHINUM/2) {
			donor[0] = donor_make_2dim(Tlist[ii-2], philist[jj-1], 	tcfref[ii-1][jj-1]);
			donor[1] = donor_make_2dim(Tlist[ii-1], philist[jj-1], 	tcfref[ii][jj-1]);
			donor[2] = donor_make_2dim(Tlist[ii-1], philist[jj], 	tcfref[ii][jj]);
			donor[3] = donor_make_2dim(Tlist[ii], 	philist[jj-1], 	tcfref[ii][jj-1]);
			donor[4] = donor_make_2dim(Tlist[ii], 	philist[jj], 	tcfref[ii][jj]);
			donor[5] = donor_make_2dim(Tlist[ii], 	philist[jj+1], 	tcfref[ii][jj+1]);
		} else {
			donor[0] = donor_make_2dim(Tlist[ii-2], philist[jj], 	tcfref[ii-2][jj]);
			donor[1] = donor_make_2dim(Tlist[ii-1], philist[jj-1], 	tcfref[ii-1][jj-1]);
			donor[2] = donor_make_2dim(Tlist[ii-1], philist[jj], 	tcfref[ii-1][jj]);
			donor[3] = donor_make_2dim(Tlist[ii], 	philist[jj-2], 	tcfref[ii][jj-2]);
			donor[4] = donor_make_2dim(Tlist[ii], 	philist[jj-1], 	tcfref[ii][jj-1]);
			donor[5] = donor_make_2dim(Tlist[ii], 	philist[jj], 	tcfref[ii][jj]);
		}
	}

	x0[0] = T;
	x0[1] = phi;
	ans = least_squares_approximation(x0, donor);

	if (ans<0) ans = 1.0e-16; 
	return ans;
}