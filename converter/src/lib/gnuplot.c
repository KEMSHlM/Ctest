#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <sys/select.h>
#include "common.h"
#include "gnuplot.h"
#include "alm.h"

#define SEP 10

FILE *gp; 
static int NQ[DROPLET_NUMBER], NR[DROPLET_NUMBER];
static int initflg = 0;
static double PI;
static double R[DROPLET_NUMBER];
static double rmax[DROPLET_NUMBER];
static double *r[DROPLET_NUMBER], *rm[DROPLET_NUMBER];

struct timeval tv;

void plot_init() {
    if (initflg) return;

    int i, j;

    NQ[0] = NQs;
    NR[0] = Ngs;
    R[0] = Rss0;
    rmax[0] = rsmax;
    r[0] = dvec(NR[0]+1);
    rm[0] = dvec(NR[0]+1);
    for (j=0; j<=NR[0]; j++) {
        r[0][j] = R[0]*pow(rmax[0]/R[0],((double)j)/((double)NR[0]));
        rm[0][j] = R[0]*pow(rmax[0]/R[0],((double)j+0.5)/((double)NR[0]));
    }
    NQ[1] = NQl;
    NR[1] = Ngl;
    R[1] = Rsl0;
    rmax[1] = rlmax;
    r[1] = dvec(NR[1]+1);
    rm[1] = dvec(NR[1]+1);
    for (j=1; j<=NR[1]; j++) {
        r[1][j] = R[1]*pow(rmax[1]/R[1],((double)j)/((double)NR[1]));
        rm[1][j] = R[1]*pow(rmax[1]/R[1],((double)j+1.5)/((double)NR[1]));
    }
    tv.tv_sec = SLEEP_TIME;
    tv.tv_usec = 0;

    PI = 4.0*atan(1.0);

    initflg = 1;    
}

int check_droplet_number(Coordinate mode) {
    switch (mode) {
        case SCALAR_SPHERICAL_S:
        case VECTOR_SPHERICAL_QS:
        case VECTOR_SPHERICAL_RS:
            return 0;
            break;
        case SCALAR_SPHERICAL_L:
        case VECTOR_SPHERICAL_QL:
        case VECTOR_SPHERICAL_RL:
            return 1;
            break;
        default:
            return 2;
            break;
    }
}

void view_select(Viewmode mode, char *view_opt) {
    switch (mode) {
        case VIEW_2D_COLOR:
            fprintf(gp, "set pm3d map\n");
            fprintf(gp, "set size ratio -1\n");
            strcpy(view_opt, "with pm3d");
            break;
        case VIEW_2D_COLOR_LOG:
            fprintf(gp, "set pm3d map\n");
            fprintf(gp, "set logscale zcb\n");
            fprintf(gp, "set size ratio -1\n");
            strcpy(view_opt, "with pm3d");
            break;
        case VIEW_2D_GREY:
            fprintf(gp, "set palette gray\n");
            fprintf(gp, "set pm3d map\n");
            fprintf(gp, "set size ratio -1\n");
            strcpy(view_opt, "with pm3d");
            break;
        case VIEW_3D:
            fprintf(gp, "set pm3d at b\n");
            fprintf(gp, "set view equal xy\n");
            strcpy(view_opt, "with lines");
            break;
        case VIEW_3D_COLOR:
            fprintf(gp, "set view equal xy\n");
            strcpy(view_opt, "with pm3d");
            break;
        case VIEW_3D_COLOR_LOG:
            fprintf(gp, "set logscale zcb\n");
            fprintf(gp, "set view equal xy\n");
            strcpy(view_opt, "with pm3d");
            break;
        case VIEW_3D_GREY:
            fprintf(gp, "set palette gray\n");
            fprintf(gp, "set view equal xy\n");
            strcpy(view_opt, "with pm3d");
            break;
        default:
            exit(1);
            break;
    }
}

void coord_select(Coordinate mode) {
    int k = check_droplet_number(mode);
    switch (mode) {
        case SCALAR_SPHERICAL_S:
        case VECTOR_SPHERICAL_QS:
        case VECTOR_SPHERICAL_RS:
        case SCALAR_SPHERICAL_L:
        case VECTOR_SPHERICAL_QL:
        case VECTOR_SPHERICAL_RL:
            fprintf(gp, "unset xtics\n");
            fprintf(gp, "unset ytics\n");
            fprintf(gp, "set xrange[-%lf:%lf]\n", rmax[k], rmax[k]);
            fprintf(gp, "set yrange[0:%lf]\n", rmax[k]);
            break;
        case SCALAR_CYLINDRICAL:
        case VECTOR_CYLINDRICAL_Z:
        case VECTOR_CYLINDRICAL_X:
            fprintf(gp, "set xlabel 'NZ'\n");
            fprintf(gp, "set ylabel 'NX'\n");
            fprintf(gp, "set xrange [0:%d]\n", NZ);
            fprintf(gp, "set yrange [0:%d]\n", NX);
            break;
        default: 
            printf("Usage: Coord\n");
            printf("\tSCALAR_SHERICAL_S\n");
            printf("\tSCALAR_SHERICAL_L\n");
            printf("\tSCALAR_CYLINDRICAL\n");
            printf("\tVECTOR_SPHERICAL_QS\n");
            printf("\tVECTOR_SPHERICAL_RS\n");
            printf("\tVECTOR_SPHERICAL_QL\n");
            printf("\tVECTOR_SPHERICAL_RL\n");
            printf("\tVECTOR_CYLINDRICAL_Z\n");
            printf("\tVECTOR_CYLINDRICAL_X\n");
            exit(1);
            break;
    }
}

void para_plot(double **p, Coordinate coord, char *view_opt, char *filename) {
    int i, j, k;
    int i_start=0, j_start=0;
    int i_end, j_end;
    double p_min = 1.0e-10, p_max = 1;
    FILE *fp, *fp2;

    k = check_droplet_number(coord);

    if      (coord == SCALAR_SPHERICAL_S || coord == SCALAR_SPHERICAL_L)    { i_end = NQ[k]-1;  j_end = NR[k]+1; }
    else if (coord == SCALAR_CYLINDRICAL)   { i_end = NZ-1;    j_end = NX-1;  }
    else if (coord == VECTOR_SPHERICAL_QS || coord == VECTOR_SPHERICAL_QL)  { i_start = 1, j_start = 1, i_end = NQ[k]-1;  j_end = NR[k]; }
    else if (coord == VECTOR_SPHERICAL_RS || coord == VECTOR_SPHERICAL_RL)  { i_end = NQ[k]-1;  j_end = NR[k]; }
    else if (coord == VECTOR_CYLINDRICAL_Z) { i_end = NZ;      j_end = NX-1;  }
    else if (coord == VECTOR_CYLINDRICAL_X) { i_end = NZ-1;    j_end = NX;    }

    char plotfile[40], outputfile[40]; 

    sprintf(plotfile, "./output/");
    sprintf(outputfile, "./output/");
    strcat(outputfile, filename);
    strcat(plotfile, filename);
    strcat(outputfile, "_output.dat");
    strcat(plotfile, "_plot.dat");

    fp = fopen(plotfile, "w");
    if (!fp) { printf("error!\tcan't open plotdata\n"); exit(1); } 
    fp2 = fopen(outputfile, "w");
    if (!fp2) { printf("error!\tcan't open outputdata\n"); exit(1); } 

    for (i=i_start; i<=i_end; i++) {
        if (i%OUP_INTERVAL != 0) continue;
            if (i%SEP == 0 && i%OUP_INTERVAL == 0) { 
                fprintf(fp2, "     ");
                for (j=0;j<=j_end;j++) {
                    fprintf(fp2,"%9d\t",j);	
                }
                fprintf(fp2,"\n");
                fprintf(fp2, "     ");
                for (j=0;j<j_end;j++) {
                    fprintf(fp2,"================");	
                }
                fprintf(fp2,"\n");
            }
        fprintf(fp2, "%2d|| ", i);
        for (j=j_start; j<=j_end; j++) {
            if (j%OUP_INTERVAL!=0) continue;
            if (coord == SCALAR_SPHERICAL_L || coord == SCALAR_SPHERICAL_S) {
                if (j == 0) {
                    fprintf(fp, "%lf %lf %lf\n", r[k][j]*cos(((double)i+0.5)*PI/NQ[k]), r[k][j]*sin(((double)i+0.5)*PI/NQ[k]), p[i][j]);
                } else {
                    fprintf(fp, "%lf %lf %lf\n", rm[k][j-1]*cos(((double)i+0.5)*PI/NQ[k]), rm[k][j-1]*sin(((double)i+0.5)*PI/NQ[k]), p[i][j]);
                }
            } else if (coord == VECTOR_SPHERICAL_QL || coord == VECTOR_SPHERICAL_QS) {
                fprintf(fp, "%lf %lf %lf\n", rm[k][j-1]*cos(((double)i)*PI/NQ[k]), rm[k][j-1]*sin(((double)i)*PI/NQ[k]), p[i][j]);
            } else if (coord == VECTOR_SPHERICAL_RL || coord == VECTOR_SPHERICAL_RS) {
                fprintf(fp, "%lf %lf %lf\n", r[k][j]*cos(((double)i+0.5)*PI/NQ[k]), r[k][j]*sin(((double)i+0.5)*PI/NQ[k]), p[i][j]);
            } else {
                fprintf(fp, "%d %d %lf\n", i, j, p[i][j]);
            }
            fprintf(fp2, "%lf\t", p[i][j]);
            if (p_max < p[i][j]) p_max = p[i][j];
            if (p_min > p[i][j]) p_min = p[i][j];
        } 
        fprintf(fp, "\n");
        fprintf(fp2, "\n");
    }

    fflush(fp);
    fprintf(gp, "set zrange[%lf:%lf]\n", p_min, p_max);
    fprintf(gp, "set cbrange[%3.0lf:%3.0lf]\n", p_min, p_max);
    // fprintf(gp, "set zrange[%lf:%lf]\n", 1.0e-5, 10000);
    // fprintf(gp, "set cbrange[%3.0lf:%3.0lf]\n", 1.0e-5, 10000);
    fprintf(gp, "splot '%s' %s\n", plotfile , view_opt);
    fclose(fp);
    fclose(fp2);
}

void plot(double **p, Coordinate coord, Viewmode view, char *paraname, int k) {
    if (!initflg) { printf("error! before plot_init()\n"); exit(1); }

    gp = popen("gnuplot -persist","w");
    if (!gp)  { printf("error! failed to pipe gnuplot\n"); exit(1); }

    char filename[20], charf[10] ,chark[10], view_opt[10];
    view_opt[0] = '\0';

    sprintf(chark,"%03d", k);
    strcpy(filename, paraname);
    strcpy(charf, filename);
    strcat(filename, chark);
    strcat(filename, ".png");

    fprintf(gp, "set title '%s (%d)'\n",paraname, k);
    fprintf(gp, "set terminal png\n");
    fprintf(gp, "set output '%s'\n", filename);
    fprintf(gp, "unset key\n");

    view_select(view, view_opt);
    coord_select(coord);
    para_plot(p, coord, view_opt, charf);
    select(0, NULL, NULL, NULL, &tv);
    filename[0] = '\0';
    fprintf(gp, "exit\n");
    fflush(gp);
    pclose(gp);
}

