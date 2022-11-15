#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <dirent.h>
#include <string.h>
#include "common.h" 
#include "lagrange_inter.h" 
#include "alm.h" 
#include "tcf_ref.h" 
#include "gnuplot.h"

#define P 5.0e5
#define LDSET 5.0

int NODE,ELM,OUP;
double tref;
double ld[2];
double **zc, **zsl, **zss;
double **rc, **rsl, **rss;
double **damc, **damsl, **damss;
double **Tgc, **Tgsl, **Tgss;
double **phic, **phisl, **phiss;
double **speedc, **speedsl, **speedss; 

bool isDatFile (const char *filename) {
    const char *ext = strrchr(filename, '.');
    return strcmp(".dat", ext) == 0;
}

void damfile_def(char *outputfile, int OUP) {
    char num[10];
    sprintf(num, "%06d", OUP);
    strcpy(outputfile, "./out/dam");
    strcat(outputfile, num);
    strcat(outputfile, ".dat");
}

void tcf_alm () {
    zc       = dmat(NZ+1, NX+1);
    zsl      = dmat(NQl+1, Nl+Ngl+1);
    zss      = dmat(NQs+1, Nl+Ngs+1);
    rc       = dmat(NZ+1, NX+1);
    rsl      = dmat(NQl+1, Nl+Ngl+1);
    rss      = dmat(NQs+1, Nl+Ngs+1);
    damc     = dmat(NZ, NX);
    damsl    = dmat(NQl, Ngl+1);
    damss    = dmat(NQs, Ngs+1);
    Tgc      = dmat(NZ, NX);
    Tgsl     = dmat(NQl, Ngl+1);
    Tgss     = dmat(NQs, Ngs+1);
    phic     = dmat(NZ, NX);
    phisl    = dmat(NQl, Ngl+1);
    phiss    = dmat(NQs, Ngs+1);
    speedc   = dmat(NZ, NX);
    speedsl  = dmat(NQl, Ngl+1);
    speedss  = dmat(NQs, Ngs+1);
}

void para_read (const char *readfile) {
    int i, j;
    FILE *fp;
    char str[40];
	char c;

    fp = fopen(readfile, "r"); 
    if (!fp) {
        fprintf(stderr, "error in para_read 1\n");
        exit(1);
    }

    fgets(str, 40, fp);
    fgets(str, 40, fp);
	fscanf(fp, "ZONE N=%d, E=%d, STRANDID=%d, SOLUTIONTIME=%le, DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL\n", &NODE, &ELM, &OUP, &tref);
	while ((c = fgetc(fp)) != '#');
	fgets(str, 40, fp);
    for (i=0;i<=NZ;i++) {
		for(j=0;j<=NX;j++) {
			fscanf(fp,"%le\t", &zc[i][j]);
		}
	}
	fscanf(fp,"%le\t", &ld[0]);
    for (i=0; i<=NQl; i++) {
        for (j=0; j<Nl+Ngl; j++) {
			fscanf(fp,"%le\t", &zsl[i][j]);
        }
    }
	fscanf(fp,"%le\t", &ld[1]);
    for (i=0; i<=NQs; i++) {
        for (j=0; j<Nl+Ngs; j++) {
			fscanf(fp,"%le\t", &zss[i][j]);
        }
    }
	while ((c = fgetc(fp)) != '#');
	fgets(str, 40, fp);
    for (i=0;i<=NZ;i++) {
		for(j=0;j<=NX;j++) {
			fscanf(fp,"%le\t", &rc[i][j]);
		}
	}
	fscanf(fp,"%*le\t");
    for (i=0; i<=NQl; i++) {
        for (j=0; j<Nl+Ngl; j++) {
			fscanf(fp,"%le\t", &rsl[i][j]);
        }
    }
	fscanf(fp,"%*le\t");
    for (i=0; i<=NQs; i++) {
        for (j=0; j<Nl+Ngs; j++) {
			fscanf(fp,"%le\t", &rss[i][j]);
        }
    }
    
	while ((c = fgetc(fp)) != '#');
	fgets(str, 40, fp);
    for (i=0;i<NZ;i++) {
		for(j=0;j<NX;j++) {
			fscanf(fp,"%le\t", &Tgc[i][j]);
		}
	}
    for (i=0; i<NQl; i++) {
        for (j=0; j<Nl+Ngl; j++) {
			if (j<Nl) fscanf(fp,"%*le\t");
			else fscanf(fp,"%le\t", &Tgsl[i][j+1-Nl]);
        }
    }
    for (i=0; i<NQs; i++) {
        for (j=0; j<Nl+Ngs; j++) {
			if (j<Nl) fscanf(fp,"%*le\t");
			else fscanf(fp,"%le\t", &Tgss[i][j+1-Nl]);
        }
    }
	while ((c = fgetc(fp)) != '#');
	fgets(str, 40, fp);
    for (i=0;i<NZ;i++) {
		for(j=0;j<NX;j++) {
			fscanf(fp,"%le\t", &phic[i][j]);
		}
	}
    for (i=0; i<NQl; i++) {
        for (j=0; j<Nl+Ngl; j++) {
			if (j<Nl) fscanf(fp,"%*le\t");
			else fscanf(fp,"%le\t", &phisl[i][j+1-Nl]);
        }
    }
    for (i=0; i<NQs; i++) {
        for (j=0; j<Nl+Ngs; j++) {
			if (j<Nl) fscanf(fp,"%*le\t");
			else fscanf(fp,"%le\t", &phiss[i][j+1-Nl]);
        }
    }
	while ((c = fgetc(fp)) != '#');
	fgets(str, 40, fp);
    for (i=0;i<NZ;i++) {
		for(j=0;j<NX;j++) {
			fscanf(fp,"%le\t", &speedc[i][j]);
		}
	}
    for (i=0; i<NQl; i++) {
        for (j=0; j<Nl+Ngl; j++) {
			if (j<Nl) fscanf(fp,"%*le\t");
			else fscanf(fp,"%le\t", &speedsl[i][j+1-Nl]);
        }
    }
    for (i=0; i<NQs; i++) {
        for (j=0; j<Nl+Ngs; j++) {
			if (j<Nl) fscanf(fp,"%*le\t");
			else fscanf(fp,"%le\t", &speedss[i][j+1-Nl]);
        }
    }
	while ((c = fgetc(fp)) != '#');
	fgets(str, 40, fp);
    for (i=0;i<NZ;i++) {
		for(j=0;j<NX;j++) {
			fscanf(fp,"%le\t", &damc[i][j]);
		}
	}
    for (i=0; i<NQl; i++) {
        for (j=0; j<Nl+Ngl; j++) {
			if (j<Nl) fscanf(fp,"%*le\t");
			else fscanf(fp,"%le\t", &damsl[i][j+1-Nl]);
        }
    }
    for (i=0; i<NQs; i++) {
        for (j=0; j<Nl+Ngs; j++) {
			if (j<Nl) fscanf(fp,"%*le\t");
			else fscanf(fp,"%le\t", &damss[i][j+1-Nl]);
        }
    }
    fclose(fp);
}

void dam_cal() {
    int i, j;

    for (i=0; i<NZ; i++) {
        for (j=0; j<NX; j++) {
            if (!speedc[i][j]) {
				damc[i][j] = 0.0;
			} else {
				damc[i][j] = 2*Rsl0/(speedc[i][j]*tcf_referance(phic[i][j], Tgc[i][j]));
			}
        }
    }
    for (i=0; i<NQl; i++) {
        for (j=1; j<=Ngl; j++) {
            if (!speedsl[i][j]) {
				damsl[i][j] = 0.0;
			} else {
				damsl[i][j] = 2*Rsl0/(speedsl[i][j]*tcf_referance(phisl[i][j], Tgsl[i][j]));
			}
        }
    }
    for (i=0; i<NQs; i++) {
        for (j=1; j<=Ngs; j++) {
            if (!speedss[i][j]) {
				damss[i][j] = 0;
			} else {
				damss[i][j] = 2*Rsl0/(speedss[i][j]*tcf_referance(phiss[i][j], Tgss[i][j]));
			}
        }
    }
}

void dam_write() {
    FILE *fp;
    int i, j, g;
    char outputfile[100];

    damfile_def(outputfile, OUP);
	g = (NX+1)*(NZ+1)+1;

    fp = fopen(outputfile, "w"); 
    if (!fp) {
        fprintf(stderr, "error in dam_write 1\n");
        exit(1);
    }
    fprintf(fp,"TITLE = ""HEAT MASS TRANSFER""\n");
	fprintf(fp,"VARIABLES = ""z"", ""r"", ""dam""\n");
	fprintf(fp,"ZONE N=%d, E=%d, STRANDID=%d, SOLUTIONTIME=%le, DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL\n",NODE,ELM,OUP,tref);
	fprintf(fp,"VARLOCATION=([3]=CELLCENTERED)\n");
	fprintf(fp, "# z value\n");
	for (i = 0; i <= NZ; i++) {
		for (j = 0; j <= NX; j++) {
			fprintf(fp, "%le\t", 2.0*zc[i][j]);
		//	fprintf(fp, "%le\t", zc[i][j]/(2.0*Rsl0));
		}
		fprintf(fp, "\n");
	}
	fprintf(fp,"%le\n", LDSET);
//	fprintf(fp,"%le\n", ld[0]/(2.0*Rsl0));
	/*large*/
	for (i = 0; i <= NQl; i++) {
		for (j=0; j< Ngl+Nl; j++) {
			fprintf(fp, "%le\t", 2.0*zsl[i][j]);
		//	fprintf(fp, "%le\t", zsl[i][j]/(2.0*Rsl0));
		}
		fprintf(fp, "\n");
	}
	fprintf(fp,"%le\n", ld[1]);
//	fprintf(fp,"%le\n", ld[1]/(2.0*Rsl0));
	/*small*/
	for (i = 0; i <= NQs; i++) {
		for (j=0; j< Ngs+Nl; j++) {
			fprintf(fp, "%le\t", 2.0*zss[i][j]);
		//	fprintf(fp, "%le\t", zss[i][j]/(2.0*Rsl0));
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "# r value\n");
	for (i = 0; i <= NZ; i++) {
		for (j = 0; j <= NX; j++) {
			fprintf(fp, "%le\t", 2.0*rc[i][j]);
		//	fprintf(fp, "%le\t", rc[i][j]/(2.0*Rsl0));
		}
		fprintf(fp, "\n");
	}
	fprintf(fp,"0.0\n");
	/*large*/
	for (i = 0; i <= NQl; i++) {
		for (j=0; j< Ngl+Nl; j++) {
			fprintf(fp, "%le\t", 2.0*rsl[i][j]);
		//	fprintf(fp, "%le\t", rsl[i][j]/(2.0*Rsl0));
		}
		fprintf(fp, "\n");
	}
	fprintf(fp,"0.0\n");
	/*small*/
	for (i = 0; i <= NQs; i++) {
		for (j=0; j< Ngs+Nl; j++) {
			fprintf(fp, "%le\t", 2.0*rss[i][j]);
		//	fprintf(fp, "%le\t", rss[i][j]/(2.0*Rsl0));
		}
		fprintf(fp, "\n");
	}

	fprintf(fp, "# dam value\n");
	for (i = 0; i < NZ; i++) {
		for (j = 0; j < NX; j++) {
			fprintf(fp, "%le\t", damc[i][j]);
		}
		fprintf(fp, "\n");
	}
	/*large*/
	for (i = 0; i < NQl; i++) {
		for (j=0; j< Nl; j++) {
			fprintf(fp, "%le\t", 0.0);
		}
		for (j = 1; j <= Ngl; j++) {
			fprintf(fp, "%le\t", damsl[i][j]);
		}
		fprintf(fp, "\n");
	}
	/*small*/
	for (i = 0; i < NQs; i++) {
		for (j=0; j< Nl; j++) {
			fprintf(fp, "%le\t", 0.0);
		}
		for (j = 1; j <= Ngs; j++) {
			fprintf(fp, "%le\t", damss[i][j]);
		}
		fprintf(fp, "\n");
	}

	fprintf(fp,"# Connectivity list\n");
	for(i=0;i<NZ;i++) {
		for(j=0;j<NX;j++) {
			fprintf(fp,"%d\t%d\t%d\t%d\n",i*(NX+1)+j+1,i*(NX+1)+j+2,(i+1)*(NX+1)+j+2,(i+1)*(NX+1)+j+1);
		}
	}
	/*large*/
	for(i=0;i<NQl;i++) {
		for(j=0;j<Nl+Ngl;j++) {
			if (j==0) {
				fprintf(fp,"%d\t%d\t%d\t%d\n",g,g,g+i*(Nl+Ngl)+1,g+(i+1)*(Nl+Ngl)+1);
			} else {
				fprintf(fp,"%d\t%d\t%d\t%d\n",g+i*(Nl+Ngl)+j,g+i*(Nl+Ngl)+j+1,g+(i+1)*(Nl+Ngl)+j+1,g+(i+1)*(Nl+Ngl)+j);
			}
		}
	}
	g+=(NQl+1)*(Nl+Ngl)+1;
	/*small*/
	for(i=0;i<NQs;i++) {
		for(j=0;j<Nl+Ngs;j++) {
			if (j==0) {
				fprintf(fp,"%d\t%d\t%d\t%d\n",g,g,g+i*(Nl+Ngs)+1,g+(i+1)*(Nl+Ngs)+1);
			} else {
				fprintf(fp,"%d\t%d\t%d\t%d\n",g+i*(Nl+Ngs)+j,g+i*(Nl+Ngs)+j+1,g+(i+1)*(Nl+Ngs)+j+1,g+(i+1)*(Nl+Ngs)+j);
			}
		}
	}
    fclose(fp);
}

void output()
{
	int i,j;
	FILE *fp;

	if ((fp=fopen("test.dat","w"))==NULL) {
		fprintf(stderr, "Can't open datafile\n");
		exit(1);
	}
	fprintf(fp,"# T value (cell-centered)\n");
	for(i=0;i<NZ;i++) {
		for(j=0;j<NX;j++) {
			fprintf(fp,"%le\t",Tgc[i][j]);
		}
		fprintf(fp,"\n");
	}
	/*large*/
	for(i=0;i<NQl;i++) {
		for(j=0;j<Nl+Ngl;j++) {
			if (j<Nl) fprintf(fp,"%le\t",0.0);
			else fprintf(fp,"%le\t",Tgsl[i][j+1-Nl]);
		}
		fprintf(fp,"\n");
	}
	/*small*/
	for(i=0;i<NQs;i++) {
		for(j=0;j<Nl+Ngs;j++) {
			if (j<Nl) fprintf(fp,"%le\t",0.0);
			else fprintf(fp,"%le\t",Tgss[i][j+1-Nl]);
		}
		fprintf(fp,"\n");
	}

	fprintf(fp,"# phi value\n");
	for(i=0;i<NZ;i++) {
		for(j=0;j<NX;j++) {
			fprintf(fp,"%le\t",phic[i][j]);
		}
		fprintf(fp,"\n");
	}
	/*large*/
	for(i=0;i<NQl;i++) {
		for(j=0;j<Nl+Ngl;j++) {
			if (j<Nl) fprintf(fp,"%le\t",0.0);
			else fprintf(fp,"%le\t",phisl[i][j]);
		}
		fprintf(fp,"\n");
	}
	/*small*/
	for(i=0;i<NQs;i++) {
		for(j=0;j<Nl+Ngs;j++) {
			if (j<Nl) fprintf(fp,"%le\t",0.0);
			else fprintf(fp,"%le\t",phiss[i][j]);
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"# speed value\n");
	for(i=0;i<NZ;i++) {
		for(j=0;j<NX;j++) {
			fprintf(fp,"%le\t",speedc[i][j]);
		}
		fprintf(fp,"\n");
	}
	/*large*/
	for(i=0;i<NQl;i++) {
		for(j=0;j<Nl+Ngl;j++) {
			if (j<Nl) fprintf(fp,"%le\t",0.0);
			else fprintf(fp,"%le\t",speedsl[i][j]);
		}
		fprintf(fp,"\n");
	}
	/*small*/
	for(i=0;i<NQs;i++) {
		for(j=0;j<Nl+Ngs;j++) {
			if (j<Nl) fprintf(fp,"%le\t",0.0);
			else fprintf(fp,"%le\t",speedss[i][j]);
		}
		fprintf(fp,"\n");
	}

	fprintf(fp, "# dam value\n");
	for (i = 0; i < NZ; i++) {
		for (j = 0; j < NX; j++) {
			fprintf(fp, "%le\t", damc[i][j]);
		}
		fprintf(fp, "\n");
	}
	/*large*/
	for (i = 0; i < NQl; i++) {
		for (j=0; j< Nl; j++) {
			fprintf(fp, "%le\t", 0.0);
		}
		for (j = 1; j <= Ngl; j++) {
			fprintf(fp, "%le\t", damsl[i][j]);
		}
		fprintf(fp, "\n");
	}
	/*small*/
	for (i = 0; i < NQs; i++) {
		for (j=0; j< Nl; j++) {
			fprintf(fp, "%le\t", 0.0);
		}
		for (j = 1; j <= Ngs; j++) {
			fprintf(fp, "%le\t", damss[i][j]);
		}
		fprintf(fp, "\n");
	}

	fclose(fp);
}

void dam_output (const char *filename) {
    para_read(filename);
    dam_write();
}

int main () {
    DIR *dir;
    struct dirent *dp;
    char path[64] = "./read/";
    char filename[60];
	
    tcf_alm();

    dir = opendir(path);    
    for (dp=readdir(dir); dp!=NULL; dp=readdir(dir)) {
        if (isDatFile(dp->d_name)) {
            printf("Now loading -> %s\n", dp->d_name);
			sprintf(filename, "%s%s", path, dp->d_name);
            dam_output(filename);
        }  
    }
    closedir(dir);

    return 0;
}