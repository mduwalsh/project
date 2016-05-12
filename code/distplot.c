#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

#define EPS 1

unsigned long B[] = {1, 2, 4, 8, 12, 16};                 // haploid bit length 
unsigned long G[] = {1, 2, 4};                            // generations
unsigned long R[] = {1, 3};                               // no. of runs

// plot data using lines with gnuplot
void plot(FILE *p, int wid, char *datafile, int columns, char *title, char *xlabel, char *ylabel, int logscale, char *out, int b, int g )
/*
 * p: file pointer to gnuplot
 * wid: 
 */
{
  fprintf(p, "set key outside \n");     
#if EPS
  fprintf(p, "set term postscript eps enhanced color solid font \"Helvetica,25\" \n");
#else
  fprintf(p, "set term pngcairo size 1024, 720 enhanced color solid font \"Helvetica,12\" \n");
#endif
  fprintf(p, "set output '%s' \n", out);
  if(logscale) fprintf(p, "set logscale xy \n");
  fprintf(p, "set title '%s' \n", "");
  fprintf(p, "set rmargin 2 \n");
  fprintf(p, "set bmargin 5 \n");
  fprintf(p, "set xlabel 'log{/Helvetica-Oblique N}' \n");
  fprintf(p, "set ylabel 'log{/Helvetica-Oblique d}' \n");
  fprintf(p, "f(x) = m*x + b \n");
  fprintf(p, "fit f(x) '%s' using 1:2 via m,b\n", datafile);
  fprintf(p, "plot for [col=2:%d] '%s' using 1:col with points pt 3 ps 2 title '', f(x) lt 3 lw 2 title '' \n", columns, datafile);  
  
}

double log2(double v)
{
  return log(v)/log(2);
}

void convert_data_to_logscale(char *ifname, char *ofname)
{
  int ival, c;
  double dval;
  FILE *fp1, *fp2;
  if(!(fp1 = fopen(ifname,"r"))) {
    printf("\n %s couldn't be opened for read!! \n", ifname);
    return;
  }
  fp2 = fopen(ofname, "w");
  while(c = fscanf(fp1, "%d %lf", &ival, &dval) != EOF){
    fprintf(fp2, "%.6lf  %.6lf \n", log(ival), log(dval));
  }
  fclose(fp1); fclose(fp2);
}

int main()
{
  int i, j, k, nb, ng, nr;
  char datafile[200], title[200], ofname[200];

  nb = sizeof(B)/sizeof(unsigned long);
  ng = sizeof(G)/sizeof(unsigned long);
  nr = sizeof(R)/sizeof(unsigned long);

  for( i = 0; i < nb; i++){
    for( j = 0; j < ng; j++){
      for( k = 0; k < nr; k++){
        sprintf(datafile, "b%lug%lur%d_dist.dat", B[i], G[j], R[k]);
	sprintf(ofname, "b%lug%lur%d_dist_log.dat", B[i], G[j], R[k]);
	convert_data_to_logscale(datafile, ofname);
      }
    }
  }
  
  FILE *gp = popen ("gnuplot -persistent", "w"); // open gnuplot in persistent mode
  for( i = 0; i < nb; i++){
    for( j = 0; j < ng; j++){
      for( k = 0; k < nr; k++){
        sprintf(datafile, "b%lug%lur%lu_dist_log.dat", B[i], G[j], R[k]);
#if EPS
	sprintf(ofname, "b%02lug%04lur%02lu_dist.eps", B[i], G[j], R[k]);
#else
	sprintf(ofname, "b%02lug%04lur%02lu_dist.png", B[i], G[j], R[k]);
#endif
	sprintf(title, "b = %lu, g = %lu, r = %lu", B[i], G[j], R[k]);
        sprintf(title, "l = %lu", B[i]);
	plot(gp, 0, datafile, 2, title, "log(N)", "log(d)", 0, ofname, B[i], G[j]);
      }
    }
  }
  fflush(gp);
  pclose(gp);
  return 0;
}

