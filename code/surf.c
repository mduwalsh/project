#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

#define EPS 1

unsigned long B[] = {4, 6, 8, 10, 12, 14};
unsigned long G[] = {1, 2, 4, 8, 16, 32, 64, 128};
unsigned long R[] = {3};


// plot data using lines with gnuplot
void plot(FILE *p, int wid, char *datafile, int columns, char *title, char *xlabel, char *ylabel, int logscale, char *out )
{
  int i;
  fprintf(p, "set key outside \n");     
#if EPS
  fprintf(p, "set term postscript eps enhanced color solid font \"Helvetica,20\" \n");
#else
  fprintf(p, "set term pngcairo size 1024, 720 enhanced color solid font \"Helvetica,12\" \n");
#endif
  fprintf(p, "set output '%s' \n", out);
  if(logscale) fprintf(p, "set logscale xy \n");
  fprintf(p, "set title '%s' \n", "");
  fprintf(p, "set rmargin 0 \n");
  fprintf(p, "set lmargin 8 \n");
  fprintf(p, "set offset-0, -2, 0, 0\n");
  fprintf(p, "set xlabel 'log{/Helvetica-Oblique N}' rotate parallel offset 0,-1\n");
  fprintf(p, "set zlabel 'log{/Helvetica-Oblique d}'\n");
  fprintf(p, "set ylabel 'log{/Helvetica-Oblique n}' rotate parallel offset 0,-1\n");
  fprintf(p, "set ticslevel 0 \n");
  fprintf(p, "set ztics 2 \n");
  fprintf(p, "set yrange [0:6] \n");
  fprintf(p, "set ytics 1\n");
  fprintf(p, "set pm3d \n");
  fprintf(p, "set hidden3d \n");
  fprintf(p, "set palette rgbformulae 22,13,-31\n");
  //fprintf(p, "set yrange [0:160] \n");
  //fprintf(p, "f(x) = m*x + b \n");
  //fprintf(p, "fit f(x) '%s' using 1:2 via m,b\n", datafile);
  /*fprintf(p, "splot ");
  for(i = 0; i < sizeof(G)/sizeof(unsigned long); i++){
  fprintf(p, "'%s' using 1:3:2 index %d with points pt 3 ps 1 lc %d title '%lu' ", datafile, i, i+1,  G[i]); 
    if( i < sizeof(G)/sizeof(unsigned long)-1) fprintf(p, ", ");
  } 
  fprintf(p, "\n");*/
  fprintf(p, "splot '%s'  notitle with lines\n", datafile);
  //fprintf(p, "splot exp(-x*x)*exp(-y*y)");
}

int main()
{
  int i, j, k, nb, ng, nr, c;
  double dval, ival;
  char datafile[200], title[200], ofname[200];
  FILE *fp1, *fp2;
  

  nb = sizeof(B)/sizeof(unsigned long);
  ng = sizeof(G)/sizeof(unsigned long);
  nr = sizeof(R)/sizeof(unsigned long);

  for( i = 0; i < nb; i++){
    for( k = 0; k < nr; k++){
      sprintf(ofname, "b%lur%lu_dist_log_gen.dat", B[i], R[k]);
      fp2 = fopen(ofname, "w");
      for( j = 0; j < ng; j++){
        sprintf(datafile, "b%lug%lur%lu_dist_log.dat", B[i], G[j], R[k]);
        fp1 = fopen(datafile, "r");
	while(c = fscanf(fp1, "%lf %lf", &ival, &dval) != EOF){
          fprintf(fp2, "%.6lf  %.6lf  %.6lf \n", ival, log(G[j]), dval);
        }	
        fclose(fp1);
        fprintf(fp2, "\n");
      }
      fclose(fp2);
    }
  }
  
  FILE *gp = popen ("gnuplot -persistent", "w"); // open gnuplot in persistent mode
  for( i = 0; i < nb; i++){
    
      for( k = 0; k < nr; k++){
        sprintf(datafile, "b%lur%lu_dist_log_gen.dat", B[i], R[k]);
#if EPS
	sprintf(ofname, "b%lur%lu_dist_surf.eps", B[i], R[k]);
#else
	sprintf(ofname, "b%lur%lu_dist_surf.png", B[i], R[k]);
#endif
        sprintf(title, "l = %lu", B[i]);
	plot(gp, 0, datafile, 2, title, "log(N)", "log(d)", 0, ofname);
        printf("%d  ", k);
      }
    
  }
  fflush(gp);
  pclose(gp);
  return 0;
}

