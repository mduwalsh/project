#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<math.h>

#define EPS 1
#define MuOrChi 1
/* sets of parameters to be simulated on */
double E[] = {0.01, 0.1, 0.5};

// plot data using lines with gnuplot
void plot(FILE *p, int wid, char *datafile, int columns, char *title, char *xlabel, char *ylabel, int logscale, char *out )
{
  int i;
  fprintf(p, "set key autotitle columnheader\n");
  fprintf(p, "set key outside \n");     
#if EPS
  fprintf(p, "set term postscript eps size 10, 5 enhanced color solid font \"Helvetica,15\" \n");
#else
  fprintf(p, "set term pngcairo size 1024, 720 enhanced color solid font \"Helvetica,12\" \n");
#endif
  fprintf(p, "set output '%s' \n", out);
  if(logscale) fprintf(p, "set logscale xy \n");
  
  fprintf(p, "set multiplot layout 2,3 title '%s' font ',20' \n", title);    // set subplots layout
  fprintf(p, "set xlabel '{/Helvetica-Oblique N}'\n");
  fprintf(p, "set ylabel 'd'\n");
  fprintf(p, "set xrange [0:81920]\n");
  fprintf(p, "set yrange [0.002:0.018]\n");
  fprintf(p, "xr = 81920 - 0\n");
  
  fprintf(p, "set xtics 40000\n");
  //haploid
  for(i = 0; i < 3; i++){    
    fprintf(p, "set title '{/Symbol e} = %g (Haploid) ' \n", E[i]);  
    fprintf(p, "plot for [col=2:%d] '%s' index %d  using 1:col with lines lw 2 notitle , '%s' index %d using 1:6 with lines lt 0 lw 4 notitle \n", 5, datafile, i, datafile, i);  
  }
  //diploid
  for(i = 0; i < 3; i++){    
    fprintf(p, "set title '{/Symbol e} = %g (Diploid) ' \n", E[i]);  
    fprintf(p, "plot for [col=2:%d] '%s' index %d  using 1:col with lines lw 2 notitle , '%s' index %d using 1:6 with lines lt 0 lw 4 notitle \n", 5, datafile, 3+i, datafile, 3+i);  
  }
  fprintf(p, "unset multiplot \n");
}


int main(int argc, char **argv)
{	
  unsigned long j;
  unsigned long i = 1;
  
  char fname[200], title[200];
  sprintf(fname, "%s", argv[1]);  
  FILE *gp = popen ("gnuplot -persistent", "w"); // open gnuplot in persistent mode
  sprintf(title, "Distance");
  plot(gp, 0, fname, 2, title, "", "", 0, "dist_mu.eps"); 
     
  fflush(gp);
  pclose(gp);
  
  return 0;
}

/** Usage:
    compile : gcc -o amp amp.c -lm
    run     : ./amp
**/

 
