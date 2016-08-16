#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<math.h>

#define EPS 1

/* sets of parameters to be simulated on */
unsigned long L[] = {2, 4, 6, 8, 10, 12, 14}; // bits
unsigned long G[] = {1000};
unsigned long GS = 50;                        // number of generations to show in graphs
unsigned long N0 = 8;
unsigned long N;
unsigned long Ni[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 40, 100, 200};  // 2^n
double  Epsilon[] = {0.01, 0.02, 0.05, 0.1, 0.2, 0.5};
unsigned long Seed = 0;
unsigned long Runs = 16;

/* end of sets of parameters to be simulated on */

// plot data using lines with gnuplot
void plot(char *datafile, int columns, char *title, char *xlabel, char *ylabel, int logscale, char *out)
{
  FILE *p = popen("gnuplot -persistent", "w"); // open gnuplot in persistent mode

  fprintf(p, "set key outside \n");     
#if EPS
  fprintf(p, "set term postscript eps enhanced color solid font \"Helvetica,25\" \n");
#else
  fprintf(p, "set term pngcairo size 1024, 720 enhanced color solid font \"Helvetica,12\" \n");
#endif
  fprintf(p, "set output '%s' \n", out);
  if(logscale) fprintf(p, "set logscale xy \n");
  fprintf(p, "set title '%s' \n", title);
  //fprintf(p, "set rmargin 2 \n");
  //fprintf(p, "set bmargin 5 \n");
  fprintf(p, "set xlabel '{/Helvetica-Oblique %s}' \n", xlabel);
  fprintf(p, "set ylabel '{/Helvetica-Oblique %s}' \n", ylabel);
  fprintf(p, "plot for [col=2:%d] '< head -%lu %s' using 1:col with lines title '' \n", columns, GS, datafile);  
  
  fflush(p);
  pclose(p);
}

void plot_data_all_infinite(unsigned long b, unsigned long g, double epsilon, unsigned long run)
{
  char str[200], ifn[300], ofn[300], title[200];
  //create file name string
  sprintf(str, "b%02lu_g%04lu_eps%.6lf", b, g, epsilon);
  // haploid infinite
  sprintf(ifn, "%s_osc_inf_haploid_%02lu.dat", str, run);
  sprintf(ofn, "b%02lu_g%04lu_eps%.6lf_osc_inf_hap_%02lu.eps", b, GS, epsilon, run);
  sprintf(title, "infinite haploid l:%lu, g:%lu, eps:%.2f", b, GS, epsilon);
  plot(ifn, 3, title, "g", "d", 0, ofn);
  // diploid infinite
  sprintf(ifn, "%s_osc_inf_diploid_%02lu.dat", str, run);
  sprintf(ofn, "b%02lu_g%04lu_eps%.6lf_osc_inf_dip_%02lu.eps", b, GS, epsilon, run);
  sprintf(title, "infinite diploid l:%lu, g:%lu, eps:%.2f", b, GS, epsilon);
  plot(ifn, 3, title, "g", "d", 0, ofn);  
}

void plot_data_all_finite(unsigned long b, unsigned long g, unsigned long n, double epsilon, unsigned long run)
{
  char str[200], ifn[300], ofn[300], title[200];
  //create file name string
  sprintf(str, "b%02lu_g%04lu_n%06lu_eps%.6lf", b, g, n, epsilon);
  
  // haploid finite
  sprintf(ifn, "%s_osc_haploid_%02lu.dat", str, run);
  sprintf(ofn, "b%02lu_g%04lu_n%6lu_eps%.6lf_osc_fin_hap_%02lu.eps", b, GS, n, epsilon, run);
  sprintf(title, "finite haploid l:%lu, g:%lu, n:%lu, eps:%.2f", b, GS, n, epsilon);
  plot(ifn, 3, title, "g", "d", 0, ofn);
  // diploid finite
  sprintf(ifn, "%s_osc_diploid_%02lu.dat", str, run);
  sprintf(ofn, "b%02lu_g%04lu_n%06lu_eps%.6lf_osc_fin_dip_%02lu.eps", b, GS, n, epsilon, run);
  sprintf(title, "finite diploid l:%lu, g:%lu, n:%lu, eps:%.2f", b, GS, n, epsilon);
  plot(ifn, 3, title, "g", "d", 0, ofn);
}

int main()
{
  unsigned long li, ni, gi, ei, j;
  unsigned long ls, gs, es;
  ls  = sizeof(L)/sizeof(unsigned long);
  gs  = sizeof(G)/sizeof(unsigned long);
  es  = sizeof(Epsilon)/sizeof(double);
  
  unsigned long i = 1;
  for(li = 0; li < ls; li++)    
      for(gi = 0; gi < gs; gi++)
        for(ei = 0; ei < es; ei++, i++){
          for(j = 0; j < Runs; j++){
	  plot_data_all_infinite(L[li], G[gi], Epsilon[ei], j);
	  for(ni = 0; ni < sizeof(Ni)/sizeof(unsigned long); ni++){                        // through all sizes of finite population
	    N = N0*Ni[ni];
	    plot_data_all_finite(L[li], G[gi], N, Epsilon[ei], j);            
	  }     
        }      
      }
  printf("\n%lu\n", i);
  return 0;
}

/** Usage:
    compile : gcc -o vioplot vioplot.c -lm
    run     : ./vioplot
**/



