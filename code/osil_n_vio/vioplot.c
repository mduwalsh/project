#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<math.h>

#define EPS 1

/* sets of parameters to be simulated on */
unsigned long L[] = {2, 4, 6, 8, 10, 12, 14}; // bits
unsigned long G[] = {100};
unsigned long N0 = 8;
unsigned long N;
unsigned long Ni[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 40, 100};  // 2^n
double  Epsilon[] = {0.01, 0.02, 0.05, 0.1, 0.2, 0.5};
unsigned long Seed = 0;

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
  fprintf(p, "plot for [col=2:%d] '%s' using 1:col with lines title '' \n", columns, datafile);  
  
  fflush(p);
  pclose(p);
}

void plot_data_all_infinite(unsigned long b, unsigned long g, double epsilon)
{
  char str[200], ifn[300], ofn[300], title[200];
  //create file name string
  sprintf(str, "b%lug%lueps%.6lf", b, g, epsilon);
  // haploid infinite
  sprintf(ifn, "%s_osc_inf_haploid.dat", str);
  sprintf(ofn, "b%lug%lueps%.2lf_osc_inf_hap.eps", b, g, epsilon);
  sprintf(title, "infinite haploid l:%lu, g:%lu, eps:%.2f", b, g, epsilon);
  plot(ifn, 3, title, "g", "d", 0, ofn);
  // diploid infinite
  sprintf(ifn, "%s_osc_inf_diploid.dat", str);
  sprintf(ofn, "b%lug%lueps%.2lf_osc_inf_dip.eps", b, g, epsilon);
  sprintf(title, "infinite diploid l:%lu, g:%lu, eps:%.2f", b, g, epsilon);
  plot(ifn, 3, title, "g", "d", 0, ofn);  
}

void plot_data_all_finite(unsigned long b, unsigned long g, unsigned long n, double epsilon)
{
  char str[200], ifn[300], ofn[300], title[200];
  //create file name string
  sprintf(str, "b%lug%lun%lueps%.6lf", b, g, n, epsilon);
  
  // haploid finite
  sprintf(ifn, "%s_osc_haploid.dat", str);
  sprintf(ofn, "b%lug%lun%lueps%.2lf_osc_fin_hap.eps", b, g, n, epsilon);
  sprintf(title, "finite haploid l:%lu, g:%lu, n:%lu, eps:%.2f", b, g, n, epsilon);
  plot(ifn, 3, title, "g", "d", 0, ofn);
  // diploid finite
  sprintf(ifn, "%s_osc_diploid.dat", str);
  sprintf(ofn, "b%lug%lun%lueps%.2lf_osc_fin_dip.eps", b, g, n, epsilon);
  sprintf(title, "finite diploid l:%lu, g:%lu, n:%lu, eps:%.2f", b, g, n, epsilon);
  plot(ifn, 3, title, "g", "d", 0, ofn);
}

int main()
{
  unsigned long li, ni, gi, ei;
  unsigned long ls, gs, es;
  ls  = sizeof(L)/sizeof(unsigned long);
  gs  = sizeof(G)/sizeof(unsigned long);
  es  = sizeof(Epsilon)/sizeof(double);
  
  unsigned long i = 1;
  for(li = 0; li < ls; li++)    
      for(gi = 0; gi < gs; gi++)
        for(ei = 0; ei < es; ei++, i++){
	  plot_data_all_infinite(L[li], G[gi], Epsilon[ei]);
	  for(ni = 0; ni < sizeof(Ni)/sizeof(unsigned long); ni++){                        // through all sizes of finite population
	    N = N0*Ni[ni];
	    plot_data_all_finite(L[li], G[gi], N, Epsilon[ei] );            
	  }           
      }
  printf("\n%lu\n", i);
  return 0;
}

/** Usage:
    compile : gcc -o vioplot vioplot.c -lm
    run     : ./vioplot
**/



