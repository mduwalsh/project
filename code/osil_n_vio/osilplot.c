#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<math.h>

/* sets of parameters to be simulated on */
unsigned long L[] = {2, 4, 6, 8, 10, 12, 14}; // bits
unsigned long G[] = {50};
unsigned long N0 = 8;
unsigned long N;
unsigned long Ni[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 40, 100};  // 2^n

#define EPS 1

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

void plot_data_all_infinite(unsigned long b, unsigned long g)
{
  char str[200], ifn[300], ofn[300], title[200];
  //create file name string
  sprintf(str, "b%lug%lu", b, g);
  // haploid infinite
  sprintf(ifn, "%s_osc_inf_haploid.dat", str);
  sprintf(ofn, "%s_osc_inf_hap.eps", str);
  sprintf(title, "infinite haploid l:%lu, g:%lu", b, g);
  plot(ifn, 3, title, "g", "d", 0, ofn);
  // diploid infinite
  sprintf(ifn, "%s_osc_inf_diploid.dat", str);
  sprintf(ofn, "%s_osc_inf_dip.eps", str);
  sprintf(title, "infinite diploid l:%lu, g:%lu", b, g);
  plot(ifn, 3, title, "g", "d", 0, ofn);  
}

void plot_data_all_finite(unsigned long b, unsigned long g, unsigned long n)
{
  char str[200], ifn[300], ofn[300], title[200];
  //create file name string
  sprintf(str, "b%lug%lun%lu", b, g, n);
  
  // haploid finite
  sprintf(ifn, "%s_osc_haploid.dat", str);
  sprintf(ofn, "%s_osc_fin_hap.eps", str);
  sprintf(title, "finite haploid l:%lu, g:%lu, n:%lu", b, g, n);
  plot(ifn, 3, title, "g", "d", 0, ofn);
  // diploid finite
  sprintf(ifn, "%s_osc_diploid.dat", str);
  sprintf(ofn, "%s_osc_fin_dip.eps", str);
  sprintf(title, "finite diploid l:%lu, g:%lu, n:%lu", b, g, n);
  plot(ifn, 3, title, "g", "d", 0, ofn);
}


int main()
{
  unsigned long li, ni, gi;
  unsigned long ls, gs;
  ls  = sizeof(L)/sizeof(unsigned long);  
  gs  = sizeof(G)/sizeof(unsigned long);
  
  unsigned long i = 1;
  for(li = 0; li < ls; li++)    
      for(gi = 0; gi < gs; gi++, i++){  
	plot_data_all_infinite(L[li], G[gi]);
	for(ni = 0; ni < sizeof(Ni)/sizeof(unsigned long); ni++){                        // through all sizes of finite population
	  N = N0*Ni[ni];
	  plot_data_all_finite(L[li], G[gi], N );            
	}
      }
  printf("\n%lu\n", i);
  return 0;
}

/** Usage:
    compile : gcc -o osilplot osilplot.c -lm
    run     : ./osilplot
**/

