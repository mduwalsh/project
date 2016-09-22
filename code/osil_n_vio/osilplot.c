#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<math.h>


#define EPS 1
#define SPEC_RUN 0                // 1: plot from specific run ; 0: plot for all runs

/* sets of parameters to be simulated on */
unsigned long L[] = { 8, 10, 12, 14}; // bits
unsigned long G[] = {100};
unsigned long GS = 50;                        // number of generations to show in graphs
unsigned long N0 = 64;
unsigned long N;
unsigned long Ni[] = {30, 40, 50};  // 2^n
unsigned long Runs = 1;
unsigned long SpecialRun = 12;           // when SPEC_RUN = 1


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


void plot_data_all_infinite(unsigned long b, unsigned long g, unsigned long run)
{
  char str[200], ifn[300], ofn[300], title[200];
  //create file name string
  sprintf(str, "b%02lu_g%04lu", b, g);
  // haploid infinite
  sprintf(ifn, "%s_osc_inf_haploid_%02lu.dat", str, run);
  sprintf(ofn, "b%02lu_g%04lu_osc_inf_hap_%02lu.eps", b, GS, run);
  sprintf(title, "infinite haploid {/Helvetica-Oblique l}:%lu, g:%lu", b, GS);
  plot(ifn, 3, title, "g", "d", 0, ofn);
  // diploid infinite
  sprintf(ifn, "%s_osc_inf_diploid_%02lu.dat", str, run);
  sprintf(ofn, "b%02lu_g%04lu_osc_inf_dip_%02lu.eps", b, GS, run);
  sprintf(title, "infinite diploid {/Helvetica-Oblique l}:%lu, g:%lu", b, GS);
  plot(ifn, 3, title, "g", "d", 0, ofn);  
}

void plot_data_all_finite(unsigned long b, unsigned long g, unsigned long n, unsigned long run)
{
  char str[200], ifn[300], ofn[300], title[200];
  //create file name string
  sprintf(str, "b%02lu_g%04lu_n%06lu", b, g, n);
  
  // haploid finite
  sprintf(ifn, "%s_osc_haploid_%02lu.dat", str, run);
  sprintf(ofn, "b%02lu_g%04lu_n%06lu_osc_fin_hap_%02lu.eps", b, GS, n, run);
  sprintf(title, "finite haploid {/Helvetica-Oblique l}:%lu, g:%lu, n:%lu", b, GS, n);
  plot(ifn, 3, title, "g", "d", 0, ofn);
  // diploid finite
  sprintf(ifn, "%s_osc_diploid_%02lu.dat", str, run);
  sprintf(ofn, "b%02lu_g%04lu_n%06lu_osc_fin_dip_%02lu.eps", b, GS, n, run);
  sprintf(title, "finite diploid {/Helvetica-Oblique l}:%lu, g:%lu, n:%lu", b, GS, n);
  plot(ifn, 3, title, "g", "d", 0, ofn);
}

void plot_dist(FILE *p, char *datafile, int columns, char *title, char *xlabel, char *ylabel, int logscale, char *out)
{
  fprintf(p, "set border 3\n");
  fprintf(p, "set key outside \n");     
#if EPS
  fprintf(p, "set term postscript eps enhanced color solid font \"Helvetica,25\" \n");
#else
  fprintf(p, "set term pngcairo size 1024, 720 enhanced color solid font \"Helvetica,12\" \n");
#endif
  fprintf(p, "set output '%s' \n", out);
  if(logscale) fprintf(p, "set logscale xy \n");
  fprintf(p, "set title '%s' \n", title);  
  fprintf(p, "set xtics nomirror \n");
  fprintf(p, "set ytics nomirror \n");
  fprintf(p, "set xlabel '{/Helvetica-Oblique %s}' \n", xlabel);
  fprintf(p, "set ylabel '{/Helvetica-Oblique %s}' \n", ylabel);
  fprintf(p, "plot for [col=2:%d] '< head -%lu %s' using 1:col with lines title '' \n", columns, GS+1, datafile);    
  fprintf(p, "unset border\n");
  fprintf(p, "show border\n");
}

void plot_data_all_dist_fin_inf(unsigned long b, unsigned long g, unsigned long run)
{
  char str[200], ifn[300], ofn[300], title[200];
  unsigned long ni, c, n;
  
  FILE *p = popen("gnuplot -persistent", "w"); // open gnuplot in persistent mode
  fprintf(p, "ymax_h = 0.0\n");
  fprintf(p, "ymax_d = 0.0\n");
  for(c = 0, ni = 0; ni < sizeof(Ni)/sizeof(unsigned long); ni++, c++){                        // through all sizes of finite population
      n = pow(N0,2)*Ni[ni];
      //create file name string
      sprintf(str, "b%02lu_g%04lu_n%06lu", b, g, n);
      // haploid finite
      sprintf(ifn, "%s_osc_haploid_dist_%02lu.dat", str, run);
      fprintf(p, "stats '%s' using 2 prefix 'A%ld' nooutput\n", ifn, c);
      
      fprintf(p, "if(A%ld_max > ymax_h) {ymax_h = A%ld_max}\n", c, c);
      // diploid finite
      sprintf(ifn, "%s_osc_diploid_dist_%02lu.dat", str, run);
      fprintf(p, "stats '%s' using 2 prefix 'B%ld' nooutput\n", ifn, c);
      fprintf(p, "if(B%ld_max > ymax_d) {ymax_d = B%ld_max}\n", c, c);
  }
  for(ni = 0; ni < sizeof(Ni)/sizeof(unsigned long); ni++){                        // through all sizes of finite population
    n = pow(N0,2)*Ni[ni];
    //create file name string
    sprintf(str, "b%02lu_g%04lu_n%06lu", b, g, n);
    
    // haploid finite    
    fprintf(p, "set yrange [0:ymax_h] \n");
    sprintf(ifn, "%s_osc_haploid_dist_%02lu.dat", str, run);
    sprintf(ofn, "b%02lu_g%04lu_n%06lu_osc_fin_hap_dist_%02lu.eps", b, GS, n, run);
    sprintf(title, "distance haploid {/Helvetica-Oblique l}:%lu, g:%lu, n:%lu", b, GS, n);
    plot_dist(p, ifn, 2, title, "g", "d", 0, ofn);
    // diploid finite
    fprintf(p, "set yrange [0:ymax_d] \n");
    sprintf(ifn, "%s_osc_diploid_dist_%02lu.dat", str, run);
    sprintf(ofn, "b%02lu_g%04lu_n%06lu_osc_fin_dip_dist_%02lu.eps", b, GS, n, run);
    sprintf(title, "distance diploid {/Helvetica-Oblique l}:%lu, g:%lu, n:%lu", b, GS, n);
    plot_dist(p, ifn, 2, title, "g", "d", 0, ofn);
  }
  fprintf(p, "set autoscale y \n");
  fflush(p);
  pclose(p);
}

int main()
{
  unsigned long li, ni, gi, j;
  unsigned long ls, gs;
  ls  = sizeof(L)/sizeof(unsigned long);  
  gs  = sizeof(G)/sizeof(unsigned long);
  
  unsigned long i = 1;
  for(li = 0; li < ls; li++)    
      for(gi = 0; gi < gs; gi++, i++){  
#if SPEC_RUN == 0
        for(j = 0; j < Runs; j++){
	  plot_data_all_infinite(L[li], G[gi], j);
	  for(ni = 0; ni < sizeof(Ni)/sizeof(unsigned long); ni++){                        // through all sizes of finite population
	    N = pow(N0,2)*Ni[ni];
	    plot_data_all_finite(L[li], G[gi], N, j);   	    
	  }
	  plot_data_all_dist_fin_inf(L[li], G[gi], j);
        }
#else
	plot_data_all_infinite(L[li], G[gi], SpecialRun);
	for(ni = 0; ni < sizeof(Ni)/sizeof(unsigned long); ni++){                        // through all sizes of finite population
	  N = pow(N0,2)*Ni[ni];
	  plot_data_all_finite(L[li], G[gi], N, SpecialRun);    	  
	}
	plot_data_all_dist_fin_inf(L[li], G[gi], SpecialRun);
#endif
      }
  printf("\n%lu\n", i);
  return 0;
}

/** Usage:
    compile : gcc -o osilplot osilplot.c -lm
    run     : ./osilplot
**/

