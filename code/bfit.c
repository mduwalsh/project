#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

#define EPS 1

unsigned long B[] = {};
unsigned long G[] = {128, 64, 32, 16, 8, 4, 2, 1};
unsigned long R[] = {3};


// plot data using lines with gnuplot
void plot(FILE *p, int wid, char *datafile, int columns, char *title, char *xlabel, char *ylabel, int logscale, char *out )
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
  fprintf(p, "set xlabel '{/Hershey/Complex_Italic l}' \n");
  fprintf(p, "set ylabel '{/Helvetica_Oblique e^b}' \n");
  fprintf(p, "set yrange [0:] \n");
  fprintf(p, "f(x) = a/(x-x0) + c \n");
  fprintf(p, "fit f(x) '%s' using 1:2 via a,c,x0\n", datafile);
  fprintf(p, "plot for [col=2:%d] '%s' using 1:col with points pt 3 ps 2 title '', f(x) lt 3 lw 2 title '' \n", columns, datafile);  
}

// plot data using lines with gnuplot
void plotall(FILE *p, int wid, char *datafile, int columns, char *title, char *xlabel, char *ylabel, int logscale, char *out )
{
  int i;
  fprintf(p, "set key outside \n");     
#if EPS
  fprintf(p, "set term postscript eps enhanced color solid font \"Helvetica,25\" \n");
#else
  fprintf(p, "set term pngcairo size 1024, 720 enhanced color solid font \"Helvetica,12\" \n");
#endif
  fprintf(p, "set output '%s' \n", out);
  if(logscale) fprintf(p, "set logscale xy \n");
  fprintf(p, "set title '%s' \n", "");
  fprintf(p, "set rmargin 12 \n");
  fprintf(p, "set bmargin 5 \n");
  fprintf(p, "set xlabel '{/Hershey/Complex_Italic l}' \n");
  fprintf(p, "set ylabel '{/Helvetica_Oblique b}' rotate by 0\n");
  fprintf(p, "set yrange [0:] \n");
  fprintf(p, "set ytics 0.2\n");
  for(i = 0; i < sizeof(G)/sizeof(unsigned long); i++){
    sprintf(datafile, "eb%ld.dat", G[i]);
    fprintf(p, "f%d(x) = a%d/(x-x0%d) + c%d \n", i+1, i+1, i+1, i+1);
    fprintf(p, "fit f%d(x) '%s' using 1:3 via a%d,c%d,x0%d\n", i+1, datafile, i+1, i+1, i+1);
  }
  /*for(i = 0; i < sizeof(G)/sizeof(unsigned long); i++){
    sprintf(datafile, "eb%ld.dat", G[i]);
    if(i == 0) fprintf(p, "plot ");
    fprintf(p, "'%s' using 1:2 with points pt 6 ps 2 lc %d title '', f%d(x) lw 3 lc %d title '%lu'", datafile, i+1, i+1, i+1, G[i]);
    if(i < (sizeof(G)/sizeof(unsigned long) -1) )
      fprintf(p, ", ");
  }
  fprintf(p, "\n");*/
  
  for(i = 0; i < sizeof(G)/sizeof(unsigned long); i++){
    sprintf(datafile, "eb%ld.dat", G[i]);
    if(i == 0) fprintf(p, "plot ");
    fprintf(p, "'%s' using 1:3 with points pt 6 ps 2 lc %d title '', f%d(x) lw 3 lc %d title '%lu'", datafile, i+1, i+1, i+1, G[i]);
    if(i < (sizeof(G)/sizeof(unsigned long) -1) )
      fprintf(p, ", ");
  }
  fprintf(p, "\n");
}

// plot data using lines with gnuplot
void plotslope(FILE *p, int wid, char *datafile, int columns, char *title, char *xlabel, char *ylabel, int logscale, char *out )
{
  int i;
  fprintf(p, "set key outside \n");     
#if EPS
  fprintf(p, "set term postscript eps enhanced color solid font \"Helvetica,25\" \n");
#else
  fprintf(p, "set term pngcairo size 1024, 720 enhanced color solid font \"Helvetica,12\" \n");
#endif
  fprintf(p, "set output '%s' \n", out);
  if(logscale) fprintf(p, "set logscale xy \n");
  fprintf(p, "set title '%s' \n", "");
  fprintf(p, "set rmargin 12 \n");
  fprintf(p, "set bmargin 5 \n");
  fprintf(p, "set xlabel '{/Hershey/Complex_Italic l}' \n");
  fprintf(p, "set ylabel '{/Helvetica_Oblique m}' rotate by 0 \n");
  fprintf(p, "set yrange [-1:0] \n");
  for(i = 0; i < sizeof(G)/sizeof(unsigned long); i++){
    sprintf(datafile, "eb%ld.dat", G[i]);
    fprintf(p, "f%d(x) = a%d*x + c%d \n", i+1, i+1, i+1);
    fprintf(p, "fit f%d(x) '%s' using 1:4 via a%d,c%d\n", i+1, datafile, i+1, i+1);
  }
  // plot slopes
  for(i = 0; i < sizeof(G)/sizeof(unsigned long); i++){
    sprintf(datafile, "eb%ld.dat", G[i]);
    if(i == 0) fprintf(p, "plot ");
    fprintf(p, "'%s' using 1:4 with points pt 6 ps 2 lc %d title '', f%d(x) lw 3 lc %d title '%lu'", datafile, i+1, i+1, i+1, G[i]);
    //fprintf(p, "'%s' using 1:4 with points pt 6 ps 1 lc %d title ''", datafile, i+1);
    if(i < (sizeof(G)/sizeof(unsigned long) -1) )
      fprintf(p, ", ");
  }
  fprintf(p, "\n");
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
   int i;
   char ofname[200], datafile[200];
   FILE *gp = popen ("gnuplot -persistent", "w"); // open gnuplot in persistent mode
   /*for(i = 0; i < sizeof(G)/sizeof(unsigned long); i++){
     sprintf(ofname,"eb%ld.eps", G[i]);  
     sprintf(datafile, "eb%ld.dat", G[i]);
     plot(gp, 0, datafile, 2, "", "", "", 0, ofname);
   }*/
   sprintf(ofname, "b.eps");
   plotall(gp, 0, datafile, 2, "", "", "", 0, ofname);
   
  //sprintf(ofname, "m.eps");
  //plotslope(gp, 0, datafile, 2, "", "", "", 0, ofname);
  fflush(gp);
  pclose(gp);
  return 0;
}

