#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<math.h>

/* sets of parameters to be simulated on */
unsigned long L[] = {2, 4, 6, 8, 10, 12, 14}; // bits
unsigned long G[] = {1000};
double  Epsilon[] = {0.01, 0.05, 0.1, 0.2, 0.5};
unsigned long Seed = 0;
unsigned long Runs = 16;

/* end of sets of parameters to be simulated on */

void create_script(int i, unsigned long b, unsigned long g, double epsilon)
{
  char sfile[20];
  //char path[256];
  FILE *fp;
  sprintf(sfile, "script.%d.sh", i);
  fp = fopen(sfile, "w");
  fprintf(fp, "#PBS -l nodes=1\n");
  fprintf(fp, "#PBS -N vio.b%lug%lue%.2f\n", b, g, epsilon);
  fprintf(fp, "#PBS -l walltime=24:00:00\n" );
  fprintf(fp, "cd $PBS_O_WORKDIR\n");
  //getcwd(path, 256);
  fprintf(fp, "./vio %lu %lu %lu %lf %lu\n", b, Seed, g, epsilon, Runs);
  fclose(fp);  
}

int main()
{
  unsigned long li, gi, ei;
  unsigned long ls, gs, es;
  char scall[200];
  ls  = sizeof(L)/sizeof(unsigned long);
  //ns  = sizeof(N)/sizeof(unsigned long);
  gs  = sizeof(G)/sizeof(unsigned long);
  es  = sizeof(Epsilon)/sizeof(double);
  
  unsigned long i = 1;
  for(li = 0; li < ls; li++)
    //for(ni = 0; ni < ns; ni++)	
      for(gi = 0; gi < gs; gi++)
        for(ei = 0; ei < es; ei++, i++){
	  create_script(i, L[li], G[gi], Epsilon[ei]);	   
	  sprintf(scall, "qsub script.%lu.sh", i);
	  if(system(scall) == -1) printf("\nError executing!!\n");           
	}
  printf("\n%lu\n", i);
  return 0;
}

/** Usage:
    compile : gcc -o viojob viojob.c -lm
    run     : ./viojob
**/



