#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>

/* sets of parameters to be simulated on */
unsigned long L[] = {2, 4, 6, 8, 12}; // bits
unsigned long g[] = {2, 5, 6, 10, 12, 15, 24};
unsigned long Seed[] = {0, 1, 5, 10, 20, 50, 100, 500, 1000, 1234, 2345, 5000, 10000, 50000, 100000};

/* end of sets of parameters to be simulated on */

void create_script(int i, unsigned long b, unsigned long gg, unsigned long seed)
{
  char sfile[20];
  char path[256];
  FILE *fp;
  sprintf(sfile, "script.%d.sh", i);
  fp = fopen(sfile, "w");
  fprintf(fp, "#PBS -l nodes=1\n");
  fprintf(fp, "#PBS -N opt.b%lug%lus%lu\n", b, gg, seed);
  fprintf(fp, "#PBS -l walltime=24:00:00\n" );
  fprintf(fp, "cd $PBS_O_WORKDIR\n");
  //getcwd(path, 256);
  fprintf(fp, "./opt %lu %lu %lu\n", b, seed, gg);
  fclose(fp);  
}

int main()
{
  unsigned long li, si, gi;
  unsigned long ls, ss, gs;
  char scall[200];
  ls  = sizeof(L)/sizeof(unsigned long);
  ss  = sizeof(Seed)/sizeof(unsigned long);
  gs  = sizeof(g)/sizeof(unsigned long);
  
  unsigned long i = 1;
  for(li = 0; li < ls; li++)
    for(si = 0; si < ss; si++)	
      for(gi = 0; gi < gs; gi++, i++){  
            if(g[gi] < (1ul<<L[li])){       
              create_script(i, L[li], g[gi], Seed[si]);	   
              sprintf(scall, "qsub script.%lu.sh", i);
              system(scall);   
            }
      }
  printf("\n%lu\n", i);
  return 0;
}

/** Usage:
    compile : gcc -o optjob optjob.c -lm
    run     : ./optjob
**/



