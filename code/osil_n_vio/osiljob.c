#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<math.h>

/* sets of parameters to be simulated on */
unsigned long L[] = {2, 4, 6, 8, 10, 12, 14}; // bits
unsigned long G[] = {50, 100, 200, 500};
unsigned long N[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};  // 2^n
unsigned long Seed = 0;

/* end of sets of parameters to be simulated on */

void create_script(int i, unsigned long b, unsigned long g, unsigned long n)
{
  char sfile[20];
  //char path[256];
  FILE *fp;
  sprintf(sfile, "script.%d.sh", i);
  fp = fopen(sfile, "w");
  fprintf(fp, "#PBS -l nodes=1\n");
  fprintf(fp, "#PBS -N osil.b%lug%lun%lu\n", b, g, n);
  fprintf(fp, "#PBS -l walltime=24:00:00\n" );
  fprintf(fp, "cd $PBS_O_WORKDIR\n");
  //getcwd(path, 256);
  fprintf(fp, "./osil %lu %lu %lu %lu\n", b, Seed, g, n);
  fclose(fp);  
}

int main()
{
  unsigned long li, ni, gi;
  unsigned long ls, ns, gs;
  char scall[200];
  ls  = sizeof(L)/sizeof(unsigned long);
  ns  = sizeof(N)/sizeof(unsigned long);
  gs  = sizeof(G)/sizeof(unsigned long);
  
  unsigned long i = 1;
  for(li = 0; li < ls; li++)
    for(ni = 0; ni < ns; ni++)	
      for(gi = 0; gi < gs; gi++, i++){  
                 
      create_script(i, L[li], G[gi], (unsigned long) pow(2, N[ni]));	   
      sprintf(scall, "qsub script.%lu.sh", i);
      if(system(scall) == -1) printf("\nError executing!!\n");   
            
      }
  printf("\n%lu\n", i);
  return 0;
}

/** Usage:
    compile : gcc -o osiljob osiljob.c -lm
    run     : ./osiljob
**/



