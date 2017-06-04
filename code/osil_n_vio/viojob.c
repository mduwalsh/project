#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<math.h>

/* sets of parameters to be simulated on */
unsigned long L[] = {4, 6, 8, 10, 12, 14}; // bits
unsigned long G[] = {100};
unsigned long Seed = 0;
unsigned long Runs = 16;

/* end of sets of parameters to be simulated on */

void create_script(int i, unsigned long b, unsigned long g)
{
  char sfile[20];
  //char path[256];
  FILE *fp;
  sprintf(sfile, "script.%d.sh", i);
  fp = fopen(sfile, "w");
  fprintf(fp, "#PBS -l nodes=1\n");
  fprintf(fp, "#PBS -N vio.b%lug%lu\n", b, g);
  fprintf(fp, "#PBS -l walltime=24:00:00\n" );
  fprintf(fp, "cd $PBS_O_WORKDIR\n");
  //getcwd(path, 256);
  fprintf(fp, "./vio1 %lu %lu %lu %lu\n", b, Seed, g, Runs);
  fclose(fp);  
}

int main()
{
  unsigned long li, gi;
  unsigned long ls, gs;
  char scall[200];
  ls  = sizeof(L)/sizeof(unsigned long);
  //ns  = sizeof(N)/sizeof(unsigned long);
  gs  = sizeof(G)/sizeof(unsigned long);  
  
  unsigned long i = 1;
  for(li = 0; li < ls; li++)
    //for(ni = 0; ni < ns; ni++)	
      for(gi = 0; gi < gs; gi++, i++){
	  create_script(i, L[li], G[gi]);	   
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



 
