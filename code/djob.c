// to run on cluster

#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>

/* sets of parameters to be simulated on */
int Runs[]        = {1, 3, 30};
unsigned long B[]     = {2, 4, 6, 8, 12};
unsigned long G[] = {1, 5, 10, 15, 20, 25, 30};
/* end of sets of parameters to be simulated on */

/* other basic parameters */
unsigned Seed  = 10;     
double   Mu    = 0.001;   
double   Sigma = 0.01;     
/* end of other basic parameters */

void create_script(int i, int runs, unsigned long b, unsigned long gen );

int main()
{
  int ri, bi, gi;
  int rs, bs, gs;
  char scall[50];
  rs  = sizeof(Runs)/sizeof(int);
  bs  = sizeof(B)/sizeof(unsigned long);
  gs  = sizeof(G)/sizeof(unsigned long);
  
  int i = 1;
  for(ri = 0; ri < rs; ri++)
    for(bi = 0; bi < bs; bi++)	
    for(gi = 0; gi < gs; gi++, i++){      
      create_script(i, Runs[ri], B[bi], G[gi]);	   
      sprintf(scall, "qsub script.%d.sh", i);
      system(scall);   
    }
  printf("\n%d\n", i);
  return 0;
}

void create_script(int i, int runs, unsigned long b, unsigned long gen)
{
  char sfile[20];
  char path[256];
  FILE *fp;
  sprintf(sfile, "script.%d.sh", i);
  fp = fopen(sfile, "w");
  fprintf(fp, "#PBS -l nodes=1\n");
  fprintf(fp, "#PBS -N dst.r%db%lug%lu\n", runs, b, gen);
  fprintf(fp, "#PBS -l walltime=100:00:00\n" );
  fprintf(fp, "cd $PBS_O_WORKDIR\n");
  //getcwd(path, 256);
  fprintf(fp, "./dist %lu %d %d %lu %lf %lf\n", b, Seed, runs, gen, Mu, Sigma);
  fclose(fp);  
}


/** Usage:
    compile : gcc -o djob djob.c
    run     : ./djob
**/

