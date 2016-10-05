#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include <unistd.h>
#include<stdbool.h>
#include "rand.c"
#include "sort.c"

#define ERR 0.00000001 // error for oscillation
#define PREC "7.4"/* printing precision */

unsigned long L;  // no. of bits to represent chromosome
int Runs;         // no. of runs of simulation
unsigned long N;  // no. of population in finite population
unsigned long G;  // no. of generations to simulate
unsigned long *Pop[2]; // population generation 0 & 1
double Chi_n;      // crossover probability
double Mu_n;       // mutation probability
double *Chi;       // crossover probability array
double *Mu;        // mutation probability array
dist *Cr;          // array of crossover probability for ith gene
dist *M;           // array of mutation probability for ith gene
unsigned long *S;  // unique set of diploids in finite population
double *P0, *P1, *P2, *P3;   // nth, n+1th, n+2th and n+3th generation of haploids in infinite population
int X;             // no. of distinct diploids in finite population
double C, W[256];  // pow(2.,L/2.), pow(2.,-L/2.)*pow(-1.,ones(index))  
double *Z;         // Z = \hat{\Mu}
unsigned long U;   // universal mask ((1ul<<L)-1)  
double *GZ, *GW;   
double **Mh;       // Mh = Mhat; mixing matrix in walsh basis

long mem_count = 0;

unsigned char Ones[256];

void *Malloc(size_t size)
{
  void *p = malloc(size);

  if (p) {
    mem_count += (long)size;
    return p;    
  }
  printf("malloc failed\n");
  _exit(2);
}

#define malloc(x) Malloc(x) 

void *Calloc(size_t nmemb, size_t size)
{
  void *p = calloc(nmemb, size);

  if (p) {
    mem_count += (long)(size*nmemb);
    return p;    
  }
  printf("calloc failed\n");
  _exit(3);
}

#define calloc(x,y) Calloc(x,y) 

void initOneCount() // initialize Ones[i] to number of 1 bits in i
{
  int i,n;  
  for (i = 0; i < 256; i++){
    n = i;
    while(n){
      Ones[i] += n&1;                             // add 1 to count if n has lowest bit '1'
      n >>= 1;
    }
  }
}

unsigned long setbitscount(int n)                 // returns the number of '1' bits in binary representation of n
{
  int count = 0;  
  while(n){
    count += Ones[n&255];                         // add count of lowes 8 bits in n
    n >>= 8;                                      // shift n right by 8 bits
  }
  return count;
}


void initpop(unsigned long *p)                    // write initial diploid population for L haploid bits
{
  int i = N;

  while (i--){
    p[i]   = rndm();
    p[i]  |= (((unsigned long)rndm()) << 32);
    p[i] >>= (64-2*L);
  }  
}

// display bits in integer
void disp_bits(FILE *f, unsigned long x, unsigned n)
{
  int l;

  for (l = 1ul<<(n-1); l; l >>= 1){        // print bits in blocks
    fputc(((x&l)? '1': '0'), f);
  }  
}

// get two haploids x0, x1 from diploid x = x0x1
static inline void get_x0x1(unsigned long x, unsigned *x0, unsigned *x1)
{
  unsigned long m;  
  m = (1ul<<L)-1;    // {000...000}L{111...111}L
  *x1 = x&m;         // 2nd part
  *x0 = x>>L;        // 1st part
}


// calculate diploid<x0,x1> population proportion in finite population
double qf_x(unsigned x0, unsigned x1)
{
  int i;
  unsigned long x;
  double sum;
  x = (x0<<L) + x1;                   // shift left by no. of bits in haploid for 1st part and then add 2nd part to get one diploid

  for(sum = 0, i = 0; i < N; i++){    // through all diploids in finite population
    if(Pop[0][i] == x){               // if population diploid matches with x, increase sum
      sum++;
    }
  }
  return sum/(double)N;
}

// calculate diploid<x0, x1> population porportion from haploids
double qi_x(unsigned x0, unsigned x1, double *p)
{
  return p[x0]*p[x1];
}

// calculates weighted count (proportions) 'p' of haploid g in finite population Pop[0]
void calc_px_from_finite_population(double *p)   
/* 
p: haploids
*/
{
  int i = N;
  unsigned x0, x1;
  bzero(p, (1ul<<L)*sizeof(double));  // reset haploids array p  
  while(i--){                         // loop through finite population
    get_x0x1(Pop[0][i], &x0, &x1);    // get haploids x0 and x1 from diploid Pop[0][i]
    p[x0] += .5/N;                    // sum weight of each haploid out of 2N haploids
    p[x1] += .5/N;
  }    

}

// reproduces two offspring diploids from two parent diploids
void reproduce(int a, int b, int j)   
/*
 *a: index of parent 1
 *b: index of parent 2
 *j: index of offspring to be produced
*/       
{
  unsigned x0, x1, y0, y1, z0, z1;   // haploids
  unsigned m1, m2;
  
  get_x0x1(Pop[0][a], &x0, &x1);          // get haploids x0 and x1 of one parent
  get_x0x1(Pop[0][b], &y0, &y1);          // get haploids y0 and y1 of 2nd parent
  m1 = drand(Cr);                         // crossover mask 1 from distribution Cr
  if(rnd(2) == 0)                         // keep one haploid from one parent randomly
    z0 = drand(M)^((x0&m1)|(x1&~m1));     // apply mutation and crossover mask
  else
    z0 = drand(M)^((x0&m1)|(x1&~m1));     // apply mutation and crossover mask
  m2 = drand(Cr);                         // crossover mask 2
  if(rnd(2) == 0)                         // keep one haploid from 2nd parent randomly
    z1 = drand(M)^((y0&m2)|(y1&~m2));     // apply mutation and crossover mask
  else
    z1 = drand(M)^((y0&m2)|(y1&~m2));     // apply mutation and crossover mask
   
  // combine two haploids from two parents to form a diploid
  Pop[1][j] = (z0<<L) + z1;               // new offspring diploid <z0,z1> at index j  
}

static inline double w(unsigned long i, unsigned long j) // pow(-1., ones(i&j))/C   // walsh function
{
  i &= j;                                                // assumes 32 bit maximum
  i ^= (i >> 16);
  i ^= (i >> 8);
  return W[i&255];
}

void walsh_v(double *x, double *y) // y = Wx              // walsh transform apply to vector x to get y.
{
  unsigned long i,j;
  for (i = 1ul<<L; i--;)
    for (y[i] = 0, j = 1ul<<L; j--; y[i] += x[j]*w(i,j));
    
}

void fwt_v(double *x, double *y) // y = Wx
{
  unsigned long i,j,k, t1, t2, m, n, z;
  double a, b;
  z = 1ul<<L;
  for(i = 0; i < z; i++){
      y[i] = x[i];
  }
  for (i = 0; i < L; i++){
    m = z>>i;       // m = pow(2, L-i);
    n = m>>1;             // n = pow(2, L-i)/2;
    for(j = 0; j < 1ul<<i; j++){       // j < pow(2,i);
      for(k = 0; k < n; k++){
	t1 = m*j+k;
	t2 = m*j + n + k;
	a = y[t1];
	b = y[t2];
	y[t1] = a + b;
	y[t2] = a - b;
      }
    }
  }
  
  a = 1./C;
  for(i = 0; i < z; i++){
      y[i] = a*y[i];
  }
}

#define bar(x) (U&~(x)) /* bar x */ 

double *prime_h(double *x, double *y) // y = x' in walsh basis; x is array of proportions in walsh basis; y is next generation of proportions in walsh basis
{
  unsigned long g,i;
  for (g = 1ul<<L; g--; y[g] *= C){                   // through all possible haploids
    for (y[g] = 0, i = 1ul<<L; i--;)                  // through all possible values of i (2^L)
      if (i == (g&i)) y[g] += x[i]*x[i^g]*Mh[i][i^g]; // sum only if i == g&i
  }
  return y;
}

static int ones(unsigned long x) // returns |x| the number of 1 bits 
{
  int a;
  for (a = x&1; x >>= 1; a += x&1);
  return(a);
}


void walsh(double **x) // x = \hat{M}; equation 8 ; // mixing matrix in walsh basis
{
  unsigned long u,v,q,k;  static double *z = NULL;

  if (!z) z = malloc((1ul<<L)*sizeof(double));
  walsh_v(Mu, z);                                  // z = \hat{\M}
  for (u = 1ul<<L; u--;)
    for (v = 1ul<<L; v--;){
      x[u][v] = 0;
      if (u&v) continue;
      q = bar(u^v);
      for (k = 1ul<<L; k--;) if (k == (k&q)) x[u][v] += Chi[k^u] + Chi[k^v];
      x[u][v] *= (1ul<<(L-1))*z[u]*z[v];
    }  
}

double *g_w(int n, double *x, double *y) // y = x^n; evolve n generations from x 
{
  unsigned long i;  static double *p;

  if (!n){                               // no generations, y = clone x
    for (i = 1ul<<L; i--; y[i] = x[i]);
    return y;
  }
  
  walsh_v(x,GW);                          // w = Wx
  while (n--){
    p = prime_h(GW,GZ); GZ = GW; GW = p;  // p = z = next = w' in walsh basis, then swap pointers
    GW[0] = 1./C;
  }
  walsh_v(GW,y);
  return y;                               // y = final (standard basis)
}

void display_p(double *p)                 // displays haploids' proportion p
{
  int i;
  for(i = 0; i < 1ul<<L; i++){
    printf("%.6lf ", p[i]);
  }
  printf("\n");
}

void display_qf()                          // displays diploids' proportion q for finite population
{
  int i;
  unsigned x0, x1;
  printf("qf: ");
  for(i = 0; i < 1ul<<(2*L); i++){
    get_x0x1(i, &x0, &x1);
    printf("%.4lf ", qf_x(x0, x1));
  }
  printf("\n");
}

// returns distance between two generations of infinite population
double dist_p1p2(double *p1, double *p2)        // p1, p2: haploids array for infinite popn
{
  int i;
  double d, tmp;
  unsigned x0, x1;
  for(d = 0, i = 0; i < 1ul<<(2*L); i++){       // through all diploids from 0 to 2^(2L)-1 (all possibilities)
    get_x0x1(i, &x0, &x1);                      // get haploids x0 and x1
    tmp = qi_x(x0, x1, p1) - qi_x(x0, x1, p2);
    d += tmp*tmp;
  }
  return sqrt(d);
}

// finds oscillation points in p0, p1 or p2 and p3 and returns status if oscillating points are found (1) or not (0)
int osc_points(double *p0, double *p1, double *p2, double *p3)
{
  int i = 0, oscillate = 0; 
  double d1, d2;
  void *tp;
  g_w(1, p0, p1);                       // get n+1 th generation
  //display_p(p0); display_p(p1);  
  while(oscillate !=1 && i<2000){       // restrict generations in case of no oscillation
    g_w(1, p1, p2);                     // n+2 th generation
    g_w(1, p2, p3);                     // n+3 th generation
    // check if oscillation is reached
    d1 = dist_p1p2(p0, p2);             // distance between nth and n+2th generation
    d2 = dist_p1p2(p1, p3);             // distance between n+1th and n+3th generation 
    if(d1 < ERR && d2 < ERR){           // if oscillation is reached; ERR = error in distance to be considered equal
      //display_p(p2); display_p(p3);
      oscillate = 1;	
    }
    else{  // move 2 gens forward
      tp = p0; p0 = p2; p2 = tp;
      tp = p1; p1 = p3; p3 = tp; i++; i++;
    }
  }    
  if(!oscillate){
    printf("\n No oscillation after 2000 generations \n");
    return 0;
  }
  return 1;
}

// plot data using lines with gnuplot
void plot(FILE *p, int wid, char *datafile, int columns, char *title, char *xlabel, char *ylabel, bool logscale, char *out )
{
  fprintf(p, "set key outside \n");        
  fprintf(p, "set term x11 %d \n", wid);
  if(logscale) fprintf(p, "set logscale xy \n");
  fprintf(p, "set title '%s' \n", title);
  fprintf(p, "set xlabel '%s' \n", xlabel);
  fprintf(p, "set ylabel '%s' \n", ylabel);
  //fprintf(p, "set yrange [0:] \n");
  fprintf(p, "plot for [col=2:%d] '%s' using 1:col with lines title '' \n", columns, datafile);
  if (*out){
    fprintf(p, "set term png size 1024,720 enhanced font \"Helvetica,20\" \n");
    fprintf(p, "set output '%s' \n", out);
    fprintf(p, "replot \n");
  }
}

// returns distance by new method that should be faster
double dist_n(double *p)     // p: haploids array for infinite popn
{
  // compute square of distance between population in finite diploid population and infinite diploid population 
  // for diploids only in finite population
  // then compute square of distance for diploids not in finite population but in infinite population
  // that is summation of q(x)^2; x-> <x0,x1> not in S
  // square distance = {summation of (qf_x - qi_x)^2} + {summation of p(x1)^2 * p(x2)^2} - {summation of p(x1)^2*p(x2)^2 <x1,x2> in S}

  double d, tmp, qfx;
  int i, c;
  unsigned x0, x1;
  
  d = 0;c = 0;
  for(i = 0; i < N; i++){    
    if(i < (N-1)){
      if(Pop[0][i] == Pop[0][i+1]){     // move through list till diploid matches
	c++;                            // increase count by one and continue moving forward
	continue;
      }      
    }
    c++;                                // increase count by one
    qfx = c/(double)N;                  // proportion of diploid in finite population
    c = 0;                              // reset count of occurrence of diploid
    
    get_x0x1(Pop[0][i], &x0, &x1);
    
    // compute square of distance between population in finite diploid population and infinite diploid population 
    // for diploids only in finite population
    tmp = qfx - qi_x(x0, x1, p);
    d += tmp*tmp;
    // now subtract summation of p(x1)^2*p(x2)^2 <x1,x2> in S
    tmp =  p[x0]*p[x1];
    d -= tmp*tmp;
  }  
  // add to square distance the summation of p(x1)^2 * p(x2)^2 = (summation of p(x)^2)^2  // refer to paper
  for(tmp = 0, i = 0; i < (1ul<<L); i++)    
    tmp += p[i]*p[i];
  d += tmp*tmp;

  return sqrt(d);
}

// generate mutation distributions for that g such that oscillation condition holds true 
void muDist(unsigned long g)
{
  unsigned long j, i;
  double s;
  for(s = 0, j = 0; j < (1ul<<L); j++){            // through all possible haploids
    i = g&j;                                       // assumes 32 bit maximum
    i ^= (i >> 16);
    i ^= (i >> 8);
    if(ones(i&255)&1){                              // if odd, assign random value
      Mu[j] = U01();
      s += Mu[j];
    }    
  }
  for(j = (1ul<<L); j--;){                          // normalize distribution so that sums up to 1
    Mu[j] /= s;
  }
}

// generate crossover distributions for that g such that oscillation condition holds true
void chiDist(unsigned long g)
{
  unsigned long h, i, k;
  double s;
  h = bar(g); 
  for(s = 0, k = 0; k < 1ul<<L; k++){               // through all possible haploids
    for(i = 0; i < 1ul<<L; i++){                    // all i in R
      if( k == (h&i) ){                             // check if k in bar(g)R
	Chi[k^g] = U01();
	Chi[k] = U01();
        s += Chi[k] + Chi[k^g];
	break;                                      // break search if match found in set bar(g)R
      }
    }
  }
  for(k = 0; k < 1ul<<L; k++){
    Chi[k] /= s;
  }
}

// allocates memory
void setup()
{
  int i;
  double sc = 0, sm = 0;
  C = pow(2., L/2.);   // constant
  U = (1ul<<L)-1;      // universe mask
  for (i = 256; i--; W[i] = ((ones(i)&1) ? -1. : 1.)/C);
  // allocate memory and calculate crossover and mutation 
  Chi = calloc((1ul<<L),sizeof(double));
  Mu = calloc((1ul<<L),sizeof(double));
  Cr = allocdist((1ul<<L));
  M = allocdist((1ul<<L));
  initOneCount();
  P0 = malloc((1ul<<L)*sizeof(double));
  P1 = malloc((1ul<<L)*sizeof(double));
  P2 = malloc((1ul<<L)*sizeof(double));
  P3 = malloc((1ul<<L)*sizeof(double));  
  
  chiDist(12); muDist(12);                         // initialize distributions for Mu and Chi with g = 12
  for(i = 0; i < 1ul<<L; i++){                     // copy mutation and crossover distributions to Cr and M distributions
    Cr->p[i] = Chi[i]; sc += Chi[i];
    M->p[i] = Mu[i]; sm += Mu[i];
  }
  initdist(Cr, sc);                                // initialize Cr distribution
  initdist(M, sm);                                 // initialize M distribution
  
  Mh  = malloc((1ul<<L)*sizeof(double *));
  for (i = 1ul<<L; i--;) Mh[i]  = malloc((1ul<<L)*sizeof(double));
  walsh(Mh);                                       // calculate and install values in mixing matrix in walsh basis (Mhat)
  GZ = malloc((1ul<<L)*sizeof(double));
  GW = malloc((1ul<<L)*sizeof(double)); 
   
}

void init()
{  
  // allocate memory and initialize population 
  Pop[0] = malloc(N*sizeof(unsigned long));
  Pop[1] = malloc(N*sizeof(unsigned long));
  initpop(Pop[0]);     
}

void deinit()
{
  free(Pop[0]); free(Pop[1]);
}

void cleanup()
{
  int i;  
  freedist(Cr); free(Chi);
  freedist(M); free(Mu);
  free(P0); free(P1);
  for (i = 1ul<<L; i--;) free(Mh[i]);
  free(Mh);
  free(Z); free(GW); free(GZ);
}


int main(int argc, char** argv)
{  
  if(argc^7){
    printf("Usage: ./osil bits seed N G Mu Chi\n");
    exit(1);
  }
  if ((sizeof(int) < 4)|| (sizeof(long int) < 8)){
    printf("Assumptions concerning bits are invalid: int has %d bits, unsigned long int has %d bits\n",
	   (int)(8*sizeof(int)), (int)(8*sizeof(unsigned long)));
    return 1;
  }
  L = atol(argv[1]);                     // bit length
  if (L > 32){
    printf("Can't have more thqan 32 bits\n");
    return 1;
  }
  int seed = atoi(argv[2]);
  N = (unsigned long)atoi(argv[3]);      // no. of individuals in finite population
  G = (unsigned long)atol(argv[4]);      // no. of generations
  Mu_n = atof(argv[5]);          
  Chi_n = atof(argv[6]);

  int i,j, n, m;
  double *d1, *d2, dst1, dst2;
  unsigned long *tmp_ptr;
  char fname[200];                       // file name
  FILE *fp;                              // data file pointer
  
  d1 = calloc(G, sizeof(double));
  d2 = calloc(G, sizeof(double));
  
  initrand(seed);
  setup();   
         
  init();                             // initializes memory for pop, M and Cr and also installs values for these  
  merge_sort(Pop[0], N);    
  // infinite population simulation   
  calc_px_from_finite_population(P0); // calculates haploids proportion using finite set
  P1 = g_w(G, P0, P1);
  display_p(P1);
  // find oscillating points for infinite population
  if( osc_points(P0, P1, P2, P3) ){     // if oscillation occurs
    {                                   // compute distance between finite population to oscillating points
      int a, b;  
      //finite population simulation
      for(i = 0; i < G; i++){           // evolve through G generations
	for(j = 0; j < N; j++){         // reproduce N offsprings      
	  a = rnd(N); b = rnd(N);
	  reproduce(a, b, j);           // randomly two parents chosen; will be proportional to proportion
	}
	merge_sort(Pop[0], N); 
	// calculate distance to oscillating points and write to file
	dst1 = dist_n(P2);               // distance to 1st oscillating point
	dst2 = dist_n(P3);               // distance to 2nd oscillating point
	d1[i] = dst1;
	d2[i] = dst2;
	// set new generation as parent generation for next generation
	tmp_ptr = Pop[0];
	Pop[0] = Pop[1];
	Pop[1] = tmp_ptr;		    
      }
    }
  }     
  deinit();        
  cleanup();
  
  sprintf(fname, "b%lug%lun%lu_osc.dat", L, G, N);
  if(!(fp = fopen(fname, "w"))){
    printf("%s could not be opened!! Error!!\n", fname);
    exit(2);
  }
   // write distances to file
   for(n = 0; n < G; n++){
    fprintf(fp, "%d  ", n);
    fprintf(fp, "%" PREC "lf  %" PREC "lf\n", d1[n], d2[n]);
   }
  fclose(fp);
  free(d1); free(d2);
  
  FILE *gp = popen ("gnuplot -persistent", "w"); // open gnuplot in persistent mode
  plot(gp, 0, fname, 3, "oscillation", "G", "d", 0, "" );
  fflush(gp);
  pclose(gp);
  
  return 0;
}


/* Compile and Run:
gcc -O2 -march=native -o dist dist.c -lm
./dist bits seed finite_populaiton_size generations mutationRate crossoverRate
*/


