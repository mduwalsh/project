#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include <time.h>
#include <unistd.h>
#include<stdbool.h>
#include "rand.c"
#include "sort.c"

#define CLUSTER 1
#define EPS 0
#define VIO_MU_CHI 1       // 0: no violation in Mu nor chi; 1: violate Mu; 2: violate Chi

#define ERR 0.00000001 // error for oscillation
#define PREC "7.16"/* printing precision */

unsigned long L;  // no. of bits to represent chromosome
int Runs;         // no. of runs of simulation

unsigned long N0 = 64;
unsigned long Ni[] = {1, 10, 20};

double Err[] = {0.01, 0.1, 0.5};

unsigned long Nh, Nd;   // size of finite population; N = Ni*N0 
unsigned long G;  // no. of generations to simulate
unsigned long *Pop[2]; // population generation 0 & 1
unsigned long *Hpop[2];     // finite haploid generation 0 and 1
double *Chi;       // crossover probability array
double *Mu;        // mutation probability array
double *Chi_0;       // crossover probability array
double *Mu_0;        // mutation probability array
dist *Cr;          // array of crossover probability for ith gene
dist *M;           // array of mutation probability for ith gene
unsigned long *S;  // unique set of diploids in finite population
double *P0, *P1;   // nth, n+1th generation of haploids in infinite population
int X;             // no. of distinct diploids in finite population
double C, W[256];  // pow(2.,L/2.), pow(2.,-L/2.)*pow(-1.,ones(index))  
double *Z;         // Z = \hat{\Mu}
unsigned long U;   // universal mask ((1ul<<L)-1)  
double *GZ, *GW;   
double **Mh;       // Mh = Mhat; mixing matrix in walsh basis

unsigned long Gg;
unsigned long Seed;
double Epsilon;     // error change in Mu or Chi distribution
unsigned long SkipDiploid;                // skip finite diploid computation
unsigned long SkipHaploid;                // skip finite haploid computation
unsigned long SkipInfinite;               // skip infinite computation
unsigned long ND;

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
  unsigned long i = Nd;

  while (i--){
    p[i]   = rndm();
    p[i]  |= (((unsigned long)rndm()) << 32);
    p[i] >>= (64-2*L);
  }  
}


// display bits in integer
void disp_bits(FILE *f, unsigned long x, unsigned n)
{
  unsigned long l;

  for (l = 1ul<<(n-1); l; l >>= 1){        // print bits in blocks
    fputc(((x&l)? '1': '0'), f);
  }  
}

// get two haploids x0, x1 from diploid x = x0x1
static inline void get_x0x1(unsigned long x, unsigned *x0, unsigned *x1)
{
  unsigned long m;  
  m = (1ul<<L)-1;    // {000...000}L{111...111}L
  *x1 =(unsigned) x&m;         // 2nd part
  *x0 = (unsigned)x>>L;        // 1st part
}


// calculate diploid<x0,x1> population proportion in finite population
double qf_x(unsigned x0, unsigned x1)
{
  unsigned long i;
  unsigned long x;
  double sum;
  x = (x0<<L) + x1;                   // shift left by no. of bits in haploid for 1st part and then add 2nd part to get one diploid

  for(sum = 0, i = 0; i < Nd; i++){    // through all diploids in finite population
    if(Pop[0][i] == x){               // if population diploid matches with x, increase sum
      sum++;
    }
  }
  return sum/(double)Nd;
}

// calculate diploid<x0, x1> population porportion from haploids
double qi_x(unsigned x0, unsigned x1, double *p)
{
  //if(x0 >= (1ul<<L) || x1 >= (1ul<<L)) printf("Error!!\n x0: %u x1: %u\n", x0, x1);  
  return p[x0]*p[x1];
}

// calculates weighted count (proportions) 'p' of haploid g in finite population Pop[0]
void calc_px_from_finite_population(double *p)   
/* 
p: haploids
*/
{
  unsigned long i = Nh;
  unsigned x0, x1;
  bzero(p, (1ul<<L)*sizeof(double));  // reset haploids array p  
  while(i--){                         // loop through finite population
    get_x0x1(Pop[0][i], &x0, &x1);    // get haploids x0 and x1 from diploid Pop[0][i]
    p[x0] += .5/Nh;                    // sum weight of each haploid out of 2N haploids
    p[x1] += .5/Nh;
  }    
}

void install_pop_distribution(double *p)
{
  unsigned long i;
  double s;
  for(i = 0, s = 0; i < 1ul<<L; i++){
    p[i] = U01();
    s += p[i];
  }
  for(i = 0; i < 1ul<<L; i++){
    p[i] = p[i]/s;    
  }
}

// generates finite haploid initial population from given distribution p
void generate_fin_hap_init_pop_from_random_px(double *p, unsigned long *pop)
{
  unsigned long i;
  double s;
  dist *pd;
  pd = allocdist(1ul<<L);
  for(s = 0, i = 0; i < 1ul<<L; i++){
    pd->p[i] = p[i];
    s += p[i];
  }
  initdist(pd, s);
  for(i = 0; i < Nh; i++){
    pop[i] = (unsigned long)drand(pd);
  }  
  freedist(pd);
}

// calculates weighted count (proportions) 'p' of haploid g in finite population pop
void calc_px_from_finite_hap_pop(double *p, unsigned long *pop)   
/* 
p: haploid proportions
pop: haploid population
*/
{
  unsigned long i = Nh;
  bzero(p, (1ul<<L)*sizeof(double));  // reset haploids array p  
  while(i--){                         // loop through finite population
    p[pop[i]] += 1.0/Nh;               // sum weight of each haploid out of 2N haploids    
  }    
}

// generate finite haploid population corresponding to population vector p
void generate_fin_hap_pop_from_pvector(double *p, unsigned long *pop)
{
  unsigned long i, j, k;
  for(k = 0, i = 0; i < (1ul<<L); i++){
    for(j = 0; j < p[i]*Nh; j++){
      pop[k++] = i;
    }
  }
}

// install/ generate finite diploid population (size N^2 ) from distribution for haploids (size N)
void generate_fin_dipop_from_px(double *p, unsigned long *pop)
{
  unsigned long i, j, k, l;
  //double sp;
  //dist *pd;
  //pd = allocdist(1ul<<(2*L));    
  //sp = 0;
  for( l = 0, i = 0; i < 1ul<<L; i++){
    for(j = 0; j < 1ul<<L; j++){
      //pd->p[(i<<L)+j] = (p[i]*p[j]*N) - (int)(p[i]*p[j]*N);    // difference in exact expected population to real population
      //printf("p member size: %lf \n ", p[i]*p[j]*Nd);
      for(k = 0; k < ( (unsigned long)(p[i]*p[j]*Nd) ); k++){
        pop[l++] = (i<<L) + j;                          // diploid filling population	
      }
      //sp += pd->p[(i<<L)+j];
    }
  }
  //initdist(pd, sp);
  //while (l < N) pop[l++] = (unsigned long)(drand(pd));  // fill up missing number of population with population supposed to be present in populaiton according to weight 
  //freedist(pd);
}
// reproduces an offspring diploid from two parent diploids
void reproduce_diploid(int a, int b, int j, unsigned long *p0, unsigned long *p1)   
/*
 *a: index of parent 1
 *b: index of parent 2
 *j: index of offspring to be produced
 *p0: finite parent generation popn 
 *p1: next generation popn
*/       
{
  unsigned x0, x1, y0, y1, z0, z1;   // haploids
  unsigned m1, m2;
  
  get_x0x1(p0[a], &x0, &x1);          // get haploids x0 and x1 of one parent
  get_x0x1(p0[b], &y0, &y1);          // get haploids y0 and y1 of 2nd parent
  m1 = drand(Cr);                         // crossover mask 1 from distribution Cr
  if(rnd(2) == 0)                         // keep one haploid from one parent randomly
    z0 = drand(M)^((x0&m1)|(x1&~m1));     // apply mutation and crossover mask
  else
    z0 = drand(M)^((x1&m1)|(x0&~m1));     // apply mutation and crossover mask
  m2 = drand(Cr);                         // crossover mask 2
  if(rnd(2) == 0)                         // keep one haploid from 2nd parent randomly
    z1 = drand(M)^((y0&m2)|(y1&~m2));     // apply mutation and crossover mask
  else
    z1 = drand(M)^((y1&m2)|(y0&~m2));     // apply mutation and crossover mask
   
  // combine two haploids from two parents to form a diploid
  p1[j] = (z0<<L) + z1;               // new offspring diploid <z0,z1> at index j  
}

// reproduces two offspring haploids from two parent haploids
void reproduce_haploid(int a, int b, int j, unsigned long *p0, unsigned long *p1)   
/*
 *a: index of parent 1
 *b: index of parent 2
 *j: index of offspring to be produced
 *p0: finite parent generation popn 
 *p1: next generation popn
*/       
{
  unsigned m = drand(Cr);                         // crossover mask 1 from distribution Cr
  if(rnd(2) == 0)                         // keep one haploid from one parent randomly
    p1[j] = (unsigned long) drand(M)^((p0[a]&m)|(p0[b]&~m));     // apply mutation and crossover mask
  else
    p1[j] =(unsigned long) drand(M)^((p0[b]&m)|(p0[a]&~m));     // apply mutation and crossover mask    
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


void walsh(double **x, int fi) // x = \hat{M}; equation 8 ; // mixing matrix in walsh basis
{
  unsigned long u,v,q,k;  static double *z = NULL;
  if(fi == -1){ free(z); return;}
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
  
  walsh_v(x,GW);                          // GW = Wx
  while (n--){
    p = prime_h(GW,GZ); GZ = GW; GW = p;  // p = z = next = w' in walsh basis, then swap pointers
    GW[0] = 1./C;
  }
  walsh_v(GW,y);
  return y;                               // y = final (standard basis)
}

void display_p(double *p, char* str)                 // displays haploids' proportion p
{
  unsigned long i;
  double s = 0;
  printf("%s\n", str);
  for(i = 0; i < 1ul<<L; i++){
    printf("%e ", p[i]);
    s += p[i];
  }
  printf("\n");
  printf("sum: %lf\n", s);
}

void display_qf()                          // displays diploids' proportion q for finite population
{
  unsigned long i;
  unsigned x0, x1;
  printf("qf: ");
  for(i = 0; i < 1ul<<(2*L); i++){
    get_x0x1(i, &x0, &x1);
    printf("%.4lf ", qf_x(x0, x1));
  }
  printf("\n");
}

void display_qix(double *p)
{
  unsigned long i;
  unsigned x0, x1;
  printf("inf qi_x: \n");
  for(i = 0; i < 1ul<<(2*L); i++){       // through all diploids from 0 to 2^(2L)-1 (all possibilities)
    get_x0x1(i, &x0, &x1);                      // get haploids x0 and x1
    printf("%.6lf ", qi_x(x0, x1, p));
  }
  printf("\n");
}

void display_Mhat()
{
  unsigned long i, j;
  printf("Mh:\n");
  for(i = 0; i < (1ul<<L); i++){
    for(j = 0; j < (1ul<<L); j++){
      printf("%.10lf ", Mh[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}

// returns distance between two generations of infinite diploid population
double dist_p1p2_diploid(double *p1, double *p2)        // p1, p2: haploids array for infinite diploid popn
{
  unsigned long i;
  double d, tmp;
  unsigned x0, x1;
  for(d = 0, i = 0; i < 1ul<<(2*L); i++){       // through all diploids from 0 to 2^(2L)-1 (all possibilities)
    get_x0x1(i, &x0, &x1);                      // get haploids x0 and x1
    tmp = qi_x(x0, x1, p1) - qi_x(x0, x1, p2);
    d += tmp*tmp;
  }
  return sqrt(d);
}

// returns distance between two generations of infinite haploid population
double dist_p1p2_haploid(double *p1, double *p2)        // p1, p2: haploids array for infinite popn
{
  unsigned long i;
  double d, tmp;
  for(d = 0, i = 0; i < (1ul<<L); i++){       // through all diploids from 0 to 2^(2L)-1 (all possibilities)
    tmp = p1[i] - p2[i];
    d += tmp*tmp;
  }
  return sqrt(d);
}

// returns distance from origin to diploid population represented by haploids
double dist_o(double *p)
{
  unsigned long i;
  unsigned x0, x1;
  double tmp, d;
  for(d = 0, i = 0; i < 1ul<<(2*L); i++){       // through all diploids from 0 to 2^(2L)-1 (all possibilities)
    get_x0x1(i, &x0, &x1);                      // get haploids x0 and x1
    tmp = qi_x(x0, x1, p);
    d += tmp*tmp;
  }
  return sqrt(d);
}

// returns distance by new method that should be faster
double dist_n(double *p, unsigned long *pop)     // p: haploids array for infinite popn
{
  // compute square of distance between population in finite diploid population and infinite diploid population 
  // for diploids only in finite population
  // then compute square of distance for diploids not in finite population but in infinite population
  // that is summation of q(x)^2; x-> <x0,x1> not in S
  // square distance = {summation of (qf_x - qi_x)^2} + {summation of p(x1)^2 * p(x2)^2} - {summation of p(x1)^2*p(x2)^2 <x1,x2> in S}

  double d, tmp, qfx;
  unsigned long i, c;
  unsigned x0, x1;
  d = 0;c = 0;
  for(i = 0; i < Nd; i++){        
    if(i < (Nd-1)){
      if(pop[i] == pop[i+1]){           // move through list till diploid matches
	c++;                            // increase count by one and continue moving forward
	continue;
      }      
    }
    c++;                                // increase count by one
    qfx = c/(double)(Nd);                  // proportion of diploid in finite population
    c = 0;                              // reset count of occurrence of diploid
    
    get_x0x1(pop[i], &x0, &x1);
    //if(x0 < 0 || x0 > 3|| x1 < 0 || x1 > 3) printf("x0:%lu x1:%lu pop[%lu]: %lu\n", x0, x1, i, pop[i]);
    
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

double dist_haploid(double *p, unsigned long *pop)     // p: haploids array for infinite popn
{
  double d, tmp, qfx;
  unsigned long i, c;
  
  d = 0;c = 0;
  for(i = 0; i < Nh; i++){    
    if(i < (Nh-1)){
      if(pop[i] == pop[i+1]){     // move through list till diploid matches
	c++;                            // increase count by one and continue moving forward
	continue;
      }      
    }
    c++;                                // increase count by one
    qfx = c/(double)Nh;                  // proportion of diploid in finite population
    c = 0;                              // reset count of occurrence of diploid
    
    // compute square of distance between population in finite haploid population and infinite haploid population 
    // for diploids only in finite population
    tmp = qfx - p[pop[i]];
    d += tmp*tmp;
    // now subtract summation of p(x1)^2 for x1 in finite popn
    d -=  p[pop[i]]*p[pop[i]];    
  }  
  // add to square distance the summation of p(x1)^2 for all x1 in R
  for(tmp = 0, i = 0; i < (1ul<<L); i++)    
    tmp += p[i]*p[i];
  d += tmp;

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
    Mu_0[j] = Mu[j];
  }  
}

void violate_mu()
{
  double s;
  unsigned long j;
  if( VIO_MU_CHI == 1){                             // violation in Mu distribution such that oscillation condition does not hold true
    Mu[0] = Epsilon;                                // set first element to Epsilon
    s = Mu[0];
    for(j = 1; j < (1ul<<L); j++){                  // remove Epsilon proportion from other elements
      Mu[j] *= (1.0 - Epsilon); 
      s += Mu[j];
    }
    for(j = 0; j < (1ul<<L); j++){                  // normalize
      Mu[j] /= s;      
    }
  }
}

// generate crossover distributions for that g such that oscillation condition holds true
void chiDist(unsigned long g)
{
  unsigned long h, k;
  double s;
  h = bar(g); 
  for(s = 0, k = 0; k < (1ul<<L); k++){               // through all possible haploids    
    if( k == (h&k) ){                             // check if k in bar(g)R
      Chi[k^g] = U01();
      Chi[k] = U01();
      s += Chi[k] + Chi[k^g];	
    }    
  }
  for(k = 0; k < (1ul<<L); k++){
    Chi[k] /= s;
    Chi_0[k] = Chi[k];
  }  
}


void violate_chi()
{
  double s;
  unsigned long k, j;
  if( VIO_MU_CHI == 2){                             // violation in Chi distribution such that oscillation condition does not hold true
    s =0;
    j = 0;
    for(k = 0; k < (1ul<<L); k++){                  
      if(!j){
	if(Chi[k] == 0.0){ 
	  Chi[k] = Epsilon;         // install value epsilon to first element in Chi distribution with zero value.
	  j = 1;
	}
      }
      else{
	Chi[k] *= (1.0 - Epsilon); 	            // remove Epsilon proportion from other elements
      }
      s += Chi[k];
    }
    
    for(k = 0; k < (1ul<<L); k++){                  // normalize
      Chi[k] /= s;      
    }
  }
}

void copy_array(double *s, double *d, unsigned long size)
{
  unsigned long i;;
  for(i = 0; i < size; i++)
    d[i] = s[i];
}

// plot data using lines with gnuplot
void plot(FILE *p, int wid, char *datafile, int columns, char *title, char *xlabel, char *ylabel, bool logscale, char *out )
{
  fprintf(p, "set key outside \n");   
  if(EPS){
    fprintf(p, "set term postscript eps enhanced color solid font \"Helvetica,25\" \n");
    fprintf(p, "set output '%s.eps' \n", datafile);
  }
  else{     
    fprintf(p, "set term x11 %d \n", wid);
  }
  if(logscale) fprintf(p, "set logscale xy \n");
  fprintf(p, "set title '%s' \n", title);
  fprintf(p, "set xlabel '%s' \n", xlabel);
  fprintf(p, "set ylabel '%s' \n", ylabel);
  //fprintf(p, "set yrange [0:] \n");
  // fprintf(p, "plot for [col=2:%d] '< tail -50 %s' using 1:col with lines title '' \n", columns, datafile);
  fprintf(p, "plot for [col=2:%d] '%s' using 1:col with lines title '' \n", columns, datafile);
  if (*out){
    fprintf(p, "set term png size 1024,720 enhanced font \"Helvetica,20\" \n");
    fprintf(p, "set output '%s' \n", out);
    fprintf(p, "replot \n");
  }
}

double x_g(unsigned long g)
{
  if(g == 0){
    return 2;
  }
  else{
    return 2*Mh[g][0];
  }
}

double y_g(double *z, unsigned long g)
{
  double s;
  unsigned long i;
  if(ones(g) == 1) return 0;
  if(g == 0) return 0;
  for(s = 0, i = 1; i < (1ul<<L); i++){                 // all possible indices in z such that i in S_g
    if(i == g ) continue;
    if(i == (g&i)){
      if(isnan(z[i]) || isnan(z[i^g])){
	printf("NaN!! Error calculating z[i] or z[i^g]\n");
	exit(1);
      }
      s += z[i]*z[i^g]*Mh[i][i^g];
    }
  }
  return C*s;
}

// sorts si indices according to values in sv
void sort_v(unsigned long *si, unsigned long *sv, unsigned long n)
{
  unsigned long i, j, t;
  // sort si according to values sv
  for(i = 0; i < n; i++){
    for(j = 0; j < n-i-1; j++){
      if(sv[j] > sv[j+1]){
	// swap values
	t = sv[j];
	sv[j] = sv[j+1];
	sv[j+1] = t;
	// swap indices
	t = si[j];
	si[j] = si[j+1];
	si[j+1] = t;
      }
    }
  }  
}

// computes periodic orbits p_str and q_str for population p
void comp_periodic_orbits(double *p_str, double *q_str, double *p)
{
  unsigned long i, *si, *sv;
  double *ph, *qh, *pw, *mpw;
  pw = malloc((1ul<<L)*sizeof(double));
  mpw = malloc((1ul<<L)*sizeof(double));
  ph = malloc((1ul<<L)*sizeof(double));
  qh = malloc((1ul<<L)*sizeof(double));
  si = malloc((1ul<<L)*sizeof(unsigned long));
  sv = malloc((1ul<<L)*sizeof(unsigned long));
  for(i = 0; i < (1ul<<L); i++){                 // initializes values of indices
    si[i] = i;
    sv[i] = ones(i);                       // count of 1's in index i
  }
  sort_v(si, sv, 1ul<<L);                 // sort si indices according to count of 1s

  walsh_v(p, pw);                         // haploid proportions p in walsh basis
  prime_h(pw, mpw);                       // next generation M(p) in walsh basis hat(M(p))
  ph[0] = 1./C;
  qh[0] = 1./C;
  for( i = 1; i < (1ul<<L); i++){         // initialize hat p_str and q_str to NaNs
    ph[i] = NAN;
    qh[i] = NAN;   
  }
  // compute gth components of p_str and q_str in increasing value of |g|
  unsigned long g; double x;
  //display_Mhat();
  //printf("x_g, y_g: ");
  for(i = 1; i < (1ul<<L); i++){
    g = si[i];
    x = x_g(g);
    //printf("%.6lf, %.6lf ", x, y_g(qh, g));
    if( fabs(x) < 1.0){
      ph[g] = ( x*y_g(ph,g) + y_g(qh,g) )/( 1.0 -x*x );
      qh[g] = ( x*y_g(qh,g) + y_g(ph,g) )/( 1.0 -x*x );
    }
    else{
      ph[g] = pw[g];
      qh[g] = mpw[g];
    }
      
  }
  //printf("\n");
  walsh_v(ph, p_str); walsh_v(qh, q_str);  // periodic orbits in standard basis
  
  free(pw); free(mpw); free(ph); free(qh); free(si); free(sv);
}

void prep_filename_hap(char *fname, char *apndstr)
{
  sprintf(fname, "b%02lu_g%04lu_n%08lu_eps%.6lf_%s", L, G, Nh, Epsilon, apndstr);
}

void prep_filename_dip(char *fname, char *apndstr)
{
  sprintf(fname, "b%02lu_g%04lu_n%08lu_eps%.6lf_%s", L, G, Nd, Epsilon, apndstr);
}


void inf_osc(double *p, double *p_str, double *q_str, unsigned long run)      //p: initial population haploids; p_str & q_str: oscillating points
{
  FILE *fp, *gp;
  char fname[200], title[200];
  unsigned long i;
  double *p1, *p2, *tptr;
  double *d1 = calloc(G, sizeof(double));
  double *d2 = calloc(G, sizeof(double));
  double *d3 = calloc(G, sizeof(double));
  double *d4 = calloc(G, sizeof(double));
  p1 = malloc((1ul<<L)*sizeof(double));
  p2 = malloc((1ul<<L)*sizeof(double));
  for(i = 0; i < 1ul<<L; i++){                // clone initial population
    p1[i] = p[i];
  }
  for(i = 0; i < G; i++){                     // move through G generations
    g_w(1, p1, p2);                           // evolution of population 1 generation at a time
    d1[i] = dist_p1p2_haploid(p2, p_str);     // distance between p2 and p_str haploid
    d2[i] = dist_p1p2_haploid(p2, q_str);     // distance between p2 and q_str haploid
    d3[i] = dist_p1p2_diploid(p2, p_str);     // distance between p2 and p_str diploid
    d4[i] = dist_p1p2_diploid(p2, q_str);     // distance between p2 and q_str diploid  
    // swap pointers to p1 and p2
    tptr = p1;
    p1 = p2;
    p2 = tptr;
  }

  // write haploid distances to file
  char str[100];
  sprintf(str,  "b%02lu_g%04lu_eps%0.6lf_osc_inf_haploid_%02lu.dat", L, G, Epsilon, run);             // add run to filename
  sprintf(fname, "%s", str);  // haploid distance data file 
  if(!(fp = fopen(fname, "w"))){
    printf("%s could not be opened!! Error!!\n", fname);
    exit(2);
  }
  for(i = 0; i < G; i++){
    fprintf(fp, "%lu  ", i);
    fprintf(fp, "%e  %e \n", d1[i], d2[i]);
  }
  fclose(fp); 
  // plot 
  if(!CLUSTER){
    sprintf(title, "b%lu g%lu gg%lu s%lu infinite haploid", L, G, Gg, Seed);    
    gp = popen ("gnuplot -persistent", "w"); // open gnuplot in persistent mode
    plot(gp, 0, fname, 3, title, "G", "d", 0, "" );
    fflush(gp);
    pclose(gp);
  }

  // write diploid distances to file
  sprintf(str,  "b%02lu_g%04lu_eps%.6lf_osc_inf_diploid_%02lu.dat", L, G, Epsilon, run);             // add run to filename
  sprintf(fname, "%s", str);  // diploid distance data file 
  if(!(fp = fopen(fname, "w"))){
    printf("%s could not be opened!! Error!!\n", fname);
    exit(2);
  }
  for(i = 0; i < G; i++){
    fprintf(fp, "%lu  ", i);
    fprintf(fp, "%e  %e \n", d3[i], d4[i]);
  }
  fclose(fp);
  free(d1); free(d2); free(d3); free(d4); free(p1); free(p2);
  // plot
  if(!CLUSTER){
    sprintf(title, "b%lu g%lu gg%lu s%lu infinite diploid", L, G, Gg, Seed);
    gp = popen ("gnuplot -persistent", "w"); // open gnuplot in persistent mode
    plot(gp, 0, fname, 3, title, "G", "d", 0, "" );
    fflush(gp);
    pclose(gp);
  }
}

void fin_osc_diploid(double *p_str, double *q_str, unsigned long run) // p_str & q_str oscillating points
{
  int a, b;  
  FILE *fp, *gp;
  char fname[200], title[200];
  unsigned long i, j, *tmp_ptr;
  double *d1 = calloc(G, sizeof(double));
  double *d2 = calloc(G, sizeof(double));
  //finite population simulation
  for(i = 0; i < G; i++){           // evolve through G generations
    for(j = 0; j < Nd; j++){         // reproduce N offsprings      
      a = rnd(Nd); b = rnd(Nd);
      reproduce_diploid(a, b, j, Pop[0], Pop[1]);           // randomly two parents chosen; will be proportional to proportion
    }
    // set new generation as parent generation for next generation
    tmp_ptr = Pop[0];
    Pop[0] = Pop[1];
    Pop[1] = tmp_ptr;		
    merge_sort(Pop[0], Nd); 
    // calculate distance to oscillating points and write to file
    d1[i] = dist_n(p_str, Pop[0]);               // distance to 1st oscillating point
    d2[i] = dist_n(q_str, Pop[0]);               // distance to 2nd oscillating point	    
  }
  
  // write distances
  char str[100];
  sprintf(str, "osc_diploid_%02lu.dat", run);             // add run to filename
  prep_filename_dip(fname, str);
  if(!(fp = fopen(fname, "w"))){
    printf("%s could not be opened!! Error!!\n", fname);
    exit(2);
  }
  // write distances to file
  for(i = 0; i < G; i++){
    fprintf(fp, "%lu  ", i);
    fprintf(fp, "%e  %e\n", d1[i], d2[i]);
  }
  fclose(fp);
  free(d1); free(d2);
  
  if(!CLUSTER){
    sprintf(title, "b%lu g%lu n%lu gg%lu s%lu finite diploid", L, G, Nd, Gg, Seed);
    gp = popen ("gnuplot -persistent", "w"); // open gnuplot in persistent mode
    plot(gp, 0, fname, 3, "finite diploid", "G", "d", 0, "" );
    fflush(gp);
    pclose(gp);
  }
}

void fin_osc_haploid(double *p_str, double *q_str, unsigned long run) // p_str & q_str oscillating points
{
  int a, b;  
  FILE *fp, *gp;
  char fname[200], title[200];
  unsigned long i, j, *tmp_ptr;
  double *d1 = calloc(G, sizeof(double));
  double *d2 = calloc(G, sizeof(double));
  //finite population simulation
  for(i = 0; i < G; i++){           // evolve through G generations
    for(j = 0; j < Nh; j++){         // reproduce N offsprings      
      a = rnd(Nh); b = rnd(Nh);
      reproduce_haploid(a, b, j, Hpop[0], Hpop[1]);           // randomly two parents chosen; will be proportional to proportion
    }
    // set new generation as parent generation for next generation
    tmp_ptr = Hpop[0];
    Hpop[0] = Hpop[1];
    Hpop[1] = tmp_ptr;		
    merge_sort(Hpop[0], Nh); 
    // calculate distance to oscillating points and write to file
    d1[i] = dist_haploid(p_str, Hpop[0]);               // distance to 1st oscillating point
    d2[i] = dist_haploid(q_str, Hpop[0]);               // distance to 2nd oscillating point	    
  }

  // write distances
  char str[100];
  sprintf(str, "osc_haploid_%02lu.dat", run);
  prep_filename_hap(fname, str);
  if(!(fp = fopen(fname, "w"))){
    printf("%s could not be opened!! Error!!\n", fname);
    exit(2);
  }
   // write distances to file
   for(i = 0; i < G; i++){
    fprintf(fp, "%lu  ", i);
    fprintf(fp, "%elf  %elf\n", d1[i], d2[i]);
   }
  fclose(fp);
  free(d1); free(d2);
  
  if(!CLUSTER){
    sprintf(title, "b%lu g%lu n%lu gg%lu s%lu finite haploid", L, G, Nh, Gg, Seed);
    gp = popen ("gnuplot -persistent", "w"); // open gnuplot in persistent mode
    plot(gp, 0, fname, 3, title, "G", "d", 0, "" );
    fflush(gp);
    pclose(gp);
  }
}



double dot_product(double *x, double *y, unsigned long size)
{
  double s = 0;
  unsigned long i;
  for(i = 0; i < size; i++){
    s += x[i]*y[i];
  }
  return s;
}

void vector_sub(double *a, double *b, double *c, unsigned long size)
{
  unsigned long i;
  for(i = 0; i < size; i++){
    c[i] = a[i] - b[i];
  }
}

void compute_unit_vector(double *x, double *y, double *n, unsigned long size )
{
  unsigned long i;
  double dn;
  dn = dist_p1p2_haploid(x, y);
  for(i = 0; i < size; i++){
    n[i] = (x[i]-y[i])/dn;
  }
}

double dot_product_diploid_n_and_pop_sub_p(double *p, double *n, unsigned long *pop)     // p: haploids array for infinite popn
{
  // compute square of distance between population in finite diploid population and infinite diploid population 
  // for diploids only in finite population
  // then compute square of distance for diploids not in finite population but in infinite population
  // that is summation of q(x)^2; x-> <x0,x1> not in S
  // square distance = {summation of (qf_x - qi_x)^2} + {summation of p(x1)^2 * p(x2)^2} - {summation of p(x1)^2*p(x2)^2 <x1,x2> in S}

  double d, tmp, qfx;
  unsigned long i, c;
  unsigned x0, x1;
  d = 0;c = 0;
  for(i = 0; i < Nd; i++){        
    if(i < (Nd-1)){
      if(pop[i] == pop[i+1]){           // move through list till diploid matches
	c++;                            // increase count by one and continue moving forward
	continue;
      }      
    }
    c++;                                // increase count by one
    qfx = c/(double)(Nd);                  // proportion of diploid in finite population
    c = 0;                              // reset count of occurrence of diploid
    
    get_x0x1(pop[i], &x0, &x1);    
    // compute square of distance between population in finite diploid population and infinite diploid population 
    // for diploids only in finite population
    tmp = qfx - qi_x(x0, x1, p);    
    d += qi_x(x0, x1, n)*tmp;
    // now add summation of n(x1)*n(x2)*p(x1)*p(x2) <x1,x2> in S
    d += qi_x(x0, x1, n)*qi_x(x0, x1, p);
  }  
  
  // add to square distance the summation of p(x1)^2 * p(x2)^2 = (summation of p(x)^2)^2  // refer to paper
  for(tmp = 0, i = 0; i < (1ul<<L); i++)    
    tmp += p[i]*n[i];
  d -= tmp*tmp;

  return d;
}

double dot_product_haploid_n_and_pop_sub_p(double *p, double *n, unsigned long *pop)     // p: haploids array for infinite popn
{
  double d, tmp, qfx;
  unsigned long i, c;
  
  d = 0;c = 0;
  for(i = 0; i < Nh; i++){    
    if(i < (Nh-1)){
      if(pop[i] == pop[i+1]){     // move through list till diploid matches
	c++;                            // increase count by one and continue moving forward
	continue;
      }      
    }
    c++;                                // increase count by one
    qfx = c/(double)Nh;                  // proportion of diploid in finite population
    c = 0;                              // reset count of occurrence of diploid
    
    // compute square of distance between population in finite haploid population and infinite haploid population 
    // for diploids only in finite population
    tmp = qfx - p[pop[i]];
    d += n[pop[i]]*tmp;
    // now add summation of n[x1]p(x1) for x1 in finite popn
     d +=  n[pop[i]]*p[pop[i]];    
  }  
  // add to square distance the summation of p(x1)^2 for all x1 in R
  for(tmp = 0, i = 0; i < (1ul<<L); i++)    
    tmp += n[i]*p[i];
  d -= tmp;

  return d;
}

void osc_all_n_dist(double *p, double *p_str, double *q_str, double *p1_str, double *q1_str, unsigned long run, double *n_str)
/*
 * p: initial population haploids;
 * p_str: oscillation point 1
 * q_str: oscillation point 2
 */
{
  FILE *fp, *gp, *fp2, *fp3;
  char fname[200], title[200], tfile[200];
  unsigned long i, j, a, b, *tmp_ptr;
  double *p1, *p2, *tptr;
  double *d1, *d2, *d3, *d4, *d5, *d6, *d7, *d8, *d9, *d10, *d11, *d12, *d13, *d14, *d15, *d16, *d17, *d18;
  d1 = calloc(G, sizeof(double));
  d2 = calloc(G, sizeof(double));
  if(!SkipDiploid){
    d3 = calloc(G, sizeof(double));
    d4 = calloc(G, sizeof(double));
  }
  d5 = calloc(G, sizeof(double));
  d6 = calloc(G, sizeof(double));
  if(!SkipDiploid){
    d7 = calloc(G, sizeof(double));
    d8 = calloc(G, sizeof(double));
  }
  d9 = calloc(G+1, sizeof(double));
  if(!SkipDiploid){
    d10 = calloc(G+1, sizeof(double));
  }
  d11 = calloc(G, sizeof(double));
  d12 = calloc(G, sizeof(double));
  if(!SkipDiploid){
    d13 = calloc(G, sizeof(double));
    d14 = calloc(G, sizeof(double));
  }
  d15 = calloc(G, sizeof(double));
  d16 = calloc(G, sizeof(double));
  if(!SkipDiploid){
    d17 = calloc(G, sizeof(double));
    d18 = calloc(G, sizeof(double));
  }
  p1 = malloc((1ul<<L)*sizeof(double));
  p2 = malloc((1ul<<L)*sizeof(double));
  d9[0] = dist_haploid(p, Hpop[0]);           // distance between initial population point finite and infinite haploid
  if(!SkipDiploid){
    d10[0] = dist_n(p, Pop[0]);                 // distance between initial population point finite and infinite diploid
  }
  
  for(i = 0; i < 1ul<<L; i++){                // clone initial population
    p1[i] = p[i];
  }
  
  // open file for writing data about in between p_str and q_str or not
  sprintf(tfile, "b%02lu_n%06lu_eps%0.6lf_bet_hap_%02lu.dat", L, Nh, Epsilon, run);
  fp = fopen(tfile, "w");
  sprintf(tfile, "b%02lu_n%06lu_eps%0.6lf_bet_dip_%02lu.dat", L, Nd, Epsilon, run);
  fp2 = fopen(tfile, "w");
  sprintf(tfile, "b%02lu_eps%0.6lf_bet_inf_hap_%02lu.dat", L, Epsilon, run);
  fp3 = fopen(tfile, "w");
  double *dz1_str = calloc(1ul<<L, sizeof(double));
  double *dz2_str = calloc(1ul<<L, sizeof(double));
  vector_sub(p_str, p1_str, dz1_str, 1ul<<L);
  vector_sub(p_str, q1_str, dz2_str, 1ul<<L);
  if(dot_product(n_str, dz1_str, 1ul<<L) < 0.0 && dot_product(n_str, dz2_str, 1ul<<L) > 0.0){
    fprintf(fp, "%e\n\n", 1.0);
    fprintf(fp2, "%e\n\n", 1.0);
  }
  else{
    fprintf(fp, "%e\n\n", 0.0);
    fprintf(fp2, "%e\n\n", 0.0);
  }
  
  
  
  for(i = 0; i < G; i++){                     // move through G generations
    g_w(1, p1, p2);                           // evolution of population 1 generation at a time
    d1[i] = dist_p1p2_haploid(p2, p_str);     // distance between p2 and p_str haploid
    d2[i] = dist_p1p2_haploid(p2, q_str);     // distance between p2 and q_str haploid
    d11[i] = dist_p1p2_haploid(p2, p1_str);     // distance between p2 and p1_str haploid
    d12[i] = dist_p1p2_haploid(p2, q1_str);     // distance between p2 and q1_str haploid
    
    if(!SkipDiploid){
      d3[i] = dist_p1p2_diploid(p2, p_str);     // distance between p2 and p_str diploid
      d4[i] = dist_p1p2_diploid(p2, q_str);     // distance between p2 and q_str diploid  
      d13[i] = dist_p1p2_diploid(p2, p1_str);     // distance between p2 and p1_str diploid
      d14[i] = dist_p1p2_diploid(p2, q1_str);     // distance between p2 and q1_str diploid  
    }
    vector_sub(p2, p1_str, dz1_str, 1ul<<L);
    vector_sub(p2, q1_str, dz2_str, 1ul<<L);
    if(dot_product(n_str, dz1_str, 1ul<<L) < 0.0 && dot_product(n_str, dz2_str, 1ul<<L) > 0.0){
      fprintf(fp3, "%e\n\n", 1.0);      
    }
    else{
      fprintf(fp3, "%e\n\n", 0.0);      
    }
    
    // swap pointers to p1 and p2
    tptr = p1;
    p1 = p2;
    p2 = tptr;
    
    // finite haploid
    for(j = 0; j < Nh; j++){         // reproduce N offsprings      
      a = rnd(Nh); b = rnd(Nh);
      reproduce_haploid(a, b, j, Hpop[0], Hpop[1]);           // randomly two parents chosen; will be proportional to proportion
    }
    if(dot_product_haploid_n_and_pop_sub_p(p1_str, n_str, Hpop[1]) < 0.0 && dot_product_haploid_n_and_pop_sub_p(q1_str, n_str, Hpop[1]) > 0.0){
      fprintf(fp, "%e\n", 1.0);
    }
    else{
      fprintf(fp, "%e\n", 0.0);
    }      
  
    // set new generation as parent generation for next generation
    tmp_ptr = Hpop[0];
    Hpop[0] = Hpop[1];
    
    Hpop[1] = tmp_ptr;		
    merge_sort(Hpop[0], Nh); 
    // calculate distance to oscillating points and write to file
    d5[i] = dist_haploid(p_str, Hpop[0]);               // distance to 1st oscillating point
    d6[i] = dist_haploid(q_str, Hpop[0]);               // distance to 2nd oscillating point	 
    d15[i] = dist_haploid(p1_str, Hpop[0]);               // distance to 1st oscillating point without violation
    d16[i] = dist_haploid(q1_str, Hpop[0]);               // distance to 2nd oscillating point without violation
    
    if(!SkipDiploid){
      // finite diploid
      for(j = 0; j < Nd; j++){         // reproduce N offsprings      
	a = rnd(Nd); b = rnd(Nd);
	reproduce_diploid(a, b, j, Pop[0], Pop[1]);           // randomly two parents chosen; will be proportional to proportion
      }
      
      if(dot_product_diploid_n_and_pop_sub_p(p1_str, n_str, Pop[1]) < 0.0 && dot_product_diploid_n_and_pop_sub_p(q1_str, n_str, Pop[1]) > 0.0){
	fprintf(fp2, "%e\n", 1.0);
      }
      else{
	fprintf(fp2, "%e\n", 0.0);
      }
     
      // set new generation as parent generation for next generation
      tmp_ptr = Pop[0];
      Pop[0] = Pop[1];
      Pop[1] = tmp_ptr;		
      merge_sort(Pop[0], Nd); 
      
      // calculate distance to oscillating points and write to file
      d7[i] = dist_n(p_str, Pop[0]);               // distance to 1st oscillating point
      
      d8[i] = dist_n(q_str, Pop[0]);               // distance to 2nd oscillating point	  
      d17[i] = dist_n(p1_str, Pop[0]);               // distance to 1st oscillating point without violation
      d18[i] = dist_n(q1_str, Pop[0]);               // distance to 2nd oscillating point without violation
    }
    
    // distance between infinte and finite population after each generation
    d9[i+1] = dist_haploid(p1, Hpop[0]);           // distance between population point finite and infinite haploid
    if(!SkipDiploid){
      d10[i+1] = dist_n(p1, Pop[0]);                 // distance between population point finite and infinite diploid
    }
  }
  fclose(fp); fclose(fp2); fclose(fp3);
  free(dz1_str); free(dz2_str);
  
  // write haploid distances to file
  char str[100];
  sprintf(str,  "b%02lu_g%04lu_eps%0.6lf_osc_inf_haploid_%02lu.dat", L, G, Epsilon, run);             // add run to filename
  sprintf(fname, "%s", str);  // haploid distance data file 
  if(!(fp = fopen(fname, "w"))){
    printf("%s could not be opened!! Error!!\n", fname);
    exit(2);
  }
  for(i = 0; i < G; i++){
    fprintf(fp, "%lu  ", i);
    fprintf(fp, "%e  %e \n", d1[i], d2[i]);
  }
  fclose(fp); 
  // distance to p1_str and q1_str without violation
  sprintf(str,  "b%02lu_g%04lu_eps%0.6lf_osc_inf_haploid_wovio_%02lu.dat", L, G, Epsilon, run);             // add run to filename
  sprintf(fname, "%s", str);  // haploid distance data file 
  if(!(fp = fopen(fname, "w"))){
    printf("%s could not be opened!! Error!!\n", fname);
    exit(2);
  }
  for(i = 0; i < G; i++){
    fprintf(fp, "%lu  ", i);
    fprintf(fp, "%e  %e \n", d11[i], d12[i]);
  }
  fclose(fp); 
  
  // plot 
  if(!CLUSTER){
    sprintf(title, "b%lu g%lu gg%lu s%lu infinite haploid", L, G, Gg, Seed);    
    gp = popen ("gnuplot -persistent", "w"); // open gnuplot in persistent mode
    plot(gp, 0, fname, 3, title, "G", "d", 0, "" );
    fflush(gp);
    pclose(gp);
  }

  if(!SkipDiploid){
    // write inf diploid distances to file
    sprintf(str,  "b%02lu_g%04lu_eps%0.6lf_osc_inf_diploid_%02lu.dat", L, G, Epsilon, run);             // add run to filename
    sprintf(fname, "%s", str);  // diploid distance data file 
    if(!(fp = fopen(fname, "w"))){
      printf("%s could not be opened!! Error!!\n", fname);
      exit(2);
    }
    for(i = 0; i < G; i++){
      fprintf(fp, "%lu  ", i);
      fprintf(fp, "%e  %e \n", d3[i], d4[i]);
    }
    fclose(fp);
    // write inf diploid distances without violation to file
    sprintf(str,  "b%02lu_g%04lu_eps%0.6lf_osc_inf_diploid_wovio_%02lu.dat", L, G, Epsilon, run);             // add run to filename
    sprintf(fname, "%s", str);  // diploid distance data file 
    if(!(fp = fopen(fname, "w"))){
      printf("%s could not be opened!! Error!!\n", fname);
      exit(2);
    }
    for(i = 0; i < G; i++){
      fprintf(fp, "%lu  ", i);
      fprintf(fp, "%e  %e \n", d13[i], d14[i]);
    }
    fclose(fp);
    
    
    // plot inf diploid
    if(!CLUSTER){
      sprintf(title, "b%lu g%lu gg%lu s%lu infinite diploid", L, G, Gg, Seed);
      gp = popen ("gnuplot -persistent", "w"); // open gnuplot in persistent mode
      plot(gp, 0, fname, 3, title, "G", "d", 0, "" );
      fflush(gp);
      pclose(gp);
    }
  }
  // write distance for finite haploid to oscillatins points 
  sprintf(str, "osc_haploid_%02lu.dat", run);             // add run to filename
  prep_filename_hap(fname, str);
  if(!(fp = fopen(fname, "w"))){
    printf("%s could not be opened!! Error!!\n", fname);
    exit(2);
  }
   // write distances to file
   for(i = 0; i < G; i++){
    fprintf(fp, "%lu  ", i);
    fprintf(fp, "%e  %e\n", d5[i], d6[i]);
   }
  fclose(fp);  
  // write distance for finite haploid to oscillatins points without violation
  sprintf(str, "osc_haploid_wovio_%02lu.dat", run);             // add run to filename
  prep_filename_hap(fname, str);
  if(!(fp = fopen(fname, "w"))){
    printf("%s could not be opened!! Error!!\n", fname);
    exit(2);
  }
   // write distances to file
   for(i = 0; i < G; i++){
    fprintf(fp, "%lu  ", i);
    fprintf(fp, "%e  %e\n", d15[i], d16[i]);
   }
  fclose(fp);  
  
  if(!CLUSTER){
    sprintf(title, "b%lu g%lu n%lu gg%lu s%lu finite haploid", L, G, Nh, Gg, Seed);
    gp = popen ("gnuplot -persistent", "w"); // open gnuplot in persistent mode
    plot(gp, 0, fname, 3, title, "G", "d", 0, "" );
    fflush(gp);
    pclose(gp);
  }
  
  if(!SkipDiploid){
    // write distance for finite diploid to oscillatins points 
    sprintf(str, "osc_diploid_%02lu.dat", run);             // add run to filename
    prep_filename_dip(fname, str);
    if(!(fp = fopen(fname, "w"))){
      printf("%s could not be opened!! Error!!\n", fname);
      exit(2);
    }
    // write distances to file
    for(i = 0; i < G; i++){
      fprintf(fp, "%lu  ", i);
      fprintf(fp, "%e  %e\n", d7[i], d8[i]);
    }
    fclose(fp);  
    // write distance for finite diploid to oscillatins points without violation
    sprintf(str, "osc_diploid_wovio_%02lu.dat", run);             // add run to filename
    prep_filename_dip(fname, str);
    if(!(fp = fopen(fname, "w"))){
      printf("%s could not be opened!! Error!!\n", fname);
      exit(2);
    }
    // write distances to file
    for(i = 0; i < G; i++){
      fprintf(fp, "%lu  ", i);
      fprintf(fp, "%e  %e\n", d17[i], d18[i]);
    }
    fclose(fp);  
    
    if(!CLUSTER){
      sprintf(title, "b%lu g%lu n%lu gg%lu s%lu finite diploid", L, G, Nd, Gg, Seed);
      gp = popen ("gnuplot -persistent", "w"); // open gnuplot in persistent mode
      plot(gp, 0, fname, 3, title, "G", "d", 0, "" );
      fflush(gp);
      pclose(gp);
    }
  }
  
  // write haploid population distance between finite and infinite
  sprintf(str, "osc_haploid_dist_%02lu.dat", run);             // add run to filename
  prep_filename_hap(fname, str);
  if(!(fp = fopen(fname, "w"))){
    printf("%s could not be opened!! Error!!\n", fname);
    exit(2);
  }
  // write distances to file
  for(i = 0; i < G+1; i++){
    fprintf(fp, "%lu  ", i);
    fprintf(fp, "%e \n", d9[i]);
  }
  fclose(fp);  
  
  if(!CLUSTER){
    sprintf(title, "b%lu g%lu n%lu gg%lu s%lu haploid popn distance", L, G, Nh, Gg, Seed);
    gp = popen ("gnuplot -persistent", "w"); // open gnuplot in persistent mode
    plot(gp, 0, fname, 2, title, "G", "d", 0, "" );
    fflush(gp);
    pclose(gp);
  }
  
  if(!SkipDiploid){
    // write diploid population distance between finite and infinite
    sprintf(str, "osc_diploid_dist_%02lu.dat", run);      // add run to filename
    prep_filename_dip(fname, str);
    if(!(fp = fopen(fname, "w"))){
      printf("%s could not be opened!! Error!!\n", fname);
      exit(2);
    }
    // write distances to file
    for(i = 0; i < G+1; i++){
      fprintf(fp, "%lu  ", i);
      fprintf(fp, "%e \n", d10[i]);
    }
    fclose(fp);  
    
    if(!CLUSTER){
      sprintf(title, "b%lu g%lu n%lu gg%lu s%lu diploid popn distance", L, G, Nd, Gg, Seed);
      gp = popen ("gnuplot -persistent", "w"); // open gnuplot in persistent mode
      plot(gp, 0, fname, 2, title, "G", "d", 0, "" );
      fflush(gp);
      pclose(gp);
    }
  }
  
  free(d1); free(d2); free(d5); free(d6);   free(d9);  
  free(d11); free(d12); free(d15); free(d16); 
  free(p1); free(p2);
  if(!SkipDiploid){
    free(d3); free(d4);
    free(d7); free(d8);
    free(d10);
    free(d13); free(d14);
    free(d17); free(d18);
  }
}

void copy_mutation_crossover_to_dist()
{
  unsigned long i;
  double sm = 0, sc = 0;
  for(i = 0; i < 1ul<<L; i++){                     // copy mutation and crossover distributions to Cr and M distributions
    Cr->p[i] = Chi[i]; sc += Chi[i];
    M->p[i] = Mu[i]; sm += Mu[i];
  }
  
  initdist(Cr, sc);                                // initialize Cr distribution
  initdist(M, sm);                                 // initialize M distribution
}

// allocates memory
void setup()
{
  unsigned long i, g;
  C = pow(2., L/2.);   // constant
  U = (1ul<<L)-1;      // universe mask
  for (i = 256; i--; W[i] = ((ones(i)&1) ? -1. : 1.)/C);
  // allocate memory and calculate crossover and mutation 
  Chi = calloc((1ul<<L),sizeof(double));
  Mu = calloc((1ul<<L),sizeof(double));
  Chi_0 = calloc((1ul<<L),sizeof(double));
  Mu_0 = calloc((1ul<<L),sizeof(double));
  Cr = allocdist((1ul<<L));
  M = allocdist((1ul<<L));
  initOneCount();
  P0 = malloc((1ul<<L)*sizeof(double));
  P1 = malloc((1ul<<L)*sizeof(double));

  g = 0;
  while(!g){
    g = rnd(1ul<<L);                             // random value of g
  }
  Gg = g;
  chiDist(g); muDist(g);                         // initialize distributions for Mu and Chi with g = 12
  
  Mh  = malloc((1ul<<L)*sizeof(double *));
  for (i = 1ul<<L; i--;) Mh[i]  = malloc((1ul<<L)*sizeof(double));  
  
  GZ = malloc((1ul<<L)*sizeof(double));
  GW = malloc((1ul<<L)*sizeof(double)); 
   
}

// initializes memory for finite population string array
void init()
{  
  // allocate memory and initialize population 
  //if(!SkipDiploid){ 
    //printf("\nNd: %ld \n\n", Nd);
    Pop[0] = malloc(Nd*sizeof(unsigned long));
    Pop[1] = malloc(Nd*sizeof(unsigned long));
  //}
  Hpop[0] = malloc(Nh*sizeof(unsigned long));
  Hpop[1] = malloc(Nh*sizeof(unsigned long));
  //initpop(Pop[0]);     
}

void deinit()
{
  //if(!SkipDiploid){
    free(Pop[0]); free(Pop[1]); 
 // }
  free(Hpop[0]); free(Hpop[1]);
}

void cleanup()
{
  unsigned long i;  
  freedist(Cr); free(Chi); free(Chi_0);
  freedist(M); free(Mu); free(Mu_0);
  free(P0); free(P1); 
  for (i = 1ul<<L; i--;) free(Mh[i]);
  free(Mh);
  free(Z); free(GW); free(GZ); 
}



int main(int argc, char** argv)
{  
  if(argc^5){
    printf("\n Usage: ./vbet bits seed G runs\n bits: haploid bit length \n seed: seed for initialization of random number generator \n G: number of generations to simulate for \n \n runs: number of runs \n");
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
  G = (unsigned long)atol(argv[3]);      // no. of generations
  
  Runs = (unsigned long)atol(argv[4]);
  
  SkipDiploid = 0;                // skip finite diploid computation
  SkipHaploid = 0;                // skip finite haploid computation

  
  double *p_str, *q_str, *p1_str, *q1_str, *n_str;                 // oscillating points  
  time_t now; 
  unsigned long i, j, k;
  char seedfile[100];
  sprintf(seedfile, "b%ld_g%ld_vio_seed.dat", L, G);
  FILE *fp = fopen(seedfile, "w");
  fprintf(fp, "Run  Seed\n");
  for(j = 0; j < Runs; j++){
    if(seed == 0){
      now = time(0);
      Seed = ((unsigned long)now);     
    }
    else{
      Seed = seed;
    }
    fprintf(fp, "%ld  %ld\n", j, Seed);
    initrand(Seed);
    p_str = calloc((1ul<<L),sizeof(double));  // oscillating point 1
    q_str = calloc((1ul<<L),sizeof(double));  // oscillating point 2
    p1_str = calloc((1ul<<L),sizeof(double));  // oscillating point 1 without violation
    q1_str = calloc((1ul<<L),sizeof(double));  // oscillating point 2 without violation
    n_str = calloc((1ul<<L),sizeof(double));   // unit vector
    setup();        
    // initial population calc begin
    // generate random distribution in P0
    Nh = N0;
    Nd = N0;
    init();
    install_pop_distribution(P0);                           // install haploids proportion
    generate_fin_hap_init_pop_from_random_px(P0, Hpop[0]);  // generate finite haploid initial population from P0
    merge_sort(Hpop[0], Nh);  
    calc_px_from_finite_hap_pop(P0, Hpop[0]);               // initial population vector constant for all finite size population in simulation
    deinit();  
    // initial population calc end
    //mixing matrix before violation
    walsh(Mh, 0);                                          // calculate and install values in mixing matrix in walsh basis (Mhat)
    comp_periodic_orbits(p1_str, q1_str, P0);              // computes p_str and q_str oscillating points before violation
    // compute unit vector n = (p_str - q_str)/(|| p_str - q_str ||)
    compute_unit_vector(p1_str, q1_str, n_str, 1ul<<L);    
			  
    //display_p(Chi, "before violation chi:"); 
    for(k = 0; k < sizeof(Err)/sizeof(double); k++){
      copy_array(Chi_0, Chi, 1ul<<L);
      copy_array(Mu_0, Mu, 1ul<<L);
      //display_p(Chi, "after copy from Chi_0 chi:"); 
      Epsilon = Err[k];
      // introduce violation in Mu or Chi
      if(VIO_MU_CHI == 1){
	violate_mu();                              // violates oscillating condition in Mu
      }
      else if(VIO_MU_CHI == 2){
	violate_chi();                             // violates oscillating condition in Chi
	//display_p(Chi, "after violation chi:"); 
      } 
      // calculation mixing matrix after violation
      walsh(Mh, 0);                                          // calculate and install values in mixing matrix in walsh basis (Mhat) after violation
      copy_mutation_crossover_to_dist();                     // fills up crossover distribution and mutation distribution     
      
      comp_periodic_orbits(p_str, q_str, P0);                 // computes p_str and q_str oscillating points after violation
      //display_p(p_str, "p_str:"); display_p(q_str, "q_str:");           
      
      for(i = 0; i < sizeof(Ni)/sizeof(unsigned long); i++){                        // through all sizes of finite population
	Nh =  (unsigned long)pow(N0,2)*Ni[i]; 
	Nd = Nh;
	init();                                               // initializes memory for pop, M and Cr and also installs values for these  
	
	generate_fin_hap_pop_from_pvector(P0, Hpop[0]);       // finite haploid population from P0
	
	if(SkipDiploid == 0){
	  generate_fin_dipop_from_px(P0, Pop[0]);               // finite diploid population generation from P0      
	}
	osc_all_n_dist(P0, p_str, q_str, p1_str, q1_str, j, n_str);
	
	deinit();     
	
      }
    }
    free(p_str); free(q_str); free(p1_str); free(q1_str); free(n_str);
    cleanup();
  }
  walsh(Mh, -1);  
  fclose(fp);  
  return 0;
}

/*dependent files:
rand.c
sort.c
*/

/* Compile and Run:
gcc -O2 -march=native -o vbet vbet.c -lm
./vbet bits seed generations runs
*/


 
