#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h> 
#include <math.h>
#include <strings.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include "rand.c"

unsigned long L;  // no. of bits to represent chromosome
double *Chi;       // crossover probability array
double *Mu;        // mutation probability array
double C, W[256];  // pow(2.,L/2.), pow(2.,-L/2.)*pow(-1.,ones(index))  
unsigned long U;   // universal mask ((1ul<<L)-1)  
double **Mh;       // Mh = Mhat; mixing matrix in walsh basis
double *P0, *P_str, *Q_str;
FILE *fp1, *fp2, *Log;
double D;

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

static int ones(unsigned long x) // returns |x| the number of 1 bits 
{
  int a;
  for (a = x&1; x >>= 1; a += x&1);
  return(a);
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
  for(i = 1; i < (1ul<<L); i++){
    g = si[i];
    x = x_g(g);
    if( fabs(x) < 1.0){
      ph[g] = ( x*y_g(ph,g) + y_g(qh,g) )/( 1.0 -x*x );
      qh[g] = ( x*y_g(qh,g) + y_g(ph,g) )/( 1.0 -x*x );
    }
    else{
      ph[g] = pw[g];
      qh[g] = mpw[g];
    }
      
  }
  walsh_v(ph, p_str); walsh_v(qh, q_str);  // periodic orbits in standard basis
  
  free(pw); free(mpw); free(ph); free(qh); free(si); free(sv);
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

void showBits(FILE *f, void *v, int s)
{
  int i = s;

  for (i = 0; i < s; i++)
    fprintf(f,"%02x", ((unsigned char *)v)[i]);
  fprintf(f,"\n");
}

#define showHex(f,x) showBits(f, &(x), sizeof(x));


void readBits(FILE *f, void *v, int s)
{
  int i = s;  unsigned int j;  char c;

  for (i = 0; i < s; i++){
    fscanf(f, "%02x", &j);
    ((unsigned char *)v)[i] = j;
  }
  fscanf(f,"%c", &c);
}

#define readHex(f,x) readBits(f, &(x), sizeof(x));

void readDoubleArray(FILE *f, char *s, double *a, int n)
{
  int i, j = strlen(s);  char c;

  while (j--) fscanf(f,"%c", &c); fscanf(f,"%c", &c);
  for (i = 0; i < n; i++)
    readHex(f,a[i]);
}

void showDoubleArray(FILE *f, char *s, double *a, int n)
{
  int i;

  fprintf(f, "%s\n", s);
  for (i = 0; i < n; i++)
    showHex(f,a[i]);
}

void copy_p(double *d, double *s, unsigned long size)
{
  unsigned long i;
  for(i = 0; i < size; i++){
    d[i] = s[i];
  }
}

void display_p(FILE *f, double *p, char *str, unsigned int sum)                 // displays haploids' proportion p
{
  unsigned long i;
  double s = 0;
  fprintf(f, "%s\n", str);
  for(i = 0; i < 1ul<<L; i++){
    fprintf(f, "%.16lf ", p[i]);
    if(sum == 1){
      s += p[i];
    }
  }
  fprintf(f, "\n");
  if(sum==1) fprintf(f, "sum of %s %.16lf\n", str, s);
}

void new_handler (const char * reason, const char * file, int line, int gsl_errno)
{
  printf("new_handler:\n");
  //rewind(Log);
  memset(Chi, 0, (1ul<<L)*sizeof(double));
  //readDoubleArray(Log,"\nchi:", Chi, 1ul<<L);
  display_p(stdout, Chi, "chi:", 1);
  showDoubleArray(stdout, "\nchi:", Chi, 1ul<<L);
  exit(0);
}

void setup()
{
  unsigned long i;
  C = pow(2., L/2.);   // constant
  U = (1ul<<L)-1;      // universe mask
  for (i = 256; i--; W[i] = ((ones(i)&1) ? -1. : 1.)/C);
  // allocate memory and calculate crossover and mutation 
  Chi = calloc((1ul<<L),sizeof(double));
  Mu = calloc((1ul<<L),sizeof(double)); 
  P0 = calloc((1ul<<L),sizeof(double)); 
  P_str = calloc((1ul<<L),sizeof(double)); 
  Q_str = calloc((1ul<<L),sizeof(double)); 
  initOneCount();  
  Mh  = malloc((1ul<<L)*sizeof(double *));  
  for (i = 1ul<<L; i--;) Mh[i]  = malloc((1ul<<L)*sizeof(double));
  
  // walsh(Mh);                                       // calculate and install values in mixing matrix in walsh basis (Mhat)   
  Log = fopen("LogFile","r");
  //gsl_set_error_handler(new_handler);
}

void cleanup()
{
  unsigned long i;  
  free(P0); free(P_str); free(Q_str);
  free(Mu); free(Chi);
  for (i = 1ul<<L; i--;) free(Mh[i]);
  free(Mh); 
  fclose(Log);
}

void set_chi_dist_with_g(double *chi, unsigned long g, double *y, unsigned long c)
{
  unsigned long h, i, k;
  double s;
  memset(chi, 0, (1ul<<L)*sizeof(double));  
  h = bar(g); 
  for(s = 0, i=0, k = 0; k < (1ul<<L); k++){     // through all possible haploids       
    if( k == (h&k) ){                             // check if k in bar(g)R
      if(i < c-1){
	chi[k^g] = y[i++];
	chi[k] = y[i++];     	
      }
      else{
	printf("Error i and c!!\n");
	exit(1);
      }
    } 
    s += chi[k];
  }
  //for(k = 0; k < (1ul<<L); k++) chi[k] /= s;
}

void set_mu_dist_with_g(double *mu, unsigned long g, double *x, unsigned long m)
{
  double s;
  unsigned long k, i, j;
  memset(mu, 0, (1ul<<L)*sizeof(double));
  for(s = 0, k = 0, j = 0; j < (1ul<<L); j++){            // through all possible haploids
    i = g&j;                                       // assumes 32 bit maximum
    i ^= (i >> 16);
    i ^= (i >> 8);
    if(ones(i&255)&1){                              // if odd, assign random value
      if(k < m){
	mu[j] = (double)x[k++];            
      }
      else{
	printf("Error in m and k!! \n");
	exit(1);
      }
    }    
    s += mu[j];
  }
  //for(j = 0; j < 1ul<<L; j++) mu[j] /= s;
}

double my_f (const gsl_vector *v, void *params)
{
  double *x, *y, sx, sy, d;
  unsigned long g, m, c, i, j;
  double *p_str = malloc((1ul<<L)*sizeof(double));  // oscillating point 1
  double *q_str = malloc((1ul<<L)*sizeof(double));  // oscillating point 2
  unsigned long *dp = (unsigned long *)params;
  g = dp[0];                                    // g
  m = dp[1];                                    // no. of non zero mu distributions
  c = dp[2];                                    // no. of non zero chi distributions
  x = calloc(m,sizeof(double));
  y = calloc(c,sizeof(double));
  double sm=0, sc = 0;
  for(sx=0, i=0; i<m; i++){
    x[i] = (double)gsl_vector_get(v, i);
    x[i] = pow(x[i],2);                           // squaring for non negative values
    sx += x[i];
    sm += x[i];
  }
  for(sy=0, i=0; i<c; i++){
    y[i] = (double)gsl_vector_get(v, m+i);
    y[i] = pow(y[i], 2);                           // squaring
    sy += y[i];
    sc += y[i];
  }
  // normalization of squared values
  for( i=0; i<m; i++){
    x[i] /= sx;
  }
  for(i=0; i<c; i++){
    y[i] /= sy;
  } 

  set_mu_dist_with_g(Mu, g, x, m);
  set_chi_dist_with_g(Chi, g, y, c);

#if 0
  rewind(Log);
  memset(Chi, 0, (1ul<<L)*sizeof(double));
  memset(Mu, 0, (1ul<<L)*sizeof(double));
  memset(P0, 0, (1ul<<L)*sizeof(double));
  readDoubleArray(Log,"p0", P0, 1ul<<L);
  readDoubleArray(Log,"mu", Mu, 1ul<<L);
  readDoubleArray(Log,"chi", Chi, 1ul<<L);
#endif
  //display_p(stdout, Chi, "chi:", 1);  
  
  // initialize mixing matrix
  walsh(Mh, 0);
  //display_p(stdout, Mu, "mu:", 1);
  // compute p_str and q_str
  //display_p(stdout, P0, "p0:", 1);               // displays haploids' proportion p  
  
  comp_periodic_orbits(p_str, q_str, P0);
  //display_p(stdout, p_str, "p_str:", 0); 
  //display_p(stdout, q_str, "q_str:", 0); 
  d = dist_p1p2_haploid(p_str, q_str);  
  copy_p(P_str, p_str, 1ul<<L);
  copy_p(Q_str, q_str, 1ul<<L);
  free(p_str); free(q_str);
  free(x); free(y);
  
  //printf("d: %lf a: %lf b: %lf\n", d, pow(sm-1.0, 2), pow(sc-1.0, 2));
  //return 1.0/d + pow(sm - 1.0, 2) + pow(sc-1.0, 2); 
  D = d;                                                                 // store distance
  return -d*d + (sm-1)*(sm-1) + (sc-1)*(sc-1);
}

int main(int argc, char **argv)
{
  if(argc^4){
    printf("Usage: ./opt bits seed g\n bits: haploid bit length \n seed: seed for initialization of random number generator \n g: a haploid string\n");
    exit(1);
  }
  L = atol(argv[1]);                     // bit length
  if (L > 32){
    printf("Can't have more thqan 32 bits\n");
    return 1;
  }
  unsigned long seed = atol(argv[2]);
  unsigned long g = (unsigned long)atol(argv[3]);      // no. of generations
  
  time_t now; 
  if(seed == 0){
    now = time(0);
    seed = ((unsigned long)now); 
  }
  initrand(seed);
  setup();         
  install_pop_distribution(P0);             // install haploids proportion
  
  unsigned long i, j, h, k, mc, cc;
  fp1 = fopen("Mu.dat", "w");
  fp2 = fopen("Chi.dat", "w");
  // get mu distribution non zero places count
  for(mc = 0, j = 0; j < (1ul<<L); j++){            // through all possible haploids
    i = g&j;                                       // assumes 32 bit maximum
    i ^= (i >> 16);
    i ^= (i >> 8);
    if(ones(i&255)&1){                              // if odd, assign random value
      mc++;      
    }    
  }
  // get chi distribution non zero places count
  h = bar(g); 
  for(cc = 0, k = 0; k < (1ul<<L); k++){               // through all possible haploids    
    if( k == (h&k) ){                             // check if k in bar(g)R
      cc = cc+2;
    }    
  }

  // prepare gsl function to invoke
  size_t np = mc + cc; //dimension of the problem
  unsigned long par[3] = {g, mc, cc};
  printf("%lu %lu %lu \n", g, mc, cc);
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex; //nedler-mead simplex algorithm
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;
  size_t iter = 0, ic;
  int status;
  double size;
  
  /* Initial vertex size vector */
  ss = gsl_vector_alloc (np);
  gsl_vector_set_all (ss, sqrt((double)(1ul<<(2*L))));   // set all step size to 0.001
  /* Starting point */
  x = gsl_vector_alloc (np);
  for(i = 0; i < mc; i++){          // mu distribution starting values
    gsl_vector_set (x, i, U01());
  }
  for(i = 0; i < cc; i++){          // chi distribution starting values
    gsl_vector_set (x, mc+i, U01());
  }
  
  /* Initialize method and iterate */
   minex_func.f = &my_f;
  minex_func.n = np;
  minex_func.params = (void *)&par;
  s = gsl_multimin_fminimizer_alloc (T, np);
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
  do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);
      if (status)
	break;
      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, 1e-16);

      if (status == GSL_SUCCESS)
	{
	  printf ("converged to minimum at\n");
	}
      /*for (ic = 0; ic < np; ic++)
	{
	  printf ("%10.3e ", gsl_vector_get (s->x, ic));
	}
      printf ("f() = %7.3f size = %.3f\n", s->fval, size);
      printf("\n\n iter: %ld\n\n", iter);*/
    }
  while (status == GSL_CONTINUE && iter < 1000);
  printf ("%ld \n", iter);
  char fname[200];
  sprintf(fname, "b%lug%lus%lu_maxd.dat", L, g, seed);
  FILE *f = fopen(fname, "w");
  fprintf(f, "d: %.16lf\n", D);
  fprintf (f, "f() = %.16f size = %.16f\n", s->fval, size);
  fprintf(f, "\n\n iter: %ld\n\n", iter);
  display_p(f, Mu, "Mu:", 0);
  display_p(f, Chi, "Chi:", 0);
  display_p(f, P_str, "p_str:", 0);
  display_p(f, Q_str, "q_str:", 0);
  display_p(f, P0, "P0:", 0);
  fclose(f);
  //clean up
  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);
  fclose(fp1); fclose(fp2);
  walsh(Mh, -1); // free static z memory alloc
  cleanup();
  
  return 0;
}

/* Compile and Run:
gcc -O2 -march=native -o opt opt.c -lm -lgsl -lgslcblas
./opt bits seed g 
*/



