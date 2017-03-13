/********************************************************************
Random numbers
********************************************************************/
#define TWO_32   (4294967296.0)

/* random generator */
#define rndm() ((++rndx>54)? rtab[rndx=nrndm()]: rtab[rndx])

#define U01()   (rndm()/TWO_32)     /* random in interval [0,1) */
#define U(x)    (U01()*(x))         /* random in interval [0,x) */
#define rnd(n)  ((unsigned)U((n)))  /* random from set {0..n-1} */

/* for generating random number i with probability p[i] */
typedef struct dist 
{ 
  double *p;
  int    *a;
  int     n;
  double  m1;
  double  m2;
} dist;

static unsigned rndx, rtab[55];

static int nrndm()
{
  int i;

  for (i =  0; i < 24; i++) rtab[i] -= rtab[i+31];
  for (i = 24; i < 55; i++) rtab[i] -= rtab[i-24];
  return 0;
}

/* initialize the 32 bit random number generator with seed j */
static void initrand(unsigned j) 
{
  int h,i,k;

  for (rtab[54] = j |= (k = i = 1); i < 55; i++)
    h = (21*i)%55, rtab[--h] = k, k = j - k, j = rtab[h];
  while (i--) nrndm();
  rndx = 0; 
}

static void error_rnd(char *s)
{
  printf("%s\n",s);
  exit(1);
}

/* allocate a distribution over {0..n-1} */
static dist *allocdist(unsigned n) 
{
  dist *d;

  if (!(d    = malloc(sizeof(dist))))          
    error_rnd("malloc (allocdist: d)");
  d->n = n;
  if (!(d->a = malloc(d->n * sizeof(int  )))) 
    error_rnd("malloc (allocdist: d->a)");
  if (!(d->p = malloc(d->n * sizeof(double)))) 
    error_rnd("malloc (allocdist: d->p)");
  return d;
}

static void freedist(dist *d)
{
  free(d->a);
  free(d->p);
  free(d);
}

#define getsmall { while (p[j] >= q) if ((++j) == stop) goto end; t = j++; }
#define getlarge while (p[k] <  q) if ((++k) == stop) goto cleanup;

/* initialize the distribution d  */
static dist *initdist(dist *d, double s) 
{
  /*
    d->p must have d->n elements which sum to s on entry to initdist
    d->p and d->a are overwritten by the initialization process
  */
  int j,k,t,stop,*a;  double q,*p;

  stop = d->n, q = s/stop, j = k = 0;

  d->m1 = stop/TWO_32;
  d->m2 = s/(stop * TWO_32);

  a = d->a;
  p = d->p;

  getsmall; getlarge;

 loop:
    
  a[t]  = k;
  p[k] += p[t] - q;

  if (p[k] >= q) { 
    if (j == stop) goto end;
    getsmall;
    goto loop;
  }
  t = k++;
  if (k == stop) goto cleanup;
  if (j < k) getsmall;
  getlarge;
  goto loop;

 cleanup:

  a[t] = t;
  while (j < stop) { a[j] = j; j++; }

 end: 
  return d;
}
#undef getsmall
#undef getlarge

/* returns element from {0..d->n-1} according to d */
static unsigned drand(dist *d) 
{
  unsigned j = rndm()*d->m1;

  if (rndm()*(d->m2) < d->p[j]) return j;
  return d->a[j];
}

