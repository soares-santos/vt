/* ============================================================ *
 * mathstuff.h							*
 * Martin Kilbinger, NR 2004					*
 * Contains numerical recepies functions, constants, etc.       *
 * ============================================================ */

#ifndef __MATHSTUFF_H
#define __MATHSTUFF_H

#include <stdio.h>


#ifdef NR_COMPLEX_H_
#include "nrcomplex.h"
#endif

/* ============================================================ *
 * definition of the floating type used for all float type	*
 * variables when changing, check: fscanf, nr routines.	        *
 * Note: This has been double for a long time.                  *
 * ============================================================ */
typedef double real;

#ifndef pi
#define pi     3.14159265358979323846
#endif

#ifndef pi_sqr
#define pi_sqr 9.86960440108935861883
#endif

#ifndef twopi
#define twopi  6.28318530717958647693
#endif

#define pihalf 1.57079632679489661923
#define threepi    9.4247779607693793
#define sqrt2  1.41421356237309504880
#define sqrtm2 0.70710678118654752440

/* one arc minute in rad */
#define  ARCMIN 2.90888208665721580398e-04
#define  ARCSEC 4.84813681109535984270e-06

/* the following were static, gave 'undefined' warnings */
static real darg __attribute__((unused)), maxarg1 __attribute__((unused)), 
  maxarg2 __attribute__((unused)), darg2 __attribute__((unused)), 
  absarg1 __attribute__((unused)), absarg2 __attribute__((unused));

static int iminarg1 __attribute__((unused)), iminarg2 __attribute__((unused)),
  iarg __attribute__((unused));

#define IMIN(a,b) (iminarg1=(a), iminarg2=(b), (iminarg1) < (iminarg2) ? (iminarg1) : (iminarg2))
#define IMAX(a,b) (iminarg1=(a), iminarg2=(b), (iminarg1) > (iminarg2) ? (iminarg1) : (iminarg2))
#define FMAX(a,b) (maxarg1=(a), maxarg2=(b), (maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))
#define FMIN(a,b) (maxarg1=(a), maxarg2=(b), (maxarg1) < (maxarg2) ? (maxarg1) : (maxarg2))
#define ISQR(a) ( (iarg=(a))==0.0 ? 0.0 : iarg*iarg )
#define DSQR(a) ( (darg=(a))==0.0 ? 0.0 : darg*darg )
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

/* absolute value of x-y (=|x-y|) */
#define ABS22(x, y) (sqrt( DSQR( (x)[0] - (y)[0]) + DSQR( (x)[1] - (y)[1]) ))

/* absolute square of x-y (=|x-y|^2) */
#define abs2sqr(x, y) ( absarg1=(x)[0]-(y)[0], absarg2=(x)[1]-(y)[1], absarg1*absarg1 + absarg2*absarg2 )

/* square of the absolute value of x (=|x|^2) */
#define abssqr(x) ((x)[0]*(x)[0] + (x)[1]*(x)[1])

/* (a-b)*(a-c) for a,b,c 2-d vectors */
#define scp(a, b, c) ( ((a)[0]-(b)[0])*((a)[0]-(c)[0]) + ((a)[1]-(b)[1])*((a)[1]-(c)[1]) )

/* scalar product (x,y) */
#define sp(x,y) ((x)[0]*(y)[0]+(x)[1]*(y)[1])

double Dsqr(double x);
double abs22(const double x[], const double y[]);
double abs2sqr_f(const double *x, const double *y);
double absnsqr(const double *x, const double *y, int ndim);
void vdiff(const real *x, const real *y, real *z);
void out_error(const char *);
FILE *fileopen(const char [], const char []);
void fileclose(FILE *);
real interpol(real *f, int n, real a, real b, real dx, real x, real lower,
	      real upper);
real interpol2d(real **f, int nx, real ax, real bx, real dx, real x,
		int ny, real ay, real by, real dy, real y, real lower,
		real upper);
void nrerror(const char error_text[]);
real *vector(long nl, long nh);
int *ivector(long nl, long nh);
void free_vector(real *v, long nl, long nh);
void free_ivector(int *v, long nl, long nh);
real **matrix(long nrl, long nrh, long ncl, long nch);
int  **imatrix(long nrl, long nrh, long ncl, long nch);
void free_matrix(real **m, long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void polint(real xa[], real ya[], int n, real x, real *y, real *dy);

/* Romberg integration */
real trapzd2d(real (*func)(real, real), real ax, real bx, real ay, real by, int n);

real trapzdberg(real (*func)(real, const real*), const real *intpar, real a, real b,
		int n, double *s);
real trapzd (real (*func)(real, const real*), const real *intpar, real a, real b, int n);
real trapzd1(real (*func)(real, const real*), const real *intpar, real a, real b, int n);
real trapzd2(real (*func)(real, const real*), const real *intpar, real a, real b, int n);
real trapzd3(real (*func)(real, const real*), const real *intpar, real a, real b, int n);

real qromberg(real (*func)(real, const real*), const real *intpar, real a, real b, real EPS);
real qromb (real (*func)(real, const real*), const real *intpar, real a, real b, real EPS);
real qromb1(real (*func)(real, const real*), const real *intpar, real a, real b, real EPS);
real qromb2(real (*func)(real, const real*), const real *intpar, real a, real b, real EPS);
real qromb3(real (*func)(real, const real*), const real *intpar, real a, real b, real EPS);
real qromo(real (*func)(real, const real*), const real *intpar, real a, real b, 
	   real (*choose)(real(*)(real, const real*), const real *, real, real, int),
	   real EPS);
real qrombergo(real (*func)(real, const real*), const real *intpar, real a, real b, 
	   real (*choose)(real(*)(real, const real*), const real *, real, real, int, real *),
	   real EPS);

real midpnt(real (*func)(real, const real*), const real *intpar, real a, real b, int n);
real midpntberg(real (*func)(real, const real*), const real *intpar, real a, real b, 
		int n, real *s);
real midinf(real (*funk)(real, const real*), const real *intpar, real aa, real bb, int n);
real midsql(real (*funk)(real, const real*), const real *intpar, real aa, real bb, int n);
real midsqu(real (*funk)(real, const real*), const real *intpar, real aa, real bb, int n);

/* === */

void sobseq(int *n, real x[]);
void random_init(char name[]);
real dfridr(real (*func)(real,real), real x, real h, real *err, real aa);
real dfridr1(real (*func)(real), real x, real h, real *err);
void spline(real x[], real y[], int n, real yp1, real ypn, real y2[]);
void splint(real xa[], real ya[], real y2a[], int n, real x, real *y);
void spli_pol(real *y, int n, real dx, real yp1, real ypn, real *y2);
real spli_value(real x, real *y, real *y2, int n, real xmin, real dx);

real ms_sinc(real x);
real ms_gammln(real xx);
real ms_bessj0(real x);
real bessj1(real x);
real bessj(int n, real x);
real BesselJ(int n, real x);
real bessy0(real x);
real bessy1(real x);
real bessi0(real x);
real bessi1(real x);
real bessi(int n, real x);
void choldc(real **a, int n, real p[]);
real **cholesky(real **, int);
void ludcmp(real **a, int n, int *indx, real *d);
void lubksb(real **a, int n, int *indx, real b[]);
real gasdev();
real ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_f3tensor(real ***t, long nrl, long nrh, long ncl, long nch,
		   long ndl, long ndh);
real factln(int);
real fact(int);
real bico(int, int);
real rtsec(real (*func)(real,real), real, real, real, real);
real rtbis(real (*func)(real,real), real, real, real, real);
real rtbis1(real (*func)(real), real, real, real);

real minvec(real *, int);
real maxvec(real *, int);
void permute(real *, int);

void lfit(real x[], real y[], real sig[], int ndat, real a[], int ia[],
	  int ma, real **covar, real *chisq, void (*funcs)(real, real [], int));
void covsrt(real **covar, int ma, int ia[], int mfit);
void gaussj(real **a, int n, real **b, int m);
void jacobi(real **a, int n, real d[], real **v, int *nrot);
void balanc(double **a, int n);
void elmhes(double **a, int n);
void hqr(double **a, int n, double wr[], double wi[]);

real **matrix_multiply(const real **a, const real **b, int N);
real *matrix_vector_multiply(const real **A, const real *x, int N);
real **transpose(const real **, int);
real inverse(real **, int);
real **matrix_copy(const real **, int);
void unity_test(const real **, const real **, int);
void out_matrix(char *name, const real **a, int n);
void linbcg(unsigned long n, const double **A, const double **Aguess, const double b[],
	    double x[], int itol, double tol, int itmax, int *iter, double *err,
	    int guess_is_diag);
double snrm(unsigned long n, const double sx[], int itol);
void atimes(const double** A, const double x[], double r[],
	    unsigned long n, int itrnsp);
void asolve(const double **Aguess, const double b[], double x[], unsigned long n,
	    int itrnsp, int guess_is_diag);

double gaussfunction(double xsqr, double sigmasqr);

#ifdef NR_COMPLEX_H_
dcomplex **cmatrix(long nrl, long nrh, long ncl, long nch);
void free_cmatrix(dcomplex **m, long nrl, long nrh, long ncl, long nch);
dcomplex ***c3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_c3tensor(dcomplex ***t, long nrl, long nrh, long ncl, long nch,
		   long ndl, long ndh);
void cisitab(real x, real *ci, real *si);
void cisi(real x, real *ci, real *si);
#endif

real brent(real ax, real bx, real cx, real (*f)(real), real tol,
	   real *xmin);
void avevar(real data[], unsigned long n, real *ave, real *var);
void gcf(real *gammcf, real a, real x, real *gln);
void gser(real *gamser, real a, real x, real *gln);
real gammq(real a, real x);
void mnbrak(real *ax, real *bx, real *cx, real *fa, real *fb, real *fc,
	    real (*func)(real));
real zbrent(real (*func)(real), real x1, real x2, real tol);
void fit(real x[], real y[], int ndata, real sig[], int mwt, real *a,
	 real *b, real *siga, real *sigb, real *chi2, real *q);
real chixy(real bang);


void fitexy(real x[], real y[], int ndat, real sigx[], real sigy[],
	    real *a, real *b, real *siga, real *sigb, real *chi2, real *q);

real Pleg(int l, real x);

void donothing(real x);

#endif
/* __MATHSTUFF_H */

