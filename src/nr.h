#ifndef _NR_H_
#define _NR_H_

#ifdef SINGLE
#define REAL REAL
#else
#define REAL double
#endif

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "nrutil.h"

#define FUNC(x) ((*func)(x))
REAL trapzd(REAL (*func)(REAL), REAL a, REAL b, int n);

void polint(REAL xa[], REAL ya[], int n, REAL x, REAL *y, REAL *dy);

#define JMAX 20        
#define JMAXP (JMAX+1)
#define K 5            
REAL qromb(REAL (*func)(REAL), REAL a, REAL b);

#define ITMAX 100                       // maximum number of iterations  
#define EPS 3.0e-7                      // relative accuracy  
#define FPMIN 1.0e-30                   // "smallest" floating point number 

void fit(double x[], double y[],        // fit a straight line y=a+bx 
         int ndata, double sig[],       // to the data x[1...data],y[1...data] 
         int mwt, double *a, double *b, // with standard dev sig[1..ndata] 
         double *siga, double *sigb,    // return a,b and their uncertainties 
         double *chi2, double *q);      // siga, sigb, the chi2 and the  
                                        // goodnes-of-fit probability q. 
                                        // if mwt=0 ignores sig[] and q=1.0 
 
double gammq(double a, double x);       // the incomplete gamma function Q(a,x)
double gammp(double a, double x);       // the incomplete gamma function P(a,x)

void gcf(double *gammcf, double a,      // the incomplete gamma function P(a,x)
         double x, double *gln);        // completeness fraction representation

void gser(double *gamser, double a,     // the incomplete gamma function P(a,x)
          double x, double *gln);       // series representation 
 
double gammln(double xx);               // ln(Gamma(xx)) 
 
#endif /* _NR_H_ */
