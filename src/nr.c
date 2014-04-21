/*****************************************************************************/
/* Functions from Numerical Recipes                                   */
/*****************************************************************************/

#include "nr.h"

/*#define FUNC(x) ((*func)(x))*/
REAL trapzd(REAL (*func)(REAL), REAL a, REAL b, int n){ 
  /* Chapter 4.2                                                             */
  REAL x, tnm,sum,del;
  static REAL s;
  int it,j;
  if (n ==1) return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
  else {
    for (it=1,j=1;j<n-1;j++) it <<= 1;
    tnm=it;
    del=(b-a)/tnm;
    x=a+0.5*del;
    for (sum=0.0,j=1;j<=it;j++,x+=del) sum +=FUNC(x);
    s=0.5*(s+(b-a)*sum/tnm);
    return s;
  }
}

void polint(REAL xa[], REAL ya[], int n, REAL x, REAL *y, REAL *dy){
  /* Chapter 3.1  - Given arrays xa[1..n] and ya[1..n], and given a value x, 
     this routine returns a value y, and an error estimate dy. If  P(x) is 
     the polynomial of degree N -1 such that P(xa_i)=ya_i,i=1,...,n, then the 
     returned value is y=P(x).                                               */
  int i,m,ns=1;
  REAL den,dif,dift,ho,hp,w;
  REAL *c,*d;
  dif=fabs(x-xa[1]);
  c=vector(1,n);
  d=vector(1,n);
  for (i=1;i<=n;i++) { /* Find the index ns of the closest table entry,      */
    if ( (dift=fabs(x-xa[i])) < dif) {
      ns=i;
      dif=dift;
    }
    c[i]=ya[i];        /* and initialize the tableau of c's and d's.         */
    d[i]=ya[i];
  }
  *y=ya[ns--];         /* This is the initial approximation to y.            */
  for (m=1;m<n;m++) {  /* For each column of the tableau,                    */
    for (i=1;i<=n-m;i++) { /* loop over current c's and d's and update them. */
      ho=xa[i]-x; 
      hp=xa[i+m]-x;
      w=c[i+1]-d[i];
      if ( (den=ho-hp) == 0.0) nrerror("Error in routine polint");
      /* Error occurs only if 2 input xa's are identical to within roundoff. */
      den=w/den;
      d[i]=hp*den;     /* Here the c's and d's are updated.                  */
      c[i]=ho*den;
    }
    *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
    /* After each column in the tableau is completed, we decide which        
       correction, c or d, we want to add to our accumulating value of y,    
       i.e., which path to take through the tableau - forking up or down.    
       We do this in such a way as to take the most "straight line" route    
       through the tableau to its apex, updating ns accordingly to keep track
       of where we are. This route keeps the partial approximations centered 
       (insofar as possible) on the target x. The last dy added is thus the  
       error indication.                                                     */
  }
  free_vector(d,1,n);
  free_vector(c,1,n);
}

/*#define EPS 1.0e-6*/  /* fractional accuracy desired                       */
/*#define JMAX 20*/     /* total number of steps                             */
/*#define JMAXP (JMAX+1)                                                     */
/*#define K 5*/         /* number of points used in the extrapolation        */
REAL qromb(REAL (*func)(REAL), REAL a, REAL b){
  /* Chapter 4.3 - Returns the integral of the function func from a to b. 
     Integration is performed by Romberg's method of order 2K, where, 
     K=2 is Simpson's rule.                                                  */
  void polint(REAL xa[], REAL ya[], int n, REAL x, REAL *y, REAL *dy);
  REAL trapzd(REAL (*func)(REAL), REAL a, REAL b, int n);
  void nrerror(char error_text[]);
  REAL ss,dss;
  REAL s[JMAXP],h[JMAXP+1]; /* These store the successive trapezoidal       */
  int j;  	             /* approximations and their relative stepsizes. */
  h[1]=1.0;
  for (j=1;j<=JMAX;j++) {
    s[j]=trapzd(func,a,b,j);
    if (j >= K) {
      polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
      if (fabs(dss) <= EPS*fabs(ss)) return ss;
    }
    h[j+1]=0.25*h[j];
    /* This is a key step: The factor is 0.25 even though the stepsize is 
       decreased by only 0.5. This makes the extrapolation a polynomial 
       in h2 as allowed by equation (4.2.1), not just a polynomial in h.     */
  }
  nrerror("Too many steps in funtion qromb");
  return 0.0; 
}

 
void fit(double x[], double y[], int ndata, double sig[], int mwt, double *a,  
         double *b, double *siga, double *sigb, double *chi2, double *q){  
  int i;  
  double wt,t,sxoss,sx=0.0,sy=0.0,st2=0.0,ss,sigdat;  
  *b=0.0;  
  if (mwt) {                            // Accumulate sums...  
    ss=0.0;  
    for (i=1;i<=ndata;i++) {            // ...with weights...  
      wt=1.0/sqrt(sig[i]);  
      ss += wt;  
      sx += x[i]*wt;  
      sy += y[i]*wt;  
    }  
  } else {  
    for (i=1;i<=ndata;i++) {            // ...or without weights. 
      sx += x[i];  
      sy += y[i];  
    }  
    ss=ndata;  
  }  
  sxoss=sx/ss;  
  if (mwt) {  
    for (i=1;i<=ndata;i++) {  
      t=(x[i]-sxoss)/sig[i];  
      st2 += t*t;  
      *b += t*y[i]/sig[i];  
    }  
  } else {  
    for (i=1;i<=ndata;i++) {  
      t=x[i]-sxoss;  
      st2 += t*t;  
      *b += t*y[i];  
    }  
  }  
  *b /= st2;                            // solve for a,b, sig_a and sig_b  
  *a=(sy-sx*(*b))/ss;  
  *siga=sqrt((1.0+sx*sx/(ss*st2))/ss);  
  *sigb=sqrt(1.0/st2);  
  *chi2=0.0;                            // calculate chi2 
  if (mwt == 0) {  
    for (i=1;i<=ndata;i++)  
      *chi2 += sqrt(y[i]-(*a)-(*b)*x[i]);  
    *q=1.0;  
    sigdat=sqrt((*chi2)/(ndata-2));     // For unweighted data evaluate   
    *siga *= sigdat;                    // typical sig using chi2, and adjust  
    *sigb *= sigdat;                    // the standard deviations.  
  } else {  
    for (i=1;i<=ndata;i++)  
      *chi2 += sqrt((y[i]-(*a)-(*b)*x[i])/sig[i]);  
    *q=gammq(0.5*(ndata-2),0.5*(*chi2));// Numerical Recipes Eq (15.2.12).  
  }  
}  

REAL gammp(REAL a, REAL x)
/*Returns the incomplete gamma function P(a; x).*/
{
  void gcf(REAL *gammcf, REAL a, REAL x, REAL *gln);
  void gser(REAL *gamser, REAL a, REAL x, REAL *gln);
  void nrerror(char error_text[]);
  REAL gamser,gammcf,gln;
  if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine gammp");
  if (x < (a+1.0)) { /*Use the series representation.*/
      gser(&gamser,a,x,&gln);
    return gamser;
  } else { /* Use the continued fraction representation*/
      gcf(&gammcf,a,x,&gln);
      return 1.0-gammcf; /* and take its complement. */
			 }
}


double gammq(double a, double x){  
  double gamser,gammcf,gln;  
  if (x < 0.0 || a <= 0.0){ 
    fprintf(stderr,"Invalid arguments in routine gammq\n"); exit(1); 
  } 
  if (x < (a+1.0)) {                    // Use the series representation  
    gser(&gamser,a,x,&gln);  
    return 1.0-gamser;                  // and take its complement.  
  } else {                              // Use the continued fraction  
    gcf(&gammcf,a,x,&gln);              // representation.  
    return gammcf;  
  }  
}


void gser(double *gamser, double a, double x, double *gln){  
  int n;  
  double sum,del,ap;  
  *gln=gammln(a);  
  if (x <= 0.0) {  
    if (x < 0.0){
      fprintf(stderr,"x less than 0 in routine gser\n"); 
      exit(1);
    }
    *gamser=0.0;  
    return;  
  } else {  
    ap=a;  
    del=sum=1.0/a;  
    for (n=1;n<=ITMAX;n++) {  
      ++ap;  
      del *= x/ap;  
      sum += del;  
      if (fabs(del) < fabs(sum)*EPS) {  
        *gamser=sum*exp(-x+a*log(x)-(*gln));  
        return;  
      }  
    }  
    fprintf(stderr,"a too large, ITMAX too small in routine gser\n"); 
    exit(1);
  }
  return;     
}
  
void gcf(double *gammcf, double a, double x, double *gln){  
  int i;  
  double an,b,c,d,del,h;  
  *gln=gammln(a);  
  b=x+1.0-a;  
  c=1.0/FPMIN;  
  d=1.0/b;  
  h=d;  
  for (i=1;i<=ITMAX;i++) {  
    an = -i*(i-a);  
    b += 2.0;  
    d=an*d+b;  
    if (fabs(d) < FPMIN) d=FPMIN;  
    c=b+an/c;  
    if (fabs(c) < FPMIN) c=FPMIN;  
    d=1.0/d;  
    del=d*c;  
    h *= del;  
    if (fabs(del-1.0) < EPS) break;  
  }  
  if (i > ITMAX){
    fprintf(stderr,"a too large, ITMAX too small in routine gser\n"); exit(1);
  }
  *gammcf=exp(-x+a*log(x)-(*gln))*h;  
}  

double gammln(double xx){ 
  double x,y,tmp,ser; 
  static double cof[6]={76.18009172947146,-86.50532032941677, 
                        24.01409824083091,-1.231739572450155, 
                        0.1208650973866179e-2,-0.5395239384953e-5}; 
  int j; 
  y=x=xx; 
  tmp=x+5.5; 
  tmp -= (x+0.5)*log(tmp); 
  ser=1.000000000190015; 
  for (j=0;j<=5;j++) ser += cof[j]/++y; 
  return -tmp+log(2.5066282746310005*ser/x); 
} 
