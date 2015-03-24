/* ============================================================ *
 * mathstuff.c							*
 * Martin Kilbinger, NR, 2004-2008				*
 * ============================================================ */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <sys/times.h>

#ifdef NR_COMPLEX_H_
//#include "nrcomplex.h"
#endif

#include "mathstuff.h"

#define NRANSI /* from polint */
#define NR_END 1
#define FREE_ARG char*

double Dsqr(double x)
{
   return x*x;
}

double abs22(const double x[], const double y[])
{
   return sqrt( Dsqr(x[0]-y[0]) + Dsqr(x[1]-y[1]) );
}


double abs2sqr_f(const double *x, const double *y)
{
   double d0, d1;

   d0 = x[0]-y[0];
   d1 = x[1]-y[1];
   return d0*d0 + d1*d1;
}

double absnsqr(const double *x, const double *y, int ndim)
{
   int k;
   double a;

   for (k=0,a=0.0; k<ndim; k++) {
      a += (x[k] - y[k]) * (x[k] - y[k]);
   }

   return a;
}

void vdiff(const real *x, const real *y, real *z)
{
   z[0] = x[0]-y[0];
   z[1] = x[1]-y[1];
}

void out_error(const char *s)
{
   (void)fprintf(stderr, "error: ");
   (void)fprintf(stderr, "%s", s);
   (void)fprintf(stderr, "\n");
   assert(0);
}

FILE *fileopen(const char name[], const char mode [])
{
   FILE *F;
   int l;
   char *s;

   if ((F = fopen(name, mode))!=NULL) {
#ifndef NO_FILE_MSG
      if (mode[0]=='r') {
	 fprintf(stderr, "reading %s\n", name);
      } else {
	 fprintf(stderr, "writing %s\n", name);
      }
#endif
      return F;
   } else {
      l = 50 + strlen(name);
      s = (char*)malloc(l*sizeof(char));
      sprintf(s, "Error while opening file %s", name);
      perror(s);
      exit(errno);
   }
}

void fileclose(FILE *F)
{
   if (!fclose(F)) {
      return;
   } else {
      perror("error while closing file");
   }
}

/* ============================================================ *
 * Interpolates f at the value x, where f is a real[n]		*
 * representing a function between a and b, stepwidth dx.	*
 * 'lower' and 'upper' are powers of a logarithmic power law	*
 * extrapolation. If no	extrapolation desired, set these to     *
 * > 1e30. For constant upper extrapolation, set upper=0.	*
 * ============================================================ */
real interpol(real *f, int n, real a, real b, real dx, real x,
	      real lower, real upper)
{
   real r;
   int  i;
   if (x < a) {
      assert(lower!=0.);
      return f[0] + lower*(x - a);
   }
   r = (x - a)/dx;
   i = (int)(floor(r));
   if (i+1 >= n) {
      if (upper==0.0) {
	 if (i+1 == n) {
	    return f[i];  /* constant extrapolation */
	 }
	 assert(0);
      } else {
	 assert(upper<1.e30);
	 return f[n-1] + upper*(x-b); /* linear extrapolation */
      }
   } else {
      return (r - i)*(f[i+1] - f[i]) + f[i]; /* interpolation */
   }

   assert(0); return 0;
}


/* ============================================================ *
 * like interpol, but f beeing a 2d-function			*
 * 'lower' and 'upper' are the powers of a power law extra-	*
 * polation in the second argument				*
 * ============================================================ */

real interpol2d(real **f, int nx, real ax, real bx, real dx, real x,
		int ny, real ay, real by, real dy, real y, real lower, real upper)
{
   real t, dt, s, ds;
   int i, j;

   assert(x>=ax);
   assert(x<=bx);

   t = (x - ax)/dx;
   i = (int)(floor(t));
   assert(i+1<nx); assert(i>=0);

   dt = t - i;
   if (y < ay) {
      assert(lower<1.e30);
      return ((1.-dt)*f[i][0] + dt*f[i+1][0]) + (y-ay)*lower;
   } else if (y > by) {
      assert(upper<1.e30);
      return ((1.-dt)*f[i][ny-1] + dt*f[i+1][ny-1]) + (y-by)*upper;
   }
   s = (y - ay)/dy;
   j = (int)(floor(s));
   ds = s - j;
   return (1.-dt)*(1.-ds)*f[i][j] + (1.-dt)*ds*f[i][j+1] +
     dt*(1.-ds)*f[i+1][j] + dt*ds*f[i+1][j+1];
}


/* ============================================================ *
 * nrutil.c							*
 * ============================================================ */

void nrerror(const char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

real *vector(long nl, long nh)
/* allocate a real vector with subscript range v[nl..nh] */
{
	real *v;

	v=(real *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(real)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

void free_vector(real *v, long nl, long nh)
/* free a real vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

real **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a real matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	real **m;

	/* allocate pointers to rows */
	m=(real **) malloc((size_t)((nrow+NR_END)*sizeof(real*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(real *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(real)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate an integer matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	/* allocate pointers to rows */
	m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	if (!m) nrerror("allocation failure 1 in imatrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}


void free_matrix(real **m, long nrl, long nrh, long ncl, long nch)
/* free a real matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an ineteger matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}


/* ============================================================ *
 * polint.c							*
 * Interpolation						*
 * ============================================================ */

void polint(real xa[], real ya[], int n, real x, real *y, real *dy)
{
	int i,m,ns=1;
	real den,dif,dift,ho,hp,w;
	real *c,*d;

	dif=fabs(x-xa[1]);
	c=vector(1,n);
	d=vector(1,n);
	for (i=1;i<=n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0)
			  nrerror("Error in routine polint");
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	free_vector(d,1,n);
	free_vector(c,1,n);
}


/* ============================================================ *
 * trapzd.c							*
 * Returns approximated integral over func between a and b.	*
 * It computes the nth stage of refinement of an extended	*
 * trapezoidal rule.						*
 * Must be called succesively with n = 1,2,3...			*
 * NR p.137							*
 * ============================================================ */

#define FUNC(x,y) ((*func)(x,y))

real trapzd(real (*func)(real, const real*), const real *intpar, real a, real b, int n)
{
	real x,tnm,sum,del;
	static real s;
	int it,j;

	if (n == 1) {
		return (s=0.5*(b-a)*(FUNC(a,intpar)+FUNC(b,intpar)));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) {
		   sum += FUNC(x,intpar);
		}
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}
#undef FUNC

#define FUNC(x,y) ((*func)(x,y))

/* ============================================================ *
 * Modified trapzd, static variable s now argument.		*
 * ============================================================ */

real trapzdberg(real (*func)(real, const real*), const real *intpar, real a, real b, 
		int n, double *s)
{
	real x,tnm,sum,del;
	int it,j;

	if (n == 1) {
		return (*s=0.5*(b-a)*(FUNC(a,intpar)+FUNC(b,intpar)));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) {
		   sum += FUNC(x,intpar);
		}
		*s=0.5*(*s+(b-a)*sum/tnm);
		return *s;
	}
}
#undef FUNC


/* ============================================================= *
 * trapzd2d. Computes the nth stage of refinement of an extended *
 * 2d trapezian rule. func is a pointer to a function with two   *
 * real variables, it is to be integrated			 *
 * between ax<=x<=ay and ay<=y<=by. Works analogous to the one-  *
 * dimensional trapzd (Numerical Recipes, p.137).		 *
 * ============================================================= */

#define FUNC(x,y) ((*func)(x,y))

real trapzd2d(real (*func)(real, real), real ax, real bx, real ay, real by, int n)
{
   real x, y, tnm, sum, delx, dely;
   static real s;
   int it, j, k;

   if (n==1) {
      /* 1. stage: Function is evaluated only at the corners. */
      return (s=0.25*(bx-ax)*(by-ay)*(FUNC(ax,ay)+FUNC(bx,ay)+FUNC(ax,by)+FUNC(bx,by)));
   } else {
      for (it=1, j=1; j<n-1; j++) it <<= 1;
      tnm = it;
      /* delx, dely are the spacings of the points to be added */
      delx = (bx-ax)/tnm;
      dely = (by-ay)/tnm;
      x = ax + 0.5*delx;
      /* Evaluates function at the borders of the integration area */
      for (sum=0.0, j=1; j<=it; j++, x+=delx) {
	 sum += 0.5*(FUNC(x,ay) + FUNC(x,by));
      }
      y = ay + 0.5*dely;
      for (j=1; j<=it; j++, y+=dely) {
	 sum += 0.5*(FUNC(ax,y) + FUNC(bx,y));
      }
      /* Interiour of the area */
      x = ax + 0.5*delx;
      for (j=1; j<=it; j++, x+= delx) {
	 y = ay + 0.5*dely;
	 for (k=1; k<=it; k++, y+=dely) {
	    sum += FUNC(x,y);
	 }
      }
      x = ax + delx;
      for (j=1; j<it; j++, x+=delx) {
	 y = ay + 0.5*dely;
	 for (k=1; k<=it; k++, y+=dely) {
	    sum += FUNC(x,y);
	 }
      }
      x = ax + 0.5*delx;
      for (j=1; j<=it; j++, x+=delx) {
	 y = ay + dely;
	 for (k=1; k<it; k++, y+=dely) {
	    sum += FUNC(x,y);
	 }
      }
      /* Replace s by its refined value */
      s = 0.25*(s + (bx-ax)*(by-ay)*sum/(tnm*tnm));
      return s;
   }
}

#undef FUNC


#define FUNC(x,y) ((*func)(x,y))

real trapzd2(real (*func)(real, const real *), const real *intpar, real a, real b, int n)
{
	real x,tnm,sum,del;
	static real s;
	int it,j;

	if (n == 1) {
		return (s=0.5*(b-a)*(FUNC(a,intpar)+FUNC(b,intpar)));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) {
		   sum += FUNC(x,intpar);
		}
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}
#undef FUNC


#define FUNC(x,y) ((*func)(x,y))

real trapzd3(real (*func)(real, const real*), const real *intpar, real a, real b, int n)
{
	real x,tnm,sum,del;
	static real s;
	int it,j;

	if (n == 1) {
		return (s=0.5*(b-a)*(FUNC(a, intpar)+FUNC(b, intpar)));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) {
		   sum += FUNC(x, intpar);
		}
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}
#undef FUNC


#define FUNC(x,y) ((*func)(x,y))

real trapzd1(real (*func)(real, const real *), const real *intpar, real a, real b, int n)
{
	real x,tnm,sum,del;
	static real s;
	int it,j;

	if (n == 1) {
		return (s=0.5*(b-a)*(FUNC(a,intpar)+FUNC(b,intpar)));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) {
		   sum += FUNC(x,intpar);
		}
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}
#undef FUNC


#define JMAX 35
#define JMAXP (JMAX+1)
#define K 5

real qromb3(real (*func)(real, const real*), const real *intpar, real a, real b, real EPS)
{
	real ss,dss;
	real s[JMAXP],h[JMAXP+1];
	int j;

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=trapzd3(func,intpar,a,b,j);
		if (j >= K) {
			polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) <= EPS*fabs(ss)) return ss;
		}
		h[j+1]=0.25*h[j];
	}
	nrerror("Too many steps in routine qromb3");
	return 0.0;
}
#undef JMAX
#undef JMAXP
#undef K


#define JMAX 35
#define JMAXP (JMAX+1)
#define K 5

real qromb2(real (*func)(real, const real*), const real *intpar, real a, real b, real EPS)
{
	real ss,dss;
	real s[JMAXP],h[JMAXP+1];
	int j;

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=trapzd2(func,intpar,a,b,j);
		if (j >= K) {
			polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) <= EPS*fabs(ss)) {
			   return ss;
			}
		}
		h[j+1]=0.25*h[j];
	}
	nrerror("Too many steps in routine qromb2");
	return 0.0;
}
#undef EPS
#undef JMAX
#undef JMAXP
#undef K


#define JMAX 40
#define JMAXP (JMAX+1)
#define K 5

real qromb1(real (*func)(real, const real *), const real* intpar, real a, real b, real EPS)
{
	real ss,dss;
	real s[JMAXP],h[JMAXP+1];
	int j;

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=trapzd1(func,intpar,a,b,j);
		if (j >= K) {
			polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) <= EPS*fabs(ss)) return ss;
		}
		h[j+1]=0.25*h[j];
	}
	nrerror("Too many steps in routine qromb1");
	return 0.0;
}
#undef EPS
#undef JMAX
#undef JMAXP
#undef K

/* ============================================================ *
 * Modified from qromb.c, NR p. 140.			        *
 * Now with additional parameters as argument (intpar) and      *
 * static variable defined here, given as argument strap to     *
 * trapzd.							*
 * ============================================================ */

#define JMAX 35
#define JMAXP (JMAX+1)
#define K 5

real qromberg(real (*func)(real, const real*), const real *intpar, real a, real b, real EPS)
{
	real ss,dss;
	real s[JMAXP],h[JMAXP+1], strap=0.0;
	int j;

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=trapzdberg(func,intpar,a,b,j,&strap);
		if (j >= K) {
			polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) <= EPS*fabs(ss)) return ss;
		}
		h[j+1]=0.25*h[j];
	}
	nrerror("Too many steps in routine qromb");
	return 0.0;
}
#undef JMAX
#undef JMAXP
#undef K


/* ============================================================ *
 * qromb.c							*
 * Romberg Integration. Uses trapzd. NR p. 140			*
 * ============================================================ */

#define JMAX 35
#define JMAXP (JMAX+1)
#define K 5

real qromb(real (*func)(real, const real*), const real *intpar, real a, real b, real EPS)
{
	real ss,dss;
	real s[JMAXP],h[JMAXP+1];
	int j;

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=trapzd(func,intpar,a,b,j);
		if (j >= K) {
			polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) <= EPS*fabs(ss)) return ss;
		}
		h[j+1]=0.25*h[j];
	}
	nrerror("Too many steps in routine qromb");
	return 0.0;
}
/* #undef EPS */
#undef JMAX
#undef JMAXP
#undef K

#define JMAX 35
#define JMAXP (JMAX+1)
#define K 5

/* ============================================================ *
 * Romberg integration of an open intervall [a,b]. choose is a  *
 * pointer to a routine using an open quadrature formula.	*
 * Uses polint (Neville-Aitken) to extrapolate. NR p. 143	*
 * ============================================================ */

real qromo(real (*func)(real, const real*), const real *intpar, real a, real b,
	real (*choose)(real(*)(real, const real*), const real *, real, real, int), real EPS)
{
	int j;
	real ss,dss,h[JMAXP+1],s[JMAXP];

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=(*choose)(func,intpar,a,b,j);
		if (j >= K) {
			polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) <= EPS*fabs(ss)) return ss;
		}
		h[j+1]=h[j]/9.0;
	}
	nrerror("Too many steps in routing qromo");
	return 0.0;
}
#undef JMAX
#undef JMAXP

#define JMAX 35
#define JMAXP (JMAX+1)
#define K 5

real qrombergo(real (*func)(real, const real*), const real *intpar, real a, real b,
	real (*choose)(real(*)(real, const real*), const real *, real, real, int, real *),
	       real EPS)
{
	int j;
	real ss,dss,h[JMAXP+1],s[JMAXP], schoose;

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=(*choose)(func,intpar,a,b,j,&schoose);
		if (j >= K) {
			polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) <= EPS*fabs(ss)) return ss;
		}
		h[j+1]=h[j]/9.0;
	}
	nrerror("Too many steps in routing qromo");
	return 0.0;
}
#undef JMAX
#undef JMAXP

#define FUNC(x,y) ((*func)(x,y))

real midpntberg(real (*func)(real, const real*), const real *intpar, real a, real b, 
		int n, real *s)
{
	real x,tnm,sum,del,ddel;
	int it,j;

	if (n == 1) {
		return (*s=(b-a)*FUNC(0.5*(a+b), intpar));
	} else {
		for(it=1,j=1;j<n-1;j++) it *= 3;
		tnm=it;
		del=(b-a)/(3.0*tnm);
		ddel=del+del;
		x=a+0.5*del;
		sum=0.0;
		for (j=1;j<=it;j++) {
			sum += FUNC(x, intpar);
			x += ddel;
			sum += FUNC(x, intpar);
			x += del;
		}
		*s = (*s+(b-a)*sum/tnm)/3.0;
		return *s;
	}
}
#undef FUNC

#define FUNC(x,y) ((*func)(x,y))

real midpnt(real (*func)(real, const real*), const real *intpar, real a, real b, int n)
{
	real x,tnm,sum,del,ddel;
	static real s;
	int it,j;

	if (n == 1) {
		return (s=(b-a)*FUNC(0.5*(a+b), intpar));
	} else {
		for(it=1,j=1;j<n-1;j++) it *= 3;
		tnm=it;
		del=(b-a)/(3.0*tnm);
		ddel=del+del;
		x=a+0.5*del;
		sum=0.0;
		for (j=1;j<=it;j++) {
			sum += FUNC(x, intpar);
			x += ddel;
			sum += FUNC(x, intpar);
			x += del;
		}
		s=(s+(b-a)*sum/tnm)/3.0;
		return s;
	}
}
#undef FUNC


#define FUNC(x,y) ((*funk)(1.0/(x),y)/((x)*(x)))

real midinf(real (*funk)(real, const real *), const real *intpar, real aa, real bb, int n)
{
	real x,tnm,sum,del,ddel,b,a;
	static real s;
	int it,j;

	b=1.0/aa;
	a=1.0/bb;
	if (n == 1) {
		return (s=(b-a)*FUNC(0.5*(a+b), intpar));
	} else {
		for(it=1,j=1;j<n-1;j++) it *= 3;
		tnm=it;
		del=(b-a)/(3.0*tnm);
		ddel=del+del;
		x=a+0.5*del;
		sum=0.0;
		for (j=1;j<=it;j++) {
			sum += FUNC(x, intpar);
			x += ddel;
			sum += FUNC(x, intpar);
			x += del;
		}
		return (s=(s+(b-a)*sum/tnm)/3.0);
	}
}
#undef FUNC


#define FUNC(x,y) (2.0*(x)*(*funk)(aa+(x)*(x),y))

real midsql(real (*funk)(real, const real*), const real *intpar, real aa, real bb, int n)
{
	real x,tnm,sum,del,ddel,a,b;
	static real s;
	int it,j;

	b=sqrt(bb-aa);
	a=0.0;
	if (n == 1) {
		return (s=(b-a)*FUNC(0.5*(a+b),intpar));
	} else {
		for(it=1,j=1;j<n-1;j++) it *= 3;
		tnm=it;
		del=(b-a)/(3.0*tnm);
		ddel=del+del;
		x=a+0.5*del;
		sum=0.0;
		for (j=1;j<=it;j++) {
			sum += FUNC(x,intpar);
			x += ddel;
			sum += FUNC(x,intpar);
			x += del;
		}
		s=(s+(b-a)*sum/tnm)/3.0;
		return s;
	}
}

#undef FUNC


#define FUNC(x,y) (2.0*(x)*(*funk)(bb-(x)*(x),y))

real midsqu(real (*funk)(real, const real*), const real *intpar, real aa, real bb, int n)
{
	real x,tnm,sum,del,ddel,a,b;
	static real s;
	int it,j;

	b=sqrt(bb-aa);
	a=0.0;
	if (n == 1) {
		return (s=(b-a)*FUNC(0.5*(a+b),intpar));
	} else {
		for(it=1,j=1;j<n-1;j++) it *= 3;
		tnm=it;
		del=(b-a)/(3.0*tnm);
		ddel=del+del;
		x=a+0.5*del;
		sum=0.0;
		for (j=1;j<=it;j++) {
			sum += FUNC(x,intpar);
			x += ddel;
			sum += FUNC(x,intpar);
			x += del;
		}
		s=(s+(b-a)*sum/tnm)/3.0;
		return s;
	}
}

#undef FUNC


/* ============================================================ *
 * sobseq.c							*
 * NR p.312. Sub-random number generator.			*
 * n < 0: initializes a set of MAXBIT direction numbers for	*
 * each of MAXDIM different Sobol' sequences.			*
 * n > 0: returns as x[1..n] the next values from n of these	*
 * sequences.							*
 * ============================================================ */

#define MAXBIT 30
#define MAXDIM 6

void sobseq(int *n, real x[])
{
	int k,l;
	unsigned long i,im,ipp,j;     /* changed: j was (signed) int! */
	static real fac;
	static unsigned long in,ix[MAXDIM+1],*iu[MAXBIT+1];
	static unsigned long mdeg[MAXDIM+1]={0,1,2,3,3,4,4};
	static unsigned long ip[MAXDIM+1]={0,0,1,1,2,1,4};
	static unsigned long iv[MAXDIM*MAXBIT+1]={
		0,1,1,1,1,1,1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9};

	if (*n < 0) {
		for (k=1;k<=MAXDIM;k++) ix[k]=0;
		in=0;
		if (iv[1] != 1) return;
		fac=1.0/(1L << MAXBIT);
		for (j=1,k=0;j<=MAXBIT;j++,k+=MAXDIM) iu[j] = &iv[k];
		for (k=1;k<=MAXDIM;k++) {
			for (j=1;j<=mdeg[k];j++) iu[j][k] <<= (MAXBIT-j);
			for (j=mdeg[k]+1;j<=MAXBIT;j++) {
				ipp=ip[k];
				i=iu[j-mdeg[k]][k];
				i ^= (i >> mdeg[k]);
				for (l=mdeg[k]-1;l>=1;l--) {
					if (ipp & 1) i ^= iu[j-l][k];
					ipp >>= 1;
				}
				iu[j][k]=i;
			}
		}
	} else {
		im=in++;
		for (j=1;j<=MAXBIT;j++) {
			if (!(im & 1)) break;
			im >>= 1;
		}
		if (j > MAXBIT) nrerror("MAXBIT too small in sobseq");
		im=(j-1)*MAXDIM;
		for (k=1;k<=IMIN(*n,MAXDIM);k++) {
			ix[k] ^= iv[im+k];
			x[k]=ix[k]*fac;
		}
	}
}
#undef MAXBIT
#undef MAXDIM

void random_init(char name[])
{
   struct tms time;
   clock_t cl;
   unsigned u;
   FILE *F;

   /* random init */
   cl = times(&time);
   u = (unsigned)cl;
   u = (u%10)*(u%100)*(u%1000);
   srand(u);
   if (name!=NULL) {
      F = fileopen(name, "w");
      fprintf(F, "%u\n", u);
      fileclose(F);
   }
}

/* ============================================================ *
 * dfridr.c							*
 * NR page 188. Returns derivate of func at x, initial step is  *
 * h. Error estimate in err.					*
 * Modified! func depends on two real! (like P_L)		*
 * ============================================================ */

#define CON 1.4
#define CON2 (CON*CON)
#define BIG 1.0e30
#define NTAB 10
#define SAFE 2.0


real dfridr(real (*func)(real,real), real x, real h, real *err, real aa)
{
	int i,j;
	real errt,fac,hh,**a,ans;

	ans = 1e30; /* dummy initialization */
	if (h == 0.0) nrerror("h must be nonzero in dfridr.");
	a=matrix(1,NTAB,1,NTAB);
	hh=h;
	a[1][1]=((*func)(aa,x+hh)-(*func)(aa,x-hh))/(2.0*hh);
	*err=BIG;
	for (i=2;i<=NTAB;i++) {
		hh /= CON;
		a[1][i]=((*func)(aa,x+hh)-(*func)(aa,x-hh))/(2.0*hh);
		fac=CON2;
		for (j=2;j<=i;j++) {
			a[j][i]=(a[j-1][i]*fac-a[j-1][i-1])/(fac-1.0);
			fac=CON2*fac;
			errt=FMAX(fabs(a[j][i]-a[j-1][i]),fabs(a[j][i]-a[j-1][i-1]));
			if (errt <= *err) {
				*err=errt;
				ans=a[j][i];
			}
		}
		if (fabs(a[i][i]-a[i-1][i-1]) >= SAFE*(*err)) break;
	}
	free_matrix(a,1,NTAB,1,NTAB);
	return ans;
}
#undef CON
#undef CON2
#undef BIG
#undef NTAB
#undef SAFE


/* ============================================================ *
 * dfridr.c							*
 * NR page 188. Returns derivate of func at x, initial step is  *
 * h. Error estimate in err.					*
 * Original version!						*
 * ============================================================ */

#define CON 1.4
#define CON2 (CON*CON)
#define BIG 1.0e30
#define NTAB 10
#define SAFE 2.0


real dfridr1(real (*func)(real), real x, real h, real *err)
{
	int i,j;
	real errt,fac,hh,**a,ans;

	ans = 1e30; /* dummy initialization */
	if (h == 0.0) nrerror("h must be nonzero in dfridr.");
	a=matrix(1,NTAB,1,NTAB);
	hh=h;
	a[1][1]=((*func)(x+hh)-(*func)(x-hh))/(2.0*hh);
	*err=BIG;
	for (i=2;i<=NTAB;i++) {
		hh /= CON;
		a[1][i]=((*func)(x+hh)-(*func)(x-hh))/(2.0*hh);
		fac=CON2;
		for (j=2;j<=i;j++) {
			a[j][i]=(a[j-1][i]*fac-a[j-1][i-1])/(fac-1.0);
			fac=CON2*fac;
			errt=FMAX(fabs(a[j][i]-a[j-1][i]),fabs(a[j][i]-a[j-1][i-1]));
			if (errt <= *err) {
				*err=errt;
				ans=a[j][i];
			}
		}
		if (fabs(a[i][i]-a[i-1][i-1]) >= SAFE*(*err)) break;
	}
	free_matrix(a,1,NTAB,1,NTAB);
	return ans;
}
#undef CON
#undef CON2
#undef BIG
#undef NTAB
#undef SAFE


void spline(real x[], real y[], int n, real yp1, real ypn, real y2[])
{
	int i,k;
	real p,qn,sig,un,*u;

	u=vector(1,n-1);
	if (yp1 > 0.99e30)
		y2[1]=u[1]=0.0;
	else {
		y2[1] = -0.5;
		u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
	}
	for (i=2;i<=n-1;i++) {
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	if (ypn > 0.99e30)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
	}
	y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
	for (k=n-1;k>=1;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
	free_vector(u,1,n-1);
}


void splint(real xa[], real ya[], real y2a[], int n, real x, real *y)
{
	int klo,khi,k;
	real h,b,a;

	klo=1;
	khi=n;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	if (h == 0.0) nrerror("Bad xa input to routine splint");
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	*y=a*ya[klo]+b*ya[khi]+
	  ((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}


/* ============================================================ *
 * Splines with equidistant x values				*
 * ============================================================ */

void spli_pol(real *y, int n, real dx, real yp1, real ypn, real *y2)
{
   int i;
   real *u, p, qn, un;

   u = vector(1, n-1);
   if (yp1 > .99e30) {       /* zero first derivative */
      y2[1] = u[1] = 0.;
   } else {
      y2[1] = -0.5;
      u[1]  = 3.*(y[2]-y[1])/DSQR(dx)-yp1;
   }
   for (i=2; i<n; i++) {
      p = 0.5*y2[i-1] + 2.;
      y2[i] = -0.5/p;
      u[i] = (3.*(y[i+1]-2.*y[i]+y[i-1])/DSQR(dx)-0.5*u[i-1])/p;
   }
   if (ypn > .99e-30) {
      qn = un = 0.;
   } else {
      qn = 0.5;
      un = 3./dx*(ypn-(y[n]-y[n-1])/dx);
   }
   y2[n] = (un-qn*u[n-1])/(qn*y2[n-1]+1.);
   for (i=n-1; i>=1; i--) {
      y2[i] = y2[i]*y2[i+1]+u[i];
   }
   free_vector(u, 1, n-1);
}


real spli_value(real x, real *y, real *y2, int n, real xmin, real dx)
{
   real rn, a, b, spl;
   int nn;

   rn = (x-xmin)/dx+1.;
   nn = (int)rn;
   if (nn>n-1 || nn<1) {
      fprintf(stderr, "%f requested (index %d, max is %d)\n", x, nn, n);
      out_error("extrapolation with splines is not a good idea!");
   }
   b  = rn-(real)nn;
   a = 1. - b;
   spl = a*y[nn] + b*y[nn+1] + ((a*a*a-a)*y2[nn] + (b*b*b-b)*y2[nn+1])*dx*dx/6.;
   return spl;
}


/* ============================================================ *
 * Special functions.						*
 * ============================================================ */

real ms_sinc(real x)
{
   if (x<0.001 && x>-0.001) return 1. - x*x/6.;
   else return sin(x)/x;
}

real ms_gammln(real xx)
{
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
	return (real)(-tmp+log(2.5066282746310005*ser/x));
}

real ms_bessj0(real x)
{
	real ax,z;
	double xx,y,ans,ans1,ans2;

	if ((ax=fabs(x)) < 8.0) {
		y=x*x;
		ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7
			+y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
		ans2=57568490411.0+y*(1029532985.0+y*(9494680.718
			+y*(59272.64853+y*(267.8532712+y*1.0))));
		ans=ans1/ans2;
	} else {
		z=8.0/ax;
		y=z*z;
		xx=ax-0.785398164;
		ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
			+y*(-0.2073370639e-5+y*0.2093887211e-6)));
		ans2 = -0.1562499995e-1+y*(0.1430488765e-3
			+y*(-0.6911147651e-5+y*(0.7621095161e-6
			-y*0.934935152e-7)));
		ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
	}
	return ans;
}

real bessj1(real x)
{
	real ax,z;
	double xx,y,ans,ans1,ans2;

	if ((ax=fabs(x)) < 8.0) {
		y=x*x;
		ans1=x*(72362614232.0+y*(-7895059235.0+y*(242396853.1
			+y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
		ans2=144725228442.0+y*(2300535178.0+y*(18583304.74
			+y*(99447.43394+y*(376.9991397+y*1.0))));
		ans=ans1/ans2;
	} else {
		z=8.0/ax;
		y=z*z;
		xx=ax-2.356194491;
		ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4
			+y*(0.2457520174e-5+y*(-0.240337019e-6))));
		ans2=0.04687499995+y*(-0.2002690873e-3
			+y*(0.8449199096e-5+y*(-0.88228987e-6
			+y*0.105787412e-6)));
		ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
		if (x < 0.0) ans = -ans;
	}
	return ans;
}


#define ACC 40.0
#define BIGNO 1.0e10
#define BIGNI 1.0e-10


real bessj(int n, real x)
{
	int j,jsum,m;
	real ax,bj,bjm,bjp,sum,tox,ans;

	if (n < 2) nrerror("Index n less than 2 in bessj");
	ax=fabs(x);
	if (ax == 0.0)
		return 0.0;
	else if (ax > (real) n) {
		tox=2.0/ax;
		bjm=ms_bessj0(ax);
		bj=bessj1(ax);
		for (j=1;j<n;j++) {
			bjp=j*tox*bj-bjm;
			bjm=bj;
			bj=bjp;
		}
		ans=bj;
	} else {
		tox=2.0/ax;
		m=2*((n+(int) sqrt(ACC*n))/2);
		jsum=0;
		bjp=ans=sum=0.0;
		bj=1.0;
		for (j=m;j>0;j--) {
			bjm=j*tox*bj-bjp;
			bjp=bj;
			bj=bjm;
			if (fabs(bj) > BIGNO) {
				bj *= BIGNI;
				bjp *= BIGNI;
				ans *= BIGNI;
				sum *= BIGNI;
			}
			if (jsum) sum += bj;
			jsum=!jsum;
			if (j == n) ans=bjp;
		}
		sum=2.0*sum-bj;
		ans /= sum;
	}
	return x < 0.0 && (n & 1) ? -ans : ans;
}
#undef ACC
#undef BIGNO
#undef BIGNI

real BesselJ(int n, real x)
{
   double res;

   if (n==0) res = ms_bessj0(x);
   else if (n==1 || n==-1) res = bessj1(x);
   else res = bessj(n, x);

   if (n<0) res = -res;
   return res;
}

real bessy0(real x)
{
	real z;
	double xx,y,ans,ans1,ans2;

	if (x < 8.0) {
		y=x*x;
		ans1 = -2957821389.0+y*(7062834065.0+y*(-512359803.6
			+y*(10879881.29+y*(-86327.92757+y*228.4622733))));
		ans2=40076544269.0+y*(745249964.8+y*(7189466.438
			+y*(47447.26470+y*(226.1030244+y*1.0))));
		ans=(ans1/ans2)+0.636619772*ms_bessj0(x)*log(x);
	} else {
		z=8.0/x;
		y=z*z;
		xx=x-0.785398164;
		ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
			+y*(-0.2073370639e-5+y*0.2093887211e-6)));
		ans2 = -0.1562499995e-1+y*(0.1430488765e-3
			+y*(-0.6911147651e-5+y*(0.7621095161e-6
			+y*(-0.934945152e-7))));
		ans=sqrt(0.636619772/x)*(sin(xx)*ans1+z*cos(xx)*ans2);
	}
	return ans;
}

real bessy1(real x)
{
	real z;
	double xx,y,ans,ans1,ans2;

	if (x < 8.0) {
		y=x*x;
		ans1=x*(-0.4900604943e13+y*(0.1275274390e13
			+y*(-0.5153438139e11+y*(0.7349264551e9
			+y*(-0.4237922726e7+y*0.8511937935e4)))));
		ans2=0.2499580570e14+y*(0.4244419664e12
			+y*(0.3733650367e10+y*(0.2245904002e8
			+y*(0.1020426050e6+y*(0.3549632885e3+y)))));
		ans=(ans1/ans2)+0.636619772*(bessj1(x)*log(x)-1.0/x);
	} else {
		z=8.0/x;
		y=z*z;
		xx=x-2.356194491;
		ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4
			+y*(0.2457520174e-5+y*(-0.240337019e-6))));
		ans2=0.04687499995+y*(-0.2002690873e-3
			+y*(0.8449199096e-5+y*(-0.88228987e-6
			+y*0.105787412e-6)));
		ans=sqrt(0.636619772/x)*(sin(xx)*ans1+z*cos(xx)*ans2);
	}
	return ans;
}

real bessi0(real x)
{
	real ax,ans;
	double y;

	if ((ax=fabs(x)) < 3.75) {
		y=x/3.75;
		y*=y;
		ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
			+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
	} else {
		y=3.75/ax;
		ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
			+y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
			+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
			+y*0.392377e-2))))))));
	}
	return ans;
}


real bessi1(real x)
{
	real ax,ans;
	double y;

	if ((ax=fabs(x)) < 3.75) {
		y=x/3.75;
		y*=y;
		ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
			+y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
	} else {
		y=3.75/ax;
		ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
			-y*0.420059e-2));
		ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
			+y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
		ans *= (exp(ax)/sqrt(ax));
	}
	return x < 0.0 ? -ans : ans;
}


#define ACC 40.0
#define BIGNO 1.0e10
#define BIGNI 1.0e-10

real bessi(int n, real x)
{
	int j;
	real bi,bim,bip,tox,ans;

	if (n < 2) nrerror("Index n less than 2 in bessi");
	if (x == 0.0)
		return 0.0;
	else {
		tox=2.0/fabs(x);
		bip=ans=0.0;
		bi=1.0;
		for (j=2*(n+(int) sqrt(ACC*n));j>0;j--) {
			bim=bip+j*tox*bi;
			bip=bi;
			bi=bim;
			if (fabs(bi) > BIGNO) {
				ans *= BIGNI;
				bi *= BIGNI;
				bip *= BIGNI;
			}
			if (j == n) ans=bip;
		}
		ans *= bessi0(x)/bi;
		return x < 0.0 && (n & 1) ? -ans : ans;
	}
}
#undef ACC
#undef BIGNO
#undef BIGNI

/* Cholesky decomposition A=L*L^T, NR p.97 *
 * Original version, starts counting at 1  */
void choldc(real **a, int n, real p[])
{
   int i,j,k;
   real sum;

   for (i=1;i<=n;i++) {
      for (j=i;j<=n;j++) {
	 for (sum=a[i][j],k=i-1;k>=1;k--) {
	    sum -= a[i][k]*a[j][k];
	 }
	 if (i == j) {
	    if (sum <= 0.0) {
	       fprintf(stderr, "sum=%.2e (diag=%.2e) for column %d\n",
		       sum, a[i][i], i);
	       nrerror("choldc failed");
	    }
	    p[i]=sqrt(sum);
	 } else a[j][i]=sum/p[i];
      }
   }
}

/* ============================================================ *
 * Returns L where C = L L^t. C is destroyed.			*
 * ============================================================ */

real **cholesky(real **C, int N)
{
   real **L, *diag;
   int i, j;

   L    = matrix(1, N, 1, N);
   diag = vector(1, N);

   choldc(C, N, diag);

   for(i=1; i<=N; i++) {
      for (j=1; j<i; j++) {
	 L[i][j] = C[i][j];     /* write lower triangle */
	 L[j][i] = 0.0;
      }
      L[i][i] = diag[i];
   }

   free_vector(diag, 1, N);

   return L;
}

#define TINY 1.0e-20;

void ludcmp(real **a, int n, int *indx, real *d)
{
	int i,imax=0,j,k;
	real big,dum,sum,temp;
	real *vv;

	vv=vector(1,n);
	*d=1.0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
		vv[i]=1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a[i][j];
			for (k=1;k<i;k++)
			  sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a[i][j];
			for (k=1;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) a[i][j] *= dum;
		}
	}
	free_vector(vv,1,n);
}
#undef TINY


void lubksb(real **a, int n, int *indx, real b[])
{
	int i,ii=0,ip,j;
	real sum;

	for (i=1;i<=n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n;i>=1;i--) {
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}


/* sigma*gasdev + m is a gaussian random variable with mean m and variance sigma */
real gasdev()
{
   static int iset=0;
   static real gset;
   real fac,rsq,v1,v2;

   if  (iset == 0) {
      do {
	 v1=2.0*rand()/RAND_MAX-1.0;
	 v2=2.0*rand()/RAND_MAX-1.0;
	 rsq=v1*v1+v2*v2;
      } while (rsq >= 1.0 || rsq == 0.0);
      fac=sqrt(-2.0*log(rsq)/rsq);
      gset=v1*fac;
      iset=1;
      return v2*fac;
   } else {
      iset=0;
      return gset;
   }
}


real ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a real 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
   long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
   real ***t;

   /* allocate pointers to pointers to rows */
   t=(real ***) malloc((size_t)((nrow+NR_END)*sizeof(real**)));
   if (!t) nrerror("allocation failure 1 in f3tensor()");
   t += NR_END;
   t -= nrl;

   /* allocate pointers to rows and set pointers to them */
   t[nrl]=(real **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(real*)));
   if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
   t[nrl] += NR_END;
   t[nrl] -= ncl;

   /* allocate rows and set pointers to them */
   t[nrl][ncl]=(real *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(real)));
   if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
   t[nrl][ncl] += NR_END;
   t[nrl][ncl] -= ndl;

   for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
   for(i=nrl+1;i<=nrh;i++) {
      t[i]=t[i-1]+ncol;
      t[i][ncl]=t[i-1][ncl]+ncol*ndep;
      for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
   }

   /* return pointer to array of pointers to rows */
   return t;
}

void free_f3tensor(real ***t, long nrl, long nrh, long ncl, long nch,
        long ndl, long ndh)
/* free a real f3tensor allocated by f3tensor() */
{
        free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
        free((FREE_ARG) (t[nrl]+ncl-NR_END));
        free((FREE_ARG) (t+nrl-NR_END));
}


real factln(int n)
{
   static real a[101];

   if (n < 0) {
      nrerror("Negative factorial in routine factln");
   }
   if (n <= 1) return 0.0;
   if (n <= 100) return a[n] ? a[n] : (a[n]=ms_gammln(n+1.0));
   else return ms_gammln(n+1.0);
}


real fact(int n)
{
   return exp(factln(n));
}


real bico(int n, int k)
{
   return floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
}


/* modified from nr: now with a! */

#define MAXIT 30
real rtsec(real (*func)(real,real), real x1, real x2, real xacc, real a)
{
	int j;
	real fl,f,dx,swap,xl,rts;

	fl=(*func)(a, x1);
	f=(*func)(a, x2);
	if (fabs(fl) < fabs(f)) {
		rts=x1;
		xl=x2;
		swap=fl;
		fl=f;
		f=swap;
	} else {
		xl=x1;
		rts=x2;
	}
	for (j=1;j<=MAXIT;j++) {
		dx=(xl-rts)*f/(f-fl);
		xl=rts;
		fl=f;
		rts += dx;
		f=(*func)(a, rts);
		if (fabs(dx) < xacc || f == 0.0) return rts;
	}
	nrerror("Maximum number of iterations exceeded in rtsec");
	return 0.0;
}
#undef MAXIT


#define JMAX 40
real rtbis(real (*func)(real,real), real x1, real x2, real xacc, real a)
{
	int j;
	real dx,f,fmid,xmid,rtb;

	f=(*func)(a, x1);
	fmid=(*func)(a, x2);
	if (f*fmid >= 0.0) {
	   nrerror("Root must be bracketed for bisection in rtbis");
	}
	rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
	for (j=1;j<=JMAX;j++) {
		fmid=(*func)(a, xmid=rtb+(dx *= 0.5));
		if (fmid <= 0.0) rtb=xmid;
		if (fabs(dx) < xacc || fmid == 0.0) return rtb;
	}
	nrerror("Too many bisections in rtbis");
	return 0.0;
}
#undef JMAX


#define JMAX 40
real rtbis1(real (*func)(real), real x1, real x2, real xacc)
{
   int j;
   real dx,f,fmid,xmid,rtb;

   f=(*func)(x1);
   fmid=(*func)(x2);
   if (f*fmid >= 0.0) {
      nrerror("Root must be bracketed for bisection in rtbis");
   }
   rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
   for (j=1;j<=JMAX;j++) {
      fmid=(*func)(xmid=rtb+(dx *= 0.5));
      if (fmid <= 0.0) rtb=xmid;
      if (fabs(dx) < xacc || fmid == 0.0) return rtb;
   }
   nrerror("Too many bisections in rtbis");
   return 0.0;
}
#undef JMAX

/* returns the minimum of the n-dim vector x */

real minvec(real *x, int n)
{
   real m;
   int i;

   m = 1.e30;
   for (i=0; i<n; i++) {
      if (x[i]<m) m = x[i];
   }

   return m;
}


/* returns the maximum of the n-dim vector x */

real maxvec(real *x, int n)
{
   real m;
   int i;

   m = -1.e30;
   for (i=0; i<n; i++) {
      if (x[i]>m) m = x[i];
   }

   return m;
}

/* cyclically permutes the vector x[i] i=offset..offset+3 */
void permute(real *x, int offset)
{
   int i;
   real tmp;

   for (tmp=x[2+offset],i=2+offset; i>offset; i--) {
      x[i] = x[i-1];
   }
   x[offset] = tmp;
}

#define NRANSI

void lfit(real x[], real y[], real sig[], int ndat, real a[], int ia[],
	int ma, real **covar, real *chisq, void (*funcs)(real, real [], int))
{
	int i,j,k,l,m,mfit=0;
	real ym,wt,sum,sig2i,**beta,*afunc;

	beta=matrix(1,ma,1,1);
	afunc=vector(1,ma);
	for (j=1;j<=ma;j++)
		if (ia[j]) mfit++;
	if (mfit == 0) nrerror("lfit: no parameters to be fitted");
	for (j=1;j<=mfit;j++) {
		for (k=1;k<=mfit;k++) covar[j][k]=0.0;
		beta[j][1]=0.0;
	}
	for (i=1;i<=ndat;i++) {
		(*funcs)(x[i],afunc,ma);
		ym=y[i];
		if (mfit < ma) {
			for (j=1;j<=ma;j++)
				if (!ia[j]) ym -= a[j]*afunc[j];
		}
		sig2i=1.0/DSQR(sig[i]);
		for (j=0,l=1;l<=ma;l++) {
			if (ia[l]) {
				wt=afunc[l]*sig2i;
				for (j++,k=0,m=1;m<=l;m++)
					if (ia[m]) covar[j][++k] += wt*afunc[m];
				beta[j][1] += ym*wt;
			}
		}
	}
	for (j=2;j<=mfit;j++)
		for (k=1;k<j;k++)
			covar[k][j]=covar[j][k];
	gaussj(covar,mfit,beta,1);
	for (j=0,l=1;l<=ma;l++)
		if (ia[l]) a[l]=beta[++j][1];
	*chisq=0.0;
	for (i=1;i<=ndat;i++) {
		(*funcs)(x[i],afunc,ma);
		for (sum=0.0,j=1;j<=ma;j++) sum += a[j]*afunc[j];
		*chisq += DSQR((y[i]-sum)/sig[i]);
	}
	covsrt(covar,ma,ia,mfit);
	free_vector(afunc,1,ma);
	free_matrix(beta,1,ma,1,1);
}
#undef NRANSI


#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

void covsrt(real **covar, int ma, int ia[], int mfit)
{
	int i,j,k;
	real swap;

	for (i=mfit+1;i<=ma;i++)
		for (j=1;j<=i;j++) covar[i][j]=covar[j][i]=0.0;
	k=mfit;
	for (j=ma;j>=1;j--) {
		if (ia[j]) {
			for (i=1;i<=ma;i++) SWAP(covar[i][k],covar[i][j])
			for (i=1;i<=ma;i++) SWAP(covar[k][i],covar[j][i])
			k--;
		}
	}
}

#undef SWAP


#define NRANSI
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

void gaussj(real **a, int n, real **b, int m)
{
	int *indxc,*indxr,*ipiv;
	int i,icol=0,irow=0,j,k,l,ll;
	real big,dum,pivinv,temp;

	indxc=ivector(1,n);
	indxr=ivector(1,n);
	ipiv=ivector(1,n);
	for (j=1;j<=n;j++) ipiv[j]=0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if (ipiv[j] != 1)
				for (k=1;k<=n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					} else if (ipiv[k] > 1) nrerror("gaussj: Singular Matrix-1");
				}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l])
			for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l])
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0) nrerror("gaussj: Singular Matrix-2");
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=1;l<=n;l++) a[icol][l] *= pivinv;
		for (l=1;l<=m;l++) b[icol][l] *= pivinv;
		for (ll=1;ll<=n;ll++)
			if (ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
				for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
			}
	}
	for (l=n;l>=1;l--) {
		if (indxr[l] != indxc[l])
			for (k=1;k<=n;k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}
	free_ivector(ipiv,1,n);
	free_ivector(indxr,1,n);
	free_ivector(indxc,1,n);
}
#undef SWAP


#define NRANSI
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);

void jacobi(real **a, int n, real d[], real **v, int *nrot)
{
   int j,iq,ip,i;
   real tresh,theta,tau,t,sm,s,h,g,c,*b,*z;

   b=vector(1,n);
   z=vector(1,n);
   for (ip=1;ip<=n;ip++) {
      for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
      v[ip][ip]=1.0;
   }
   for (ip=1;ip<=n;ip++) {
      b[ip]=d[ip]=a[ip][ip];
      z[ip]=0.0;
   }
   *nrot=0;
   for (i=1;i<=50;i++) {
      sm=0.0;
      for (ip=1;ip<=n-1;ip++) {
	 for (iq=ip+1;iq<=n;iq++)
	   sm += fabs(a[ip][iq]);
      }
      if (sm == 0.0) {
	 free_vector(z,1,n);
	 free_vector(b,1,n);
	 return;
      }
      if (i < 4)
	tresh=0.2*sm/(n*n);
      else
	tresh=0.0;
      for (ip=1;ip<=n-1;ip++) {
	 for (iq=ip+1;iq<=n;iq++) {
	    g=100.0*fabs(a[ip][iq]);
	    if (i > 4 && (real)(fabs(d[ip])+g) == (real)fabs(d[ip])
		&& (real)(fabs(d[iq])+g) == (real)fabs(d[iq]))
	      a[ip][iq]=0.0;
	    else if (fabs(a[ip][iq]) > tresh) {
	       h=d[iq]-d[ip];
	       if ((real)(fabs(h)+g) == (real)fabs(h))
		 t=(a[ip][iq])/h;
	       else {
		  theta=0.5*h/(a[ip][iq]);
		  t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
		  if (theta < 0.0) t = -t;
	       }
	       c=1.0/sqrt(1+t*t);
	       s=t*c;
	       tau=s/(1.0+c);
	       h=t*a[ip][iq];
	       z[ip] -= h;
	       z[iq] += h;
	       d[ip] -= h;
	       d[iq] += h;
	       a[ip][iq]=0.0;
	       for (j=1;j<=ip-1;j++) {
		  ROTATE(a,j,ip,j,iq)
		    }
	       for (j=ip+1;j<=iq-1;j++) {
		  ROTATE(a,ip,j,j,iq)
		    }
	       for (j=iq+1;j<=n;j++) {
		  ROTATE(a,ip,j,iq,j)
		    }
	       for (j=1;j<=n;j++) {
		  ROTATE(v,j,ip,j,iq)
		    }
	       ++(*nrot);
	    }
	 }
      }
      for (ip=1;ip<=n;ip++) {
	 b[ip] += z[ip];
	 d[ip]=b[ip];
	 z[ip]=0.0;
      }
   }
   nrerror("Too many iterations in routine jacobi");
}
#undef ROTATE

#define RADIX 2.0

void balanc(double **a, int n)
{
   int last,j,i;
   double s,r,g,f,c,sqrdx;

   sqrdx=RADIX*RADIX;
   last=0;
   while (last == 0) {
      last=1;
      for (i=1;i<=n;i++) {
	 r=c=0.0;
	 for (j=1;j<=n;j++)
	   if (j != i) {
	      c += fabs(a[j][i]);
	      r += fabs(a[i][j]);
	   }
	 if (c && r) {
	    g=r/RADIX;
	    f=1.0;
	    s=c+r;
	    while (c<g) {
	       f *= RADIX;
	       c *= sqrdx;
	    }
	    g=r*RADIX;
	    while (c>g) {
	       f /= RADIX;
	       c /= sqrdx;
	    }
	    if ((c+r)/f < 0.95*s) {
	       last=0;
	       g=1.0/f;
	       for (j=1;j<=n;j++) a[i][j] *= g;
	       for (j=1;j<=n;j++) a[j][i] *= f;
	    }
	 }
      }
   }
}

#define SWAP(g,h) {y=(g);(g)=(h);(h)=y;}

void elmhes(double **a, int n)
{
   int m,j,i;
   double y,x;

   for (m=2;m<n;m++) {
      x=0.0;
      i=m;
      for (j=m;j<=n;j++) {
	 if (fabs(a[j][m-1]) > fabs(x)) {
	    x=a[j][m-1];
	    i=j;
	 }
      }
      if (i != m) {
	 for (j=m-1;j<=n;j++) SWAP(a[i][j],a[m][j])
				for (j=1;j<=n;j++) SWAP(a[j][i],a[j][m])
						     }
      if (x) {
	 for (i=m+1;i<=n;i++) {
	    if ((y=a[i][m-1]) != 0.0) {
	       y /= x;
	       a[i][m-1]=y;
	       for (j=m;j<=n;j++)
		 a[i][j] -= y*a[m][j];
	       for (j=1;j<=n;j++)
		 a[j][m] += y*a[j][i];
	    }
	 }
      }
   }
}

#undef SWAP

#define NRANSI

void hqr(double **a, int n, double wr[], double wi[])
{
   int nn,m,l,k,j,its,i,mmin;
   double z,y,x,w,v,u,t,s,r=0.0,q=0.0,p=0.0,anorm;

   anorm=0.0;
   for (i=1;i<=n;i++)
     for (j=IMAX(i-1,1);j<=n;j++)
       anorm += fabs(a[i][j]);
   nn=n;
   t=0.0;
   while (nn >= 1) {
      its=0;
      do {
	 for (l=nn;l>=2;l--) {
	    s=fabs(a[l-1][l-1])+fabs(a[l][l]);
	    if (s == 0.0) s=anorm;
	    if ((double)(fabs(a[l][l-1]) + s) == s) break;
	 }
	 x=a[nn][nn];
	 if (l == nn) {
	    wr[nn]=x+t;
	    wi[nn--]=0.0;
	 } else {
	    y=a[nn-1][nn-1];
	    w=a[nn][nn-1]*a[nn-1][nn];
	    if (l == (nn-1)) {
	       p=0.5*(y-x);
	       q=p*p+w;
	       z=sqrt(fabs(q));
	       x += t;
	       if (q >= 0.0) {
		  z=p+SIGN(z,p);
		  wr[nn-1]=wr[nn]=x+z;
		  if (z) wr[nn]=x-w/z;
		  wi[nn-1]=wi[nn]=0.0;
	       } else {
		  wr[nn-1]=wr[nn]=x+p;
		  wi[nn-1]= -(wi[nn]=z);
	       }
	       nn -= 2;
	    } else {
	       if (its == 30) nrerror("Too many iterations in hqr");
	       if (its == 10 || its == 20) {
		  t += x;
		  for (i=1;i<=nn;i++) a[i][i] -= x;
		  s=fabs(a[nn][nn-1])+fabs(a[nn-1][nn-2]);
		  y=x=0.75*s;
		  w = -0.4375*s*s;
	       }
	       ++its;
	       for (m=(nn-2);m>=l;m--) {
		  z=a[m][m];
		  r=x-z;
		  s=y-z;
		  p=(r*s-w)/a[m+1][m]+a[m][m+1];
		  q=a[m+1][m+1]-z-r-s;
		  r=a[m+2][m+1];
		  s=fabs(p)+fabs(q)+fabs(r);
		  p /= s;
		  q /= s;
		  r /= s;
		  if (m == l) break;
		  u=fabs(a[m][m-1])*(fabs(q)+fabs(r));
		  v=fabs(p)*(fabs(a[m-1][m-1])+fabs(z)
			     +fabs(a[m+1][m+1]));
		  if ((double)(u+v) == v) break;
	       }
	       for (i=m+2;i<=nn;i++) {
		  a[i][i-2]=0.0;
		  if (i != (m+2)) a[i][i-3]=0.0;
	       }
	       for (k=m;k<=nn-1;k++) {
		  if (k != m) {
		     p=a[k][k-1];
		     q=a[k+1][k-1];
		     r=0.0;
		     if (k != (nn-1)) r=a[k+2][k-1];
		     if ((x=fabs(p)+fabs(q)+fabs(r)) != 0.0) {
			p /= x;
			q /= x;
			r /= x;
		     }
		  }
		  if ((s=SIGN(sqrt(p*p+q*q+r*r),p)) != 0.0) {
		     if (k == m) {
			if (l != m)
			  a[k][k-1] = -a[k][k-1];
		     } else
		       a[k][k-1] = -s*x;
		     p += s;
		     x=p/s;
		     y=q/s;
		     z=r/s;
		     q /= p;
		     r /= p;
		     for (j=k;j<=nn;j++) {
			p=a[k][j]+q*a[k+1][j];
			if (k != (nn-1)) {
			   p += r*a[k+2][j];
			   a[k+2][j] -= p*z;
			}
			a[k+1][j] -= p*y;
			a[k][j] -= p*x;
		     }
		     mmin = nn<k+3 ? nn : k+3;
		     for (i=l;i<=mmin;i++) {
			p=x*a[i][k]+y*a[i][k+1];
			if (k != (nn-1)) {
			   p += z*a[i][k+2];
			   a[i][k+2] -= p*r;
			}
			a[i][k+1] -= p*q;
			a[i][k] -= p;
		     }
		  }
	       }
	    }
	 }
      } while (l < nn-1);
   }
}

#undef NRANSI

/* ============================================================ *
 * Returns a*b.							*
 * ============================================================ */

real **matrix_multiply(const real **a, const real **b, int N)
{
   real **c, sum;
   int i, j, k;

   c = matrix(1, N, 1, N);
   for (i=1; i<=N; i++) {
      for (j=1; j<=N; j++) {
	 for (k=1,sum=0.0; k<=N; k++) {
	    sum += a[i][k]*b[k][j];
	 }
	 c[i][j] = sum;
      }
   }

   return c;
}

/* ============================================================ *
 * Returns b = A*x.						*
 * ============================================================ */

real *matrix_vector_multiply(const real **A, const real *x, int N)
{
   real *b;
   int i, j;

   b = vector(1, N);
   for (i=1,b[i]=0.0; i<=N; i++) {
      for (j=1; j<=N; j++) {
	 b[i] += A[i][j]*x[j];
      }
   }
   return b;
}

real **transpose(const real **a, int N)
{
   real **at;
   int i, j;

   at = matrix(1, N, 1, N);
   for (i=1; i<=N; i++) {
      for (j=1; j<=N; j++) {
	 at[i][j] = a[j][i];
      }
   }
   return at;
}

/* ============================================================ *
 * (See NR p.46 - 49.)						*
 * Inverts the NxN matrix C and stores the inverse in C         *
 * (old matrix gets destroyed). Returns det C.			*
 * The index goes from 1 to N !					*
 * ============================================================ */

real inverse(real **C, int N)
{
   int i, j;
   real *col, det, **Cinv;
   int *p;

   Cinv = (real**)malloc((N+1)*sizeof(real));
   for (i=1; i<=N; i++) {
      Cinv[i]   = (real*)malloc((N+1)*sizeof(real));
   }
   col = (real*)malloc((N+1)*sizeof(real));
   p = (int*)malloc((N+1)*sizeof(int));
   ludcmp(C, N, p, &det);
   for (j=1; j<=N; j++) {
      det *= C[j][j];
   }
   for (j=1; j<=N; j++) {
      for (i=1; i<=N; i++) {
	 col[i] = 0.0;
      }
      col[j] = 1.0;
      lubksb(C, N, p, col);
      for (i=1; i<=N; i++) {
	 Cinv[i][j] = col[i];
      }
   }
   for (j=1; j<=N; j++) {
      for (i=1; i<=N; i++) {
	 C[i][j] = Cinv[i][j];
      }
   }
   for (i=1; i<=N; i++) {
      free((void*)(Cinv[i]));
   }
   free((void*)Cinv);
   return det;
}

/* allocates memory for a new NxN-matrix and copies M to it */
real **matrix_copy(const real **M, int N)
{
   int i, j;
   real **C;

   C = matrix(1, N, 1, N);
   for (i=1; i<=N; i++) {
      for (j=1; j<=N; j++) {
	 C[i][j] = M[i][j];
      }
   }
   return C;
}

void unity_test(const real **a, const real **ainv, int N)
{
   int i, j;
   real max1, max2;

   real **unit;

   unit = matrix_multiply(a, ainv, N);

   for (i=1,max1=-1.,max2=-1.; i<=N; i++) {
      if (fabs(unit[i][i]-1.0)>max1) max1 = fabs(unit[i][i]-1.0);
      for (j=1; j<=N; j++) {
	 if (i!=j && fabs(unit[i][j])>max2) max2 = fabs(unit[i][j]);
      }
   }

   fprintf(stderr, "max |diag-1| = %.2e, max |offdiag| = %.2e\n", max1, max2);
}

void out_matrix(char *name, const real **a, int n)
{
   int i, j;
   FILE *F;

   F = fileopen(name, "w");
   for (i=1; i<=n; i++) {
      for (j=1; j<=n; j++) {
	 fprintf(F, "% .3e ", a[i][j]);
      }
      fprintf(F, "\n");
   }
   fileclose(F);
}

#define NRANSI
#define EPS 1.0e-14

void linbcg(unsigned long n, const double **A, const double **Aguess,
	    const double b[], double x[], int itol, double tol, int itmax,
	    int *iter, double *err, int guess_is_diag)
{
	unsigned long j;
	double ak,akden,bk,bkden=0.0,bknum,bnrm=0.0,dxnrm,xnrm,zm1nrm,znrm=0.0;
	double *p,*pp,*r,*rr,*z,*zz;

	p=vector(1,n);
	pp=vector(1,n);
	r=vector(1,n);
	rr=vector(1,n);
	z=vector(1,n);
	zz=vector(1,n);

	*iter=0;
	atimes((const double **)A, (const double *)x, r, n, 0);
	for (j=1;j<=n;j++) {
		r[j]=b[j]-r[j];
		rr[j]=r[j];
	}
	if (itol == 1) {
		bnrm=snrm(n,(const double *)b,itol);
		asolve((const double **)Aguess, (const double *)r, z,n, 0, guess_is_diag);
	}
	else if (itol == 2) {
		asolve((const double **)Aguess, (const double *)b, z, n, 0, guess_is_diag);
		bnrm=snrm(n,z,itol);
		asolve((const double **)Aguess, (const double *)r, z, n, 0, guess_is_diag);
	}
	else if (itol == 3 || itol == 4) {
		asolve((const double **)Aguess, (const double *)b, z, n, 0, guess_is_diag);
		bnrm=snrm(n,z,itol);
		asolve((const double **)Aguess, (const double *)r, z, n, 0, guess_is_diag);
		znrm=snrm(n,z,itol);
	} else nrerror("illegal itol in linbcg");
	while (*iter <= itmax) {
		++(*iter);
		asolve((const double **)Aguess, (const double *)rr, zz, n, 1, guess_is_diag);
		for (bknum=0.0,j=1;j<=n;j++) bknum += z[j]*rr[j];
		if (*iter == 1) {
			for (j=1;j<=n;j++) {
				p[j]=z[j];
				pp[j]=zz[j];
			}
		}
		else {
			bk=bknum/bkden;
			for (j=1;j<=n;j++) {
				p[j]=bk*p[j]+z[j];
				pp[j]=bk*pp[j]+zz[j];
			}
		}
		bkden=bknum;
		atimes((const double**)A, (const double *)p, z, n, 0);
		for (akden=0.0,j=1;j<=n;j++) akden += z[j]*pp[j];
		ak=bknum/akden;
		atimes((const double**)A, (const double *)pp, zz,n, 1);
		for (j=1;j<=n;j++) {
			x[j] += ak*p[j];
			r[j] -= ak*z[j];
			rr[j] -= ak*zz[j];
		}
		asolve((const double **)Aguess, (const double *)r, z, n, 0, guess_is_diag);
		if (itol == 1)
			*err=snrm(n,r,itol)/bnrm;
 		else if (itol == 2)
			*err=snrm(n,z,itol)/bnrm;
		else if (itol == 3 || itol == 4) {
			zm1nrm=znrm;
			znrm=snrm(n,z,itol);
			if (fabs(zm1nrm-znrm) > EPS*znrm) {
				dxnrm=fabs(ak)*snrm(n,p,itol);
				*err=znrm/fabs(zm1nrm-znrm)*dxnrm;
			} else {
				*err=znrm/bnrm;
				continue;
			}
			xnrm=snrm(n,x,itol);
			if (*err <= 0.5*xnrm) *err /= xnrm;
			else {
				*err=znrm/bnrm;
				continue;
			}
		}
		/* printf("iter=%4d err=%12.6f\n",*iter,*err); */
	if (*err <= tol) break;
	}

	free_vector(p,1,n);
	free_vector(pp,1,n);
	free_vector(r,1,n);
	free_vector(rr,1,n);
	free_vector(z,1,n);
	free_vector(zz,1,n);
}
#undef EPS
#undef NRANSI


/* ============================================================ *
 * Vector norm of sx. itol<=3: oo-norm, itol=4: 2-norm.         *
 * ============================================================ */

double snrm(unsigned long n, const double sx[], int itol)
{
   unsigned long i,isamax;
   double ans;

   if (itol <= 3) {
      ans = 0.0;
      for (i=1;i<=n;i++) ans += sx[i]*sx[i];
      return sqrt(ans);
   } else {
      isamax=1;
      for (i=1;i<=n;i++) {
	 if (fabs(sx[i]) > fabs(sx[isamax])) isamax=i;
      }
      return fabs(sx[isamax]);
   }
}
/* ============================================================ *
 * Calculates r = A*x (itrnsp=0) or r=A^t*x (itrnsp=1).		*
 * Used by linbcg.						*
 * ============================================================ */

void atimes(const double** A, const double x[], double r[],
	    unsigned long n, int itrnsp)
{
   int i, k;

   if (itrnsp==0) {
      for (i=1; i<=n; i++) {
	 for (k=1, r[i]=0; k<=n; k++) {
	    r[i] += A[i][k]*x[k];
	 }
      }
   } else {
      for (i=1; i<=n; i++) {
	 for (k=1, r[i]=0; k<=n; k++) {
	    r[i] += A[k][i]*x[k];
	 }
      }
   }

}

/* ============================================================ *
 * Solves Aguess*x = b for some preconditioner matrix Aguess.   *
 * Used by linbcg.						*
 * ============================================================ */

void asolve(const double **Aguess, const double b[], double x[], unsigned long n,
	    int itrnsp, int guess_is_diag)
{
   int i;

   double **A, d;
   int *indx;

   if (guess_is_diag==0) {
      indx = ivector(1, n);
      A = matrix_copy(Aguess, n);
      for (i=1; i<=n; i++) x[i] = b[i];

      ludcmp(A, n, indx, &d);
      lubksb(A, n, indx, x);

      free_matrix(A, 1, n, 1, n);
      free_ivector(indx, 1, n);
   } else {
      for (i=1; i<=n; i++) x[i] = (Aguess[i][i]!=0.0 ? b[i]/Aguess[i][i] : b[i]);
   }

}

#ifdef NR_COMPLEX_H_

dcomplex **cmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a complex matrix with subscript range m[nrl..nrh][ncl..nch] */
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        dcomplex **m;

        /* allocate pointers to rows */
        m=(dcomplex **) malloc((size_t)((nrow+NR_END)*sizeof(dcomplex*)));
        if (!m) nrerror("allocation failure 1 in matrix()");
        m += NR_END;
        m -= nrl;

        /* allocate rows and set pointers to them */
        m[nrl]=(dcomplex *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(dcomplex)));
        if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
        m[nrl] += NR_END;
        m[nrl] -= ncl;

        for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

        /* return pointer to array of pointers to rows */
        return m;
}

dcomplex ***c3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a dcomplex 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
   long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
   dcomplex ***t;

   /* allocate pointers to pointers to rows */
   t=(dcomplex ***) malloc((size_t)((nrow+NR_END)*sizeof(dcomplex**)));
   if (!t) nrerror("allocation failure 1 in f3tensor()");
   t += NR_END;
   t -= nrl;

   /* allocate pointers to rows and set pointers to them */
   t[nrl]=(dcomplex **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(dcomplex*)));
   if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
   t[nrl] += NR_END;
   t[nrl] -= ncl;

   /* allocate rows and set pointers to them */
   t[nrl][ncl]=(dcomplex *) malloc((size_t)((nrow*ncol*ndep+NR_END)
					    *sizeof(dcomplex)));
   if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
   t[nrl][ncl] += NR_END;
   t[nrl][ncl] -= ndl;

   for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
   for(i=nrl+1;i<=nrh;i++) {
      t[i]=t[i-1]+ncol;
      t[i][ncl]=t[i-1][ncl]+ncol*ndep;
      for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
   }

   /* return pointer to array of pointers to rows */
   return t;
}

void free_cmatrix(dcomplex **m, long nrl, long nrh, long ncl, long nch)
/* free a complex matrix allocated by cmatrix() */
{
        free((FREE_ARG) (m[nrl]+ncl-NR_END));
        free((FREE_ARG) (m+nrl-NR_END));
}


void free_c3tensor(dcomplex ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh)
/* free a dcomplex c3tensor allocated by c3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}

void cisitab(real x, real *ci, real *si)
{
#define Ncisi 50000
   const real logcisimin = -10.;
   const real logcisimax =  10.;
   static real tabci[Ncisi], tabsi[Ncisi],
     tabci2[Ncisi], tabsi2[Ncisi], dlogeta, c1;
   /* static real tablogeta[Ncisi]; */
   static int flagcisi = 0;
   real logeta, dummy;
   int i;

   if (flagcisi==0) {
      dlogeta = (logcisimax-logcisimin)/(Ncisi-1.);
      for (i=0,logeta=logcisimin; i<Ncisi; i++,logeta+=dlogeta) {
	 /* tablogeta[i] = logeta; */
	 cisi(exp(logeta), tabci+i, tabsi+i);
      }
      cisi(exp(logcisimin), &c1, &dummy);
      /* cisi(exp(-21.), &c2, &dummy); */

      /* (c2-c1)/(-21.-logcisimin) */
      /* spline(tablogeta-1, tabci-1, Ncisi, 1.0, 0.0, tabci2-1); */
      spli_pol(tabci-1, Ncisi, dlogeta, 1.0, 0.0, tabci2-1);
      spli_pol(tabsi-1, Ncisi, dlogeta, 0.0, 0.0, tabsi2-1);
      /* spline(tablogeta-1, tabsi-1, Ncisi, 0.0, 0.0, tabsi2-1); */
      flagcisi = 1;
   }

   logeta = log(x);
   if (logeta<logcisimin) {
      /* *ci = (logeta-logcisimin)/(-21.-logcisimin)*(c2-c1) + c1; */
      *ci = (logeta-logcisimin) + c1;
      *si = 0.0;
   } else if (logeta>logcisimax) {
      *ci = 0.0;
      *si = pi/2.;
   } else {
      /* *ci = interpol(tabci, Ncisi, logcisimin, logcisimax, dlogeta, logeta, 0., 0.); */
      /* *si = interpol(tabsi, Ncisi, logcisimin, logcisimax, dlogeta, logeta, 0., 0.); */

      /* splint(tablogeta-1, tabci-1, tabci2-1, Ncisi, logeta, ci); */
      /* splint(tablogeta-1, tabsi-1, tabsi2-1, Ncisi, logeta, si); */
      *ci = spli_value(logeta, tabci-1, tabci2-1, Ncisi, logcisimin, dlogeta);
      *si = spli_value(logeta, tabsi-1, tabsi2-1, Ncisi, logcisimin, dlogeta);
   }
}


#define EPS 6.0e-8
#define EULER 0.57721566
#define MAXIT 100
#define PIBY2 1.5707963
#define FPMIN 1.0e-30
#define TMIN 2.0
#define TRUE 1
#define ONE Complex(1.0,0.0)
void cisi(real x, real *ci, real *si)
{
	int i,k,odd;
	real a,err,fact,sign,sum,sumc,sums,t,term;
	dcomplex h,b,c,d,del;

	t=fabs(x);
	if (t == 0.0) {
		*si=0.0;
		*ci = -1.0/FPMIN;
		return;
	}
	if (t > TMIN) {
		b=Complex(1.0,t);
		c=Complex(1.0/FPMIN,0.0);
		d=h=Cdiv(ONE,b);
		for (i=2;i<=MAXIT;i++) {
			a = -(i-1)*(i-1);
			b=Cadd(b,Complex(2.0,0.0));
			d=Cdiv(ONE,Cadd(RCmul(a,d),b));
			c=Cadd(b,Cdiv(Complex(a,0.0),c));
			del=Cmul(c,d);
			h=Cmul(h,del);
			if (fabs(del.r-1.0)+fabs(del.i) < EPS) break;
		}
		if (i > MAXIT) nrerror("cf failed in cisi");
		h=Cmul(Complex(cos(t),-sin(t)),h);
		*ci = -h.r;
		*si=PIBY2+h.i;
	} else {
		if (t < sqrt(FPMIN)) {
			sumc=0.0;
			sums=t;
		} else {
			sum=sums=sumc=0.0;
			sign=fact=1.0;
			odd=TRUE;
			for (k=1;k<=MAXIT;k++) {
				fact *= t/k;
				term=fact/k;
				sum += sign*term;
				err=term/fabs(sum);
				if (odd) {
					sign = -sign;
					sums=sum;
					sum=sumc;
				} else {
					sumc=sum;
					sum=sums;
				}
				if (err < EPS) break;
				odd=!odd;
			}
			if (k > MAXIT) nrerror("maxits exceeded in cisi");
		}
		*si=sums;
		*ci=sumc+log(t)+EULER;
	}
	if (x < 0.0) *si = -(*si);
}
#undef EPS
#undef EULER
#undef MAXIT
#undef PIBY2
#undef FPMIN
#undef TMIN
#undef TRUE
#undef ONE

#endif

#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
real brent(real ax, real bx, real cx, real (*f)(real), real tol,
        real *xmin)
{
        int iter;
        real a,b,d=0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
        real e=0.0;

        a=(ax < cx ? ax : cx);
        b=(ax > cx ? ax : cx);
        x=w=v=bx;
        fw=fv=fx=(*f)(x);
        for (iter=1;iter<=ITMAX;iter++) {
                xm=0.5*(a+b);
                tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
                if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
                        *xmin=x;
                        return fx;
                }
                if (fabs(e) > tol1) {
                        r=(x-w)*(fx-fv);
                        q=(x-v)*(fx-fw);
                        p=(x-v)*q-(x-w)*r;
                        q=2.0*(q-r);
                        if (q > 0.0) p = -p;
                        q=fabs(q);
                        etemp=e;
                        e=d;
                        if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
                                d=CGOLD*(e=(x >= xm ? a-x : b-x));
                        else {
                                d=p/q;
                                u=x+d;
                                if (u-a < tol2 || b-u < tol2)
                                        d=SIGN(tol1,xm-x);
                        }
                } else {
		   d=CGOLD*(e=(x >= xm ? a-x : b-x));
                }
                u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
                fu=(*f)(u);
                if (fu <= fx) {
                        if (u >= x) a=x; else b=x;
                        SHFT(v,w,x,u)
                        SHFT(fv,fw,fx,fu)
                } else {
                        if (u < x) a=u; else b=u;
                        if (fu <= fw || w == x) {
                                v=w;
                                w=u;
                                fv=fw;
                                fw=fu;
                        } else if (fu <= fv || v == x || v == w) {
                                v=u;
                                fv=fu;
                        }
                }
        }
        nrerror("Too many iterations in brent");
        *xmin=x;
        return fx;
}
#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT

void avevar(real data[], unsigned long n, real *ave, real *var)
{
        unsigned long j;
        real s,ep;

        for (*ave=0.0,j=1;j<=n;j++) *ave += data[j];
        *ave /= n;
        *var=ep=0.0;
        for (j=1;j<=n;j++) {
                s=data[j]-(*ave);
                ep += s;
                *var += s*s;
        }
        *var=(*var-ep*ep/n)/(n-1);
}

#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30
void gcf(real *gammcf, real a, real x, real *gln)
{
        int i;
        real an,b,c,d,del,h;

        *gln=ms_gammln(a);
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
        if (i > ITMAX) nrerror("a too large, ITMAX too small in gcf");
        *gammcf=exp(-x+a*log(x)-(*gln))*h;
}
#undef ITMAX
#undef EPS
#undef FPMIN

#define ITMAX 100
#define EPS 3.0e-7
void gser(real *gamser, real a, real x, real *gln)
{
        int n;
        real sum,del,ap;

        *gln=ms_gammln(a);
        if (x <= 0.0) {
                if (x < 0.0) nrerror("x less than 0 in routine gser");
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
                nrerror("a too large, ITMAX too small in routine gser");
                return;
        }
}
#undef ITMAX
#undef EPS

#define BIG 1.0e30
real gammq(real a, real x)
{
   //   void gcf(real *gammcf, real a, real x, real *gln);
   //   void gser(real *gamser, real a, real x, real *gln);
   real gamser,gammcf,gln;

        if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine gammq");
        if (x < (a+1.0)) {
                gser(&gamser,a,x,&gln);
                return 1.0-gamser;
        } else {
                gcf(&gammcf,a,x,&gln);
                return gammcf;
        }
}
#undef BIG

#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
void mnbrak(real *ax, real *bx, real *cx, real *fa, real *fb, real *fc,
	    real (*func)(real))
{
        real ulim,u,r,q,fu,dum;

        *fa=(*func)(*ax);
        *fb=(*func)(*bx);
        if (*fb > *fa) {
                SHFT(dum,*ax,*bx,dum)
                SHFT(dum,*fb,*fa,dum)
        }
        *cx=(*bx)+GOLD*(*bx-*ax);
        *fc=(*func)(*cx);
        while (*fb > *fc) {
                r=(*bx-*ax)*(*fb-*fc);
                q=(*bx-*cx)*(*fb-*fa);
                u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
                        (2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
                ulim=(*bx)+GLIMIT*(*cx-*bx);
                if ((*bx-u)*(u-*cx) > 0.0) {
                        fu=(*func)(u);
                        if (fu < *fc) {
                                *ax=(*bx);
                                *bx=u;
                                *fa=(*fb);
                                *fb=fu;
                                return;
                        } else if (fu > *fb) {
                                *cx=u;
                                *fc=fu;
                                return;
                        }
                        u=(*cx)+GOLD*(*cx-*bx);
                        fu=(*func)(u);
                } else if ((*cx-u)*(u-ulim) > 0.0) {
                        fu=(*func)(u);
                        if (fu < *fc) {
                                SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
                                SHFT(*fb,*fc,fu,(*func)(u))
                        }
                } else if ((u-ulim)*(ulim-*cx) >= 0.0) {
                        u=ulim;
                        fu=(*func)(u);
                } else {
                        u=(*cx)+GOLD*(*cx-*bx);
                        fu=(*func)(u);
                }
                SHFT(*ax,*bx,*cx,u)
                SHFT(*fa,*fb,*fc,fu)
        }
}
#undef GOLD
#undef GLIMIT
#undef TINY
#undef SHFT

#define ITMAX 100
#define EPS 3.0e-8
real zbrent(real (*func)(real), real x1, real x2, real tol)
{
        int iter;
        real a=x1,b=x2,c=x2,d=0.0,e=0.0,min1,min2;
        real fa=(*func)(a),fb=(*func)(b),fc,p,q,r,s,tol1,xm;

        if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
                nrerror("Root must be bracketed in zbrent");
        fc=fb;
        for (iter=1;iter<=ITMAX;iter++) {
                if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
                        c=a;
                        fc=fa;
                        e=d=b-a;
                }
                if (fabs(fc) < fabs(fb)) {
                        a=b;
                        b=c;
                        c=a;
                        fa=fb;
                        fb=fc;
                        fc=fa;
                }
                tol1=2.0*EPS*fabs(b)+0.5*tol;
                xm=0.5*(c-b);
                if (fabs(xm) <= tol1 || fb == 0.0) return b;
                if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
                        s=fb/fa;
                        if (a == c) {
                                p=2.0*xm*s;
                                q=1.0-s;
                        } else {
                                q=fa/fc;
                                r=fb/fc;
                                p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
                                q=(q-1.0)*(r-1.0)*(s-1.0);
                        }
                        if (p > 0.0) q = -q;
                        p=fabs(p);
                        min1=3.0*xm*q-fabs(tol1*q);
                        min2=fabs(e*q);
                        if (2.0*p < (min1 < min2 ? min1 : min2)) {
                                e=d;
                                d=p/q;
                        } else {
                                d=xm;
                                e=d;
                        }
                } else {
                        d=xm;
                        e=d;
                }
                a=b;
                fa=fb;
                if (fabs(d) > tol1)
                        b += d;
                else
                        b += SIGN(tol1,xm);
                fb=(*func)(b);
        }
        nrerror("Maximum number of iterations exceeded in zbrent");
        return 0.0;
}
#undef ITMAX
#undef EPS

void fit(real x[], real y[], int ndata, real sig[], int mwt, real *a,
        real *b, real *siga, real *sigb, real *chi2, real *q)
{
   //real gammq(real a, real x);
        int i;
        real wt,t,sxoss,sx=0.0,sy=0.0,st2=0.0,ss,sigdat;

        *b=0.0;
        if (mwt) {
                ss=0.0;
                for (i=1;i<=ndata;i++) {
                        wt=1.0/DSQR(sig[i]);
                        ss += wt;
                        sx += x[i]*wt;
                        sy += y[i]*wt;
                }
        } else {
                for (i=1;i<=ndata;i++) {
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
        *b /= st2;
        *a=(sy-sx*(*b))/ss;
        *siga=sqrt((1.0+sx*sx/(ss*st2))/ss);
        *sigb=sqrt(1.0/st2);
        *chi2=0.0;
        *q=1.0;
        if (mwt == 0) {
                for (i=1;i<=ndata;i++)
                        *chi2 += DSQR(y[i]-(*a)-(*b)*x[i]);
                sigdat=sqrt((*chi2)/(ndata-2));
                *siga *= sigdat;
                *sigb *= sigdat;
        } else {
                for (i=1;i<=ndata;i++)
                        *chi2 += DSQR((y[i]-(*a)-(*b)*x[i])/sig[i]);
                if (ndata>2) *q=gammq(0.5*(ndata-2),0.5*(*chi2));
        }
}

int nn;
real *xx,*yy,*sx,*sy,*ww,aa,offs;

#define BIG 1.0e30
real chixy(real bang)
{
        int j;
        real ans,avex=0.0,avey=0.0,sumw=0.0,b;

        b=tan(bang);
        for (j=1;j<=nn;j++) {
                ww[j]  = DSQR(b*sx[j]);
		ww[j] += DSQR(sy[j]);
                sumw += (ww[j] = (ww[j] < 1.0/BIG ? BIG : 1.0/ww[j]));
                avex += ww[j]*xx[j];
                avey += ww[j]*yy[j];
        }
        avex /= sumw;
        avey /= sumw;
        aa=avey-b*avex;
        for (ans = -offs,j=1;j<=nn;j++)
                ans += ww[j]*DSQR(yy[j]-aa-b*xx[j]);
        return ans;
}
#undef BIG

#define POTN 1.571000
#define BIG 1.0e30
#define PI 3.14159265
#define ACC 1.0e-3
void fitexy(real x[], real y[], int ndat, real sigx[], real sigy[],
	real *a, real *b, real *siga, real *sigb, real *chi2, real *q)
{
	int j;
	real swap,amx,amn,varx,vary,ang[7],ch[7],scale,bmn,bmx,d1,d2,r2,
		dum1,dum2,dum3,dum4,dum5;

	xx=vector(1,ndat);
	yy=vector(1,ndat);
	sx=vector(1,ndat);
	sy=vector(1,ndat);
	ww=vector(1,ndat);
	avevar(x,ndat,&dum1,&varx);
	avevar(y,ndat,&dum1,&vary);
	scale=sqrt(varx/vary);
	nn=ndat;
	for (j=1;j<=ndat;j++) {
		xx[j]=x[j];
		yy[j]=y[j]*scale;
		sx[j]=sigx[j];
		sy[j]=sigy[j]*scale;
		ww[j]=sqrt(Dsqr(sx[j])<+Dsqr(sy[j]));
	}
	fit(xx,yy,nn,ww,1,&dum1,b,&dum2,&dum3,&dum4,&dum5);
	offs=ang[1]=0.0;
	ang[2]=atan(*b);
	ang[4]=0.0;
	ang[5]=ang[2];
	ang[6]=POTN;
	for (j=4;j<=6;j++) ch[j]=chixy(ang[j]);
	mnbrak(&ang[1],&ang[2],&ang[3],&ch[1],&ch[2],&ch[3],chixy);
	*chi2=brent(ang[1],ang[2],ang[3],chixy,ACC,b);
	*chi2=chixy(*b);
	*a=aa;
	*q=gammq(0.5*(nn-2),*chi2*0.5);
	for (r2=0.0,j=1;j<=nn;j++) r2 += ww[j];
	r2=1.0/r2;
	bmx=BIG;
	bmn=BIG;
	offs=(*chi2)+1.0;
	for (j=1;j<=6;j++) {
		if (ch[j] > offs) {
			d1=fabs(ang[j]-(*b));
			while (d1 >= PI) d1 -= PI;
			d2=PI-d1;
			if (ang[j] < *b) {
				swap=d1;
				d1=d2;
				d2=swap;
			}
			if (d1 < bmx) bmx=d1;
			if (d2 < bmn) bmn=d2;
		}
	}
	if (bmx < BIG) {
		bmx=zbrent(chixy,*b,*b+bmx,ACC)-(*b);
		amx=aa-(*a);
		bmn=zbrent(chixy,*b,*b-bmn,ACC)-(*b);
		amn=aa-(*a);
		*sigb=sqrt(0.5*(bmx*bmx+bmn*bmn))/(scale*Dsqr(cos(*b)));
		*siga=sqrt(0.5*(amx*amx+amn*amn)+r2)/scale;
	} else (*sigb)=(*siga)=BIG;
	*a /= scale;
	*b=tan(*b)/scale;
	free_vector(ww,1,ndat);
	free_vector(sy,1,ndat);
	free_vector(sx,1,ndat);
	free_vector(yy,1,ndat);
	free_vector(xx,1,ndat);
}
#undef POTN
#undef BIG
#undef PI
#undef ACC

/* ============================================================ *
 * Legendre Polynomials. See gsl_sf_legendre_Pl_e.		*
 * ============================================================ */

real Pleg(int l, real x)
{
   assert(l>=0 && (-1.0<=x && x<=1.0));

   if (l==0) return 1.0;
   else if (l==1) return x;
   else if (l==2) return 0.5*(3.0*x*x - 1.0);
   else if (l==3) return 0.5*x*(5.0*x*x - 3.0);

   if (x==1.0) return 1.0;
   else if (x==-1.0) {
      if (l%2==0) return 1.0;
      else return 0.0;
   } else {

      /* recurrence relation */
      int ell;
      real p_ellm2 = 1.0;       /* P_0(x) */
      real p_ellm1 = x;         /* P_1(x) */
      real p_ell   = p_ellm1;

      for (ell=2; ell<=l; ell++) {
	 p_ell = ((2*ell-1)*x*p_ellm1 - (ell-1)*p_ellm2)/ell;
	 p_ellm2 = p_ellm1;
	 p_ellm1 = p_ell;
      }

      return p_ell;

   }
      
}


/* ============================================================ *
 * Gauss, normalized correctly for 1d				*
 * ============================================================ */

double gaussfunction(double xsqr, double sigmasqr)
{
   return 1.0/sqrt(twopi*sigmasqr)*exp(-xsqr/(2.0*sigmasqr));
}

void donothing(double x)
{
   return;
}
