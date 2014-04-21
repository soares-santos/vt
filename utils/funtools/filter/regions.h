#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#ifdef __STDC__
#include <stdlib.h>
#include <stdarg.h>
#else
#include <varargs.h>
#endif

#define MASKINC 	10000
#ifndef M_PI
#define M_PI		3.14159265358979323846
#endif
#define SMALL_NUMBER	1.0E-24
#define LARGE_NUMBER	65535
#define PSTOP		-142857.142857

#ifndef SZ_LINE
#define SZ_LINE 	4096
#endif
#ifndef min
#define min(x,y)	(((x)<(y))?(x):(y))
#endif
#ifndef max
#define max(x,y)	(((x)>(y))?(x):(y))
#endif
#ifndef abs
#define abs(x)		((x)<0?(-x):(x))
#endif
#ifndef feq
#define feq(x,y)	(fabs((double)x-(double)y)<=(double)1.0E-15)
#endif
#ifndef NULL
#define NULL 		(void *)0
#endif

#ifndef TOK_EREG
#define TOK_EREG	1
#endif
#ifndef TOK_NREG
#define TOK_NREG	2
#endif
#ifndef TOK_IREG
#define TOK_IREG	4
#endif

#define PIXCEN(a)	(double)(a)
#define PIXNUM(a)	(int)((a)+0.5) 
#define PIXSTART(a)	((int)(a)+1)
#define PIXSTOP(a)	(((int)(a))==(a)?((int)(a)-1):((int)(a)))
/* to assure that geometrically adjoining regions touch but don't overlap */
/* when edge is exactly on a pixel center it goes to right or upper region. */
/* used for non-radially symetric regions instead of PIXSTART, PIXSTOP */
#define PIXINCL(a)	(int)((a)+1.0) 

/* this is the filter string for field only */
#define EVFIELDONLY "(evfield(g,1,1,1,(double)x,(double)y))"

#define XSNO 3

/* NB: these MUST match the definition in filter.h */
#ifndef __filter_h
typedef struct filtmaskrec {
  int region;
  int y;
  int xstart, xstop;
} *FilterMask, FilterMaskRec;

/* parameter structure for a scan entry */
typedef struct scanrec{
  struct scanrec *next;
  int x;
} *Scan, ScanRec;

typedef struct shaperec {
  int init;
  double ystart, ystop;
  Scan *scanlist;
  /* varargs */
  int nv;
  double *xv;
  /* circle, annulus */
  double r1sq, r2sq;
  /* ellipse */
  double angl, sinangl, cosangl;
  double cossq, sinsq;
  double xradsq, yradsq;
  double a;
  /* polygon-style shapes */
  int npt;
  double *pts;
  /* line */
  int xonly;
  double x1, x2, y1;
  double invslope;
} *Shape, ShapeRec;

/* these are global for use with special region routines */
typedef struct gfiltrec {
  int nshapes;			/* number of shapes */
  int maxshapes;		/* number of shape records we allocate */
  Shape shapes;			/* array holding range limits for one shape */
  int rid;			/* first valid region for current pixel */
  int usebinsiz;		/* whether bindizx,binsizy are used */
  char *evsect;			/* value of event section */
  double tlminx, tlminy;	/* tlmin for event section */
  double binsizx, binsizy;	/* bin sizes for event section */
  double tloff;		        /* offset for quick p2i conversion */
  int xmin, xmax, ymin, ymax;	/* section limits in original image coords */
  int block;			/* block factor */
  int x0, x1, y0, y1;		/* section limits in section coords */
  int *ybuf;			/* valid y row flags */
  int *x0s;			/* valid x start values */
  int *x1s;			/* valid x stop values */
  int nmask;			/* number of image mask record */
  int maskdim;			/* size of mask image */
  FilterMask masks;		/* mask records */
} *GFilt, GFiltRec;
#endif

/* declare image init routines */
void imannulusi(GFilt g, int rno, int sno, int flag, int type, 
		double x, double y,
		double xcen, double ycen, double iradius, double oradius);
void imboxi(GFilt g, int rno, int sno, int flag, int type,
	    double x, double y,
	    double xcen, double ycen, double xwidth, double yheight,
	    double angle);
void imcirclei(GFilt g, int rno, int sno, int flag, int type,
	       double x, double y,
	       double xcen, double ycen, double radius);
void imellipsei(GFilt g, int rno, int sno, int flag, int type,
		double x, double y,
		double xcen, double ycen, double xrad, double yrad,
		double angle);
void imfieldi(GFilt g, int rno, int sno, int flag, int type,
	      double x, double y);
void imlinei(GFilt g, int rno, int sno, int flag, int type,
	     double x, double y,
	     double x0, double y0, double x1, double y1);
void impiei(GFilt g, int rno, int sno, int flag,  int type,
	    double x, double y,
	    double xcen, double ycen, double angle1, double angle2);
void imqtpiei(GFilt g, int rno, int sno, int flag,  int type,
	      double x, double y,
	      double xcen, double ycen, double angle1, double angle2);
void impointi(GFilt g, int rno, int sno, int flag, int type,
	      double x, double y,
	      double xcen, double ycen);
void impandai(GFilt g, int rno, int sno, int flag, int type,
	      double x, double y,
	      double xcen, double ycen,
	      double anglo, double anghi, double angn,
	      double radlo, double radhi, double radn);
void imnannulusi(GFilt g, int rno, int sno, int flag, int type,
		 double x, double y,
		 double xcen, double ycen,
		 double lo, double hi, int n);
void imnboxi(GFilt g, int rno, int sno, int flag, int type,
	     double x, double y,
	     double xcen, double ycen,
	     double lox, double loy, double hix, double hiy, int n,
	     double angle);
void imnellipsei(GFilt g, int rno, int sno, int flag, int type,
		 double x, double y,
		 double xcen, double ycen,
		 double lox, double loy, double hix, double hiy, int n,
		 double angle);
void imnpiei(GFilt g, int rno, int sno, int flag, int type,
	     double x, double y,
	     double xcen, double ycen,
	     double lo, double hi, int n);

#ifdef __STDC__
void impolygoni(GFilt g, int rno, int sno, int flag, int type,
		double x, double y, ...);
void imvannulusi(GFilt g, int rno, int sno, int flag, int type,
		 double x, double y, double xcen, double ycen, ...);
void imvboxi(GFilt g, int rno, int sno, int flag, int type,
	     double x, double y, double xcen, double ycen, ...);
void imvellipsei(GFilt g, int rno, int sno, int flag, int type,
		 double x, double y, double xcen, double ycen, ...);
void imvpiei(GFilt g, int rno, int sno, int flag, int type,
	     double x, double y, double xcen, double ycen, ...);
void imvpointi(GFilt g, int rno, int sno, int flag, int type, 
	       double x, double y, ...);
#endif


/* declare image region routines */
int imannulus(GFilt g, int rno, int sno, int flag, int type,
	      double x, double y,
	      double xcen, double ycen, double iradius, double oradius);
int imbox(GFilt g, int rno, int sno, int flag, int type,
	  double x, double y,
	  double xcen, double ycen, double xwidth, double yheight,
	  double angle);
int imcircle(GFilt g, int rno, int sno, int flag, int type,
	     double x, double y,
	     double xcen, double ycen, double radius);
int imellipse(GFilt g, int rno, int sno, int flag, int type,
	      double x, double y,
	      double xcen, double ycen, double xrad, double yrad,
	      double angle);
int imfield(GFilt g, int rno, int sno, int flag, int type,
	    double x, double y);
int imline(GFilt g, int rno, int sno, int flag, int type,
	   double x, double y,
	   double x1, double y1, double x2, double y2);
int impie(GFilt g, int rno, int sno, int flag, int type,
	  double x, double y,
	  double xcen, double ycen, double angle1, double angle2);
int imqtpie(GFilt g, int rno, int sno, int flag, int type,
	    double x, double y,
	    double xcen, double ycen, double angle1, double angle2);
int impoint(GFilt g, int rno, int sno, int flag, int type,
	    double x, double y,
	    double xcen, double ycen);
int impanda(GFilt g, int rno, int sno, int flag, int type,
	     double x, double y,
	     double xcen, double ycen,
	     double anglo, double anghi, double angn,
	     double radlo, double radhi, double radn);
int imnannulus(GFilt g, int rno, int sno, int flag, int type,
	       double x, double y,
	       double xcen, double ycen,
	       double lo, double hi, int n);
int imnbox(GFilt g, int rno, int sno, int flag, int type,
	   double x, double y,
	   double xcen, double ycen,
	   double lox, double loy, double hix, double hiy, int n,
	   double angle);
int imnellipse(GFilt g, int rno, int sno, int flag, int type,
	       double x, double y,
	       double xcen, double ycen,
	       double lox, double loy, double hix, double hiy, int n,
	       double angle);
int imnpie(GFilt g, int rno, int sno, int flag, int type,
	   double x, double y,
	   double xcen, double ycen,
	   double lo, double hi, int n);
#ifdef __STDC__
int impolygon(GFilt g, int rno, int sno, int flag, int type,
	      double x, double y, ...);
int imvannulus(GFilt g, int rno, int sno, int flag, int type,
	       double x, double y, double xcen, double ycen, ...);
int imvbox(GFilt g, int rno, int sno, int flag, int type,
	   double x, double y, double xcen, double ycen, ...);
int imvellipse(GFilt g, int rno, int sno, int flag, int type,
	       double x, double y, double xcen, double ycen, ...);
int imvpie(GFilt g, int rno, int sno, int flag, int type,
	   double x, double y, double xcen, double ycen, ...);
int imvpoint(GFilt g, int rno, int sno, int flag, int type,
	     double x, double y, ...);
#endif

/* declare event region routines */
int evannulus(GFilt g, int rno, int sno, int flag, int type,
	      double x, double y,
	      double xcen, double ycen, double iradius, double oradius);
int evbox(GFilt g, int rno, int sno, int flag, int type,
	  double x, double y,
	  double xcen, double ycen, double xwidth, double yheight,
	  double angle);
int evcircle(GFilt g, int rno, int sno, int flag, int type,
	     double x, double y,
	     double xcen, double ycen, double radius);
int evellipse(GFilt g, int rno, int sno, int flag, int type,
	      double x, double y,
	      double xcen, double ycen, double xrad, double yrad,
	      double angle);
int evfield(GFilt g, int rno, int sno, int flag, int type,
	    double x, double y);
int evline(GFilt g, int rno, int sno, int flag, int type,
	   double x, double y,
	   double x1, double y1, double x2, double y2);
int evpie(GFilt g, int rno, int sno, int flag, int type,
	  double x, double y,
	  double xcen, double ycen, double angle1, double angle2);
int evqtpie(GFilt g, int rno, int sno, int flag, int type,
	    double x, double y,
	    double xcen, double ycen, double angle1, double angle2);
int evpoint(GFilt g, int rno, int sno, int flag, int type,
	    double x, double y,
	    double xcen, double ycen);
int evnannulus(GFilt g, int rno, int sno, int flag, int type,
	       double x, double y,
	       double xcen, double ycen,
	       double lo, double hi, int n);
int evnbox(GFilt g, int rno, int sno, int flag, int type,
	   double x, double y,
	   double xcen, double ycen,
	   double lox, double loy, double hix, double hiy, int n,
	   double angle);
int evnellipse(GFilt g, int rno, int sno, int flag, int type,
	       double x, double y,
	       double xcen, double ycen,
	       double lox, double loy, double hix, double hiy, int n,
	       double angle);
int evnpie(GFilt g, int rno, int sno, int flag, int type,
	   double x, double y,
	   double xcen, double ycen,
	   double lo, double hi, int n);
int evpanda(GFilt g, int rno, int sno, int flag, int type,
	    double x, double y,
	    double xcen, double ycen,
	    double anglo, double anghi, double angn,
	    double radlo, double radhi, double radn);
#ifdef __STDC__
int evpolygon(GFilt g, int rno, int sno, int flag, int type,
	      double x, double y, ...);
int evvannulus(GFilt g, int rno, int sno, int flag, int type,
	       double x, double y, double xcen, double ycen, ...);
int evvbox(GFilt g, int rno, int sno, int flag, int type,
	   double x, double y, double xcen, double ycen, ...);
int evvellipse(GFilt g, int rno, int sno, int flag, int type,
	       double x, double y, double xcen, double ycen, ...);
int evvpie(GFilt g, int rno, int sno, int flag, int type,
	   double x, double y, double xcen, double ycen, ...);
int evvpoint(GFilt g, int rno, int sno, int flag, int type,
	     double x, double y, ...);
#endif

