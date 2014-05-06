#ifndef _VORONOI_H_ 
#define _VORONOI_H_

#ifdef SINGLE
#define REAL float
#else 
#define REAL double
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "fitsio.h"
#include "triangle.h"
#include "nr.h"

#define ARRAY1(TYPE,A,N) {\
  A=(TYPE *) malloc(N * sizeof(TYPE));\
  if((A)==NULL) puts("Memory problem"), exit(1);\
  }
#define ARRAY2(TYPE,ARR,NR,NC) {\
  register int i;\
  ARR = (TYPE **)malloc(NR * sizeof(TYPE *));\
  ARR[0] = (TYPE *)malloc(NR * NC * sizeof(TYPE));\
  if(ARR[0]==NULL) puts("Memory problem"),exit(1);\
  for(i=1; i<NR; i++) ARR[i] = ARR[0] + i * NC;\
  }
#define free2(ARR) {\
  free(ARR[0]);\
  free(ARR);\
  }
#define MAXMIN(a,max,min) {if((a)>(max)) (max)=(a); if((a)<(min)) (min)=(a);}
#define MIN(a,min) {if((a)<(min)) (min)=(a);}
#define nextline(fp) {while(getc(fp) != '\n');}
#define skipcomments(fp) {while(c=getc(fp),c=='#')nextline(fp); ungetc(c,fp);}
#define ADDLIST(nl,i) {\
  list=(int *)realloc(list,sizeof(int)*(nl));\
  list[(nl)-1]=(i);}
#define CMPX(A,B) \
        v = (*(coord**)A)->x - (*(coord**)B)->x;\
	if (v>0) return 1;\
	if (v<0) return -1;
#define CMPY(A,B) \
	v = (*(coord**)A)->y - (*(coord**)B)->y;\
	if (v>0) return 1;\
        if (v<0) return -1;
#define TF(x) (gammq(alpha,beta/x)) /* density distrib */

typedef struct triangulateio T_io;
typedef struct footprint {
  REAL xc, yc,             /* Coordinates of the baricenter */
    radius,                /* max dist from baricenter */
    area,                  /* total area of the footprint */
    minflux,               /* Minimum cell density */
    sig,                   /* Significativity */
    pa,                    /* position angle */
    semi_maj,              /* semi_mag axis, ellipse centered at baricenter */
    semi_min,              /* semi_min axis */
    xgc,ygc,zgc,           /* Coords of the central glx (the smallest cell) */
    contrast,              /* density contrast */
    z_rms,                 /* rms of redshift */
    z_mean,z_median;       /* mean and median of redshits */
  int nvor,                /* Nr of pts included by Voronoi */
    nbg,                   /* Nr of background pts in the Voronoi area */
    n,                     /* Nr of pts included by Voronoi - background */ 
    edge;                  /* cluster is beyond the frame region */
  long long id;            /* cluster id */
} Footprint;
typedef struct connect {
  int n, *list;
  REAL sig;                /* Significativity */
  struct connect *link;
} Connect;
typedef struct Coord {
  REAL x,y;
} coord;

/* general declarations */

coord *ch_points, **ch_P;
T_io in, mid, vor;
char *program_name, *file_name, *file_holes;
REAL xmean, ymean, zmean, xmax, xmin, ymax, ymin,  /* Points limits */
  area_mean, area_mean_ok,                         /* Polygon mean area */
  *mag, *mag1, magmin, magmax, medmag,             /* Magnitudes */
  *redsh, *redsh1, redshmin, redshmax, medredsh,   /* Redshifts */
  *area, tot_area_covered,                         /* Polygon surfaces */
  Nbg,                                             /* Background counts */
  background_density,                              /* Background density */
  conf_lev,                                        /* Confidence level */
  reject_lev,                                      /* Rejection Conf. level */
  threshold,                                       /* detection threshold */
  *HolesCoord,                                     /* Holes coordinates */
  w_amp,w_pow,
  alpha,beta,
  c0,c1,c2,c3,
  frame[4];
int ndata, ndata_ok, nholes, *edge, *prox, Verbose, ndataused, 
  Holes, XYcoord, Ebeling, Poisson, Fits,NotFind,RawFind,
  *mask1,                                          /* Label for each galaxy */
  *central,*edge1,                                 /* flag for each galaxy  */
  nclusters;
long long *number, FieldID,*host_id_list_1;
short  *voronoi_inside;                           /* Label for voronoi edges */
Footprint *footprint;
char *ttype[4],*ttypeH[4];

/* prototypes */

void set_defaults(void);
void read_data(void);
REAL median(REAL *, int);
int real_compare(const void *i, const void *j);
void cfitsio_error(int);
void read_holes(void);
void set_TF_params(void);
REAL poly_area(int);
int intersect(REAL, REAL, REAL, REAL, REAL, REAL, REAL, REAL);
void WriteAreasFile(void);
REAL ComputeThreshold( REAL *);
REAL FitParabol(int, REAL *, REAL *, REAL *, REAL *, REAL *); 
void WriteDistributionsFile(REAL *);
void WriteCandidatesFile(int *, int *);
Connect *percola(int *, int *);
Connect *AddItem (Connect *, int, REAL, int *);
int CountItem (Connect *);
void WriteVoronoiCatalog();
void WriteDelaunayCatalog();
void ComputeRedshifts(Connect *);
REAL rms(REAL *,int ,REAL);
void WriteSimpleMembersList(void);
void WriteGalaxyCatalogAscii(void);
void WriteClusterCatalogAscii();
void WriteGalaxyCatalog(void);
void WriteClusterCatalog(void);
void SetFootprints(Connect *);
int ch2d(coord **, int);
int inpoly(REAL *, int, REAL, REAL);
int cmpl(const void *a, const void *b);
int cmph(const void *a, const void *b);
int make_chain(coord** V, int n, int (*cmp)(const void*, const void*)); 
int ccw(coord **ch_P, int i, int j, int k);
void print_init(void);
void adjust_frame(void);

#endif /* _VORONOI_H_ */
