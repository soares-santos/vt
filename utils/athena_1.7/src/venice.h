#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/times.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>

#define FAILURE 0
#define SUCCESS 1

#define PI 3.14159265358979323846
#define TWOPI 6.283185307179586476925287

//#define EPS 1.0e-10
#define INF 1.0e30
#define ODD 0
#define EVEN 1
#define LEAF 0
#define NODE 1
#define RADEC 0
#define CART 1

#define MAX(x,y) ((x) > (y)) ? (x) : (y)
#define MIN(x,y) ((x) < (y)) ? (x) : (y)
#define ABS(a) ((a) < 0 ? -(a) : (a))
#define PARITY(a) (a)%2 ? ODD : EVEN

#define NFIELD 100
#define NCHAR 20

#define getDoubleValue(array,col)  atof(array+NCHAR*(col-1))
#define getIntValue(array,col)  atoi(array+NCHAR*(col-1))
#define getCharValue(array,col) array+NCHAR*(col-1)
#define getLine(array,i) array+NFIELD*NCHAR*i

#define MAX(x,y) ((x) > (y)) ? (x) : (y)
#define MIN(x,y) ((x) < (y)) ? (x) : (y)
#define ABS(a) ((a) < 0 ? -(a) : (a))
#define SWAP(a,b) {swap = (a); (a) = (b); (b) = swap;}


char MYNAME[100];
size_t IDERR;
double EPS=99.00;

/*----------------------------------------------------------------*
 *New types                                                       *
 *----------------------------------------------------------------*/


typedef struct Config
{
  char fileRegInName[1000];
  char fileCatInName[1000];
  char fileOutName[1000];
  int nx,ny,format;
  int xcol,ycol;
  int coordType;
  size_t npart;
  double min[2];
  double max[2];
  int minDefinied[2];
  int maxDefinied[2];
} Config;


typedef struct Complex
{
    double re;
    double im;
} Complex;

typedef struct Polygon
{
  int N;
  double x[100];
  double y[100];
  double xmin;
  double xmax;
  double ymin;
  double ymax;
} Polygon;

/*----------------------------------------------------------------*
 *Global variables                                                *
 *----------------------------------------------------------------*/
FILE *FILE_REG_IN,*FILE_CAT_IN, *FILE_OUT;
int NX,NY,XCOL,YCOL,NPART,FORMAT;
double *LIMITS;
/*----------------------------------------------------------------*
 *Main routines                                                   *
 *----------------------------------------------------------------*/

gsl_histogram2d *mask2d(const Config *para);
void flagCat(const Config *para);
void randomCat(const Config *para);
int readParameters(int argc, char **argv, Config *para);

/*----------------------------------------------------------------*
 *Utils - geometric                                               *
 *----------------------------------------------------------------*/

int insidePolygon(Polygon *p,int Npoly, double x0,double y0,double x,double y, int *poly_id);
Polygon *readPolygonFile(FILE *file_in,int *Npoly);

/*----------------------------------------------------------------*
 *Utils - numeric                                                 *
 *----------------------------------------------------------------*/

double determineMachineEpsilon();
size_t determineSize_tError();
gsl_rng *randomInitialize(int time_dependent);
FILE *fopenAndCheck(const char *filename,char *mode);
int getStrings(char *line, char *strings, char *delimit, size_t *N);
void printCount(const size_t *count, const size_t *total,  const size_t step);
