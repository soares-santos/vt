/* ============================================================ *
 * 2pcf.h							*
 * Martin Kilbinger 2006, 2007, 2008				*
 * ============================================================ */

#ifndef __2PCF_H
#define __2PCF_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/times.h>
#include <time.h>
#include <string.h>
#include <math.h>

#include "mathstuff.h"
#include "nrcomplex.h"
#include "config.h"
#include "gal_cat.h"
#include "io.h"
#include "errorlist.h"
#include "stdnames.h"


/* Error codes */
#define tpcf_outofrange -10
#define tpcf_ncat       -11
#define tpcf_bintype    -12
#define tpcf_npositive  -13
#define tpcf_ndim       -14



typedef enum {LIN, LOG} bin_type;

/* Bit-coded correlation types.
   gg: shear-shear
   gn: shear-position
   nn: position-position (angular correlation)
   wn: weight-position (weighted angular correlation)

   Note: wn should not be used in combination with another type!
*/
typedef enum {gg=1, gn=2, nn=4, wn=8} wcorr_t;

typedef enum {error_none, bootstrap, jackknife} error_t;
#define serror_t(i) (		 \
   i==error_none ? "none"      : \
   i==bootstrap  ? "bootstrap" : \
   i==jackknife  ? "jackknife" : \
   "")
#define Nerror_t 3

typedef enum {wcorr_subtype_undef, nn_2d, nn_3d, nn_rp_pi, wn_all, wn_pairwise} wcorr_subtype_t;
#define swcorr_subtype_t(i) ( \
  i==wcorr_subtype_undef ? "undef" : \
  i==nn_2d               ? "nn_2d" : \
  i==nn_3d               ? "nn_3d" : \
  i==nn_rp_pi            ? "nn_rp_pi" : \
  i==wn_all              ? "wn_all" : \
  i==wn_pairwise         ? "wn_pairwise" : \
  "")
#define Nwcorr_subtype 6


typedef struct {

  /* min, max are the minimal and maximal separations that galaxy pairs can have. * 
   * So they are the smalles and largest bin corners.                             */

  double min[KTREE_MAX], max[KTREE_MAX];         /* min, max for both LIN and LOG (NEW!)    */
  double minsqr[KTREE_MAX], maxsqr[KTREE_MAX];   /* min^2, max^2 for both LIN and LOG!      */
  double logmin[KTREE_MAX], logmax[KTREE_MAX];   /* log(min), log(max) for both LIN and LOG */
  double diff[KTREE_MAX], difflog[KTREE_MAX];    /* max-min, log(max)-log(min)              */
  uint N[KTREE_MAX];				 /* Number of bins                          */
  uint ndim;					 /* Dimensionality of bins. Not to confuse  *
						  * with dimensionality of data. Typical    *
						  * value is 1.				    */
  bin_type linlog;                               /* = LIN or LOG			    */

} bin_data;

typedef struct {

  unsigned long long int *npair;
  double **npair_resample;
  double **p_resample, **m_resample, **ww_resample, **wn_resample, **gt_resample;

  double *p, *m, *x, *g11, *g22, *g12, *g21;
  double *ww, *wn, *wwww, *wnwn, *wtheta;
  int *Nnode;

  double *gt, *gx;
} xi_data;


void read_all_gal(const char *name1, const char *name2, uint *ngal1, uint *ngal2, int ncat, int radec, format_t format,
		  coord_t coord_input, galaxy_data **gal1, galaxy_data **gal2, int ndim, const char **col_names,
		  int Ncol, int quiet, error **err);
bin_data *init_bin(const double *thmin, const double *thmax, const uint *nth, int ndimbin,
		   bin_type I_BINTYPE, error **err);
uint bin_N(const bin_data *bin);
int distance2bin_index(double distance, double min, double diff, int N, int round);
int get_bin_index(const int *b, const uint *N, int ndim);
double bin_index_to_scale(const bin_data *bin, const double *wtheta, int b);

double bin_index_to_scale_low(const bin_data *bin, const double *wtheta, int b);
double bin_index_to_scale_high(const bin_data *bin, const double *wtheta, int b);

double calculate_sqrt_D(double sig_eps4, double wwww, double ww, int ncat);
double calculate_sqrt_Dcor(double sqrt_D, unsigned long long int npair, int Nnode);

xi_data *init_xi(const bin_data *bin, wcorr_t wcorr, int Nboots, int do_wtheta, error **err);
void weigh(const bin_data *bin, xi_data *xi, wcorr_t wcorr, int nresample, error **err);
double sigma_epsilon_sqr(const galaxy_data *gal, uint ngal);

double get_resample_error(const double **xi, int b, int nresample, error_t ERROR, double *mean_resample, error **err);
void get_resample_cov_xi(double *cov_ppmmpm, const xi_data *xi, int b, int c, int nresample, error_t ERROR, error **err);
double get_resample_cov_gl(const xi_data *xi, int b, int c, int nresample, error_t ERROR, error **err);
double get_boots_error(const xi_data *xi, int b, int nboots);

void fill_header(char str_header[], const char *header[], int n, int signif, int nresample, error_t ERROR, coord_t coord, int n_coord);
void out_xi(const xi_data *xi, const char *xi_name, const char *xiref_name, const char *xi_name2 , const bin_data *bin,
	    double sig_eps4, int ncat, coord_t coord, error **err);
void out_xi_resample(const xi_data *xi, const char *xi_name, const char *out_ALL_xip_resample,const char *out_ALL_xim_resample,
		     const bin_data *bin, coord_t coord_output, int nresample, error_t ERROR, error **err);
void out_cov_xi_resample(const xi_data *xi, const char *cov_name, const bin_data *bin, coord_t coord_output, int nresample,
			 error_t ERROR, error **err);
void out_cov_gl_resample(const xi_data *xi, const char *cov_name, const bin_data *bin, coord_t coord_output, int nresample,
          error_t ERROR, error **err);
void out_w(const xi_data *xi, const char *w_name, const bin_data *bin, int ngal1, int ngal2,
           coord_t coord_output, int nresample, error_t ERROR, error **err);
void out_wp_rp_pi(const xi_data *xi, const char *w_name, const bin_data *bin, int nboots, coord_t coord_output,
		  uint ngal1, uint ngal2, error **err);
void out_ww(const xi_data *xi, const char *ww_name, const bin_data *bin, int nboots, coord_t coord, double wbar, error **err);
void out_gl(const xi_data *xi, const char *gl_name, const bin_data *bin, double sig_eps2, coord_t coord, error **err);
void out_gl_resample(const xi_data *xi, const char *gl_name, const bin_data *bin,const char *out_ALL_gt_resample,
		     coord_t coord_output, int nresample, error_t ERROR, error **err);

void free_xi_data(xi_data *xi, uint N);

int check_input_format(wcorr_t wcorr, format_t format, const char **col_names, int Ncol, error **err);
int check_input_format_format(wcorr_t wcorr, format_t format, error **err);


#ifdef _WITH_FITS
void out_xi_fits(const xi_data *xi, const char *xi_name, const bin_data *bin,
                 double sig_eps4, int ncat, coord_t coord_output, int nresample, error_t ERROR, int quiet, error **err);
void out_gl_fits(const xi_data *xi, const char *gl_name, const bin_data *bin, double sig_eps2,
		 coord_t coord_output, int nresample, error_t ERROR, int quiet, error **err);
void out_w_fits(const xi_data *xi, const char *w_name, const bin_data *bin,
                 int ngal1, int ngal2, coord_t coord_output, int nresample,
                 error_t ERROR, int quiet, error **err);
int check_input_format_fits(wcorr_t wcorr, const char **col_names, int Ncol, error **err);
void write_header_keys(fitsfile *fptr, char *object, const char *object_comment,
		       coord_t coord, int nresample, error_t ERROR, int quiet, error **err);
void write_columns(fitsfile *fptr, const int *data_type, void *ptr_data[], const char *ttype[],
                    int Nbin, int Ncol, int quiet, error **err);
#endif

#endif
