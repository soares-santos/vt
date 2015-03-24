/* ============================================================ *
 * main.h							*
 * Martin Kilbinger, Christopher Bonnett 2006-2010		*
 * ============================================================ */

#ifndef __TREE_H
#define __TREE_H


#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <getopt.h>

#include "mathstuff.h"
#include "nrcomplex.h"
#include "config.h"
#include "errorlist.h"
#include "2pcf.h"
#include "kdtree.h"
#include "gal_cat.h"
#include "io.h"
#include "stdnames.h"

extern char *optarg;


/* Error codes */
#define tree_outofrange -1
#define tree_bintype    -2
#define tree_radec      -3
#define tree_ngal       -4
#define tree_id_nodes   -5
#define tree_nboots     -6
#define tree_wcorr      -7
#define tree_jack_NJ    -8
#define tree_error_t    -9

#define VERSION 1.7


#define print_err_log(LOG, my_quiet, str, ...) char msg_pel[TXT_SZ]; sprintf(msg_pel, str, __VA_ARGS__); print_err_log_func(LOG, my_quiet, msg_pel)


typedef struct {

  char GALCAT1[2048], GALCAT2[2048];
  char **COL_NAMES;
  double THMIN, THMAX, thmin[NDIMBIN_MAX], thmax[NDIMBIN_MAX], OATH;
  int RADEC, NCAT, NRESAMPLE[2], NCOL;
  uint NTH, nth[NDIMBIN_MAX];
  format_t FORMAT;
  coord_t COORD_INPUT, COORD_OUTPUT;
  char BINTYPE[64], SFORMAT[64], SCOORD_INPUT[64], SCOORD_OUTPUT[64], SERROR[64];
  char SWCORR_SUBTYPE[64], OUT_MODE[64];
  bin_type I_BINTYPE;
  wcorr_t WCORR;
  wcorr_subtype_t WCORR_SUBTYPE;
  error_t ERROR;
  double WN_ZGAP_MIN, WN_ZGAP_MAX;

  int do_wtheta;

} config_info;


void read_config_file(config_info *config, char cname[], int *ndim, int *ndimbin, error **err);
void out_config(FILE *LOG, const char *cname, config_info config, int ndimbin);
void print_err_log_func(FILE *LOG, int my_quiet, const char *str);
void print_err_log_VA(FILE *LOG, int my_quiet, const char *str, ...);

void init_tree(node **root1, node **root2, galaxy_data *gal1, galaxy_data *gal2,
	       uint ngal1, uint ngal2, int ndim, int twins_merge, const config_info *config, int quiet, FILE *LOG);

void correl_two_nodes(const node *nd1, const node *nd2, double distance, 
		      const bin_data *bin, const config_info *config, xi_data *xi);
void twopcf(const galaxy_data *gal1, const galaxy_data *gal2, uint ngal1, uint ngal2,
	    int ndim, node *nd1, node *nd2, const bin_data *bin,
	    const config_info *config,  xi_data *xi);

/* Bootstrap and Jackknife */
void fill_bootstrap_sample(int NRESAMPLE, wcorr_t WCORR, galaxy_data *gal, uint ngal, int quiet, error **err);
int update_nresample(int *NRESAMPLE, int i, int njack, int quiet);
void init_jackboot(config_info *config, galaxy_data *gal1, galaxy_data *gal2, uint ngal1, uint ngal2,
		   int quiet, error **err);
void allocate_jackboot_sample(int NRESAMPLE, galaxy_data *gal, uint ngal, error **err);
void Nxy_jackknife(uint N, double r, uint *Nx, uint *Ny);
int fill_jackknife_sample(int NJACK, galaxy_data *gal, uint ngal, const double *min, const double *max,
			  error **err);
int fill_jackknife_sample_idx(int NJACK, galaxy_data *gal, uint ngal, error **err);

void set_char_if_null(char **str, const char *name_to_set, error **err);
void some_stats(int i, int ndim, const node *root, int radec, FILE *LOG);
void write_all_correlations(const xi_data *xi, const char *names[], const bin_data *bin, const config_info *config,
			    double sig_eps2, double sig_eps4, double wbar, uint ngal1, uint ngal2, error **err);
#ifdef _WITH_FITS
void write_all_correlations_fits(const xi_data *xi, const char *names[], const bin_data *bin, const config_info *config,
			    double sig_eps2, double sig_eps4, double wbar, uint ngal1, uint ngal2, int quiet, error **err);
#endif

void usage(int ex);


#endif
