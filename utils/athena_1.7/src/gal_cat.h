/* ============================================================ *
 * gal_cat.h							*
 * Martin Kilbinger 2006-2008					*
 * ============================================================ */

#ifndef __GAL_CAT_H
#define __GAL_CAT_H

#include <sys/times.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <limits.h>

#include "stdnames.h"
#include "config.h"

#ifdef _WITH_FITS
#include <fitsio.h>

#define FITS_STATUS_ERROR(status, quiet, message, err, line, ret_val, ...) { \
	if (!quiet) fits_report_error(stderr, status); \
	testErrorRetVA(status != 0, gc_fits, message, err, line, ret_val, __VA_ARGS__); \
}

#endif


#include "mathstuff.h"
#include "errorlist.h"
#include "io.h"
#include "nrcomplex.h"
#include "kdtree.h"


#define gc_io           -20
#define gc_ngal         -21
#define gc_cat_entry    -22
#define gc_type         -23
#define gc_coord        -24
#define gc_dim          -25
#define gc_overflow     -26
#define gc_fits         -27
#define gc_format       -28


/* if modified, also change kdtree.c:swap_gal ! */
typedef struct {

  double pos[KTREE_MAX];
  double cos[KTREE_MAX], sin[KTREE_MAX];   /* Only used for RADEC=1 */
  fdcomplex ell;
  fdreal weight;
  unsigned short int idx;    /* Index, e.g. for Jackknife        */
  //double z, dzu, dzl;      /* Redshift, upper and lower errors */
  double *ind_resample;      /* Bootstrap/Jackknife index        */

} galaxy_data;

#ifndef __dsqr__
#define __dsqr__
double dsqr(double);
#endif

typedef enum {f_standard, f_hamana, f_position, f_position_jack_num, f_lensing_jack_num, f_fits} format_t;
#define sformat_t(i) ( \
  i==f_standard          ? "standard" : \
  i==f_hamana            ? "hamana"   : \
  i==f_position          ? "position" : \
  i==f_position_jack_num ? "position_jack_num" : \
  i==f_lensing_jack_num  ? "lensing_jack_num" :  \
  i==f_fits              ? "fits" : \
"")
#define Nformat_t 6

typedef enum {coord_arcsec, coord_arcmin, coord_deg, coord_rad} coord_t;
#define scoord_t(i) ( \
  i==coord_arcsec ? "arcsec" : \
  i==coord_arcmin ? "arcmin" : \
  i==coord_deg    ? "deg" : \
  i==coord_rad    ? "rad" : \
  "")
#define Ncoord_t 4

typedef enum {c_x, c_y, c_z, c_e1, c_e2, c_w, c_njk} col_names_t;
#define scol_names_t(i) ( \
  i==c_x   ? "x"   : \
  i==c_y   ? "y"   : \
  i==c_z   ? "z"   : \
  i==c_e1  ? "e1"  : \
  i==c_e2  ? "e2"  : \
  i==c_w   ? "w"   : \
  i==c_njk ? "njk" : \
  "" \
)
#define Ncol_names_t 7

#define tcol_fits_types_t(i) ( i==c_njk ? TUSHORT : TDOUBLE )

typedef enum {cd_x, cd_y, cd_z, cd_e1, cd_e2, cd_w, cd_njk} col_names_default_t;
#define scol_names_default_t(i) ( \
  i==cd_x   ? COL_NAME_DEFAULT_ra  : \
  i==cd_y   ? COL_NAME_DEFAULT_dec : \
  i==cd_z   ? COL_NAME_DEFAULT_z   : \
  i==cd_e1  ? COL_NAME_DEFAULT_e1  : \
  i==cd_e2  ? COL_NAME_DEFAULT_e2  : \
  i==cd_w   ? COL_NAME_DEFAULT_w   : \
  i==cd_njk ? COL_NAME_DEFAULT_njk : \
  "" \
)
#define Ncol_names_t 7

typedef struct {
   col_names_t key;
   char val[128];
} col_descr;

galaxy_data *read_gal(const char *name, uint *ngal, int radec, format_t format, coord_t coord_input,
		      int check_zero_weight, int ndim, const char **col_names, int Ncol,
		      int quiet, error **err);
void post_process_galaxy(galaxy_data *gal, uint ngal, coord_t coord_input, int radec, int check_zero_weight, int ndim, error **err);
void out_galaxy(const void *data, uint i, int ndim, FILE *OUT);
void out_cat_jack(const galaxy_data *gal, uint ngal, int NJACK, const char *name_in, coord_t coord_input, error **err);
void free_galaxy_data(galaxy_data *galaxy);

double get_coordinate_unit(coord_t coord, error **err);
void coordinate_to_rad(double *pos, coord_t coord, error **err);
void rad_to_coordinate(double *pos, coord_t coord, error **err);

void node_fill_galaxy_common(node *nd, const void *tmp_data, uint start, uint end, int ndim, int Nboots);
void node_fill_galaxy_xy(node *nd, const void *tmp_data, uint start, uint end, int ndim, int Nboots);
void node_fill_galaxy_radec(node *nd, const void *tmp_data, uint start, uint end, int ndim, int Nboots);
double get_pos_galaxy(void *tmp_data, uint i, int dim);
void copy_galaxy_info(galaxy_data *a, const galaxy_data *b, int ndim);
void swap_galaxy(void *tmp_data, uint i, uint j, int ndim);

double *init_ngal_boots(const galaxy_data *data, int Nboots, uint start, uint end);
dcomplex *init_well_resample(const galaxy_data *data, int Nresample, uint start, uint end, double **weight_resample);
int shift_to_branch(double alpha);

void get_extend_2d(const galaxy_data *gal, uint ngal, const galaxy_data *gal2, uint ngal2, double *min, double *max);

void get_col(const char *col_name, col_descr *col, error **err);


#ifdef _WITH_FITS
galaxy_data *read_gal_fits(const char *name, uint *ngal, int radec, coord_t coord_input,
                      int check_zero_weight, int Ndim, const char **col_names, int Ncol, int quiet, error **err);
#endif



#endif
