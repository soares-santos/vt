/* ============================================================ *
 * gal_cat.c							                               *
 * Martin Kilbinger 2006-2012					                      *
 * ============================================================ */

#include "gal_cat.h"


/* ============================================================ *
 * Reads galaxy data from a file and returns an array.		    *
 * ============================================================ */
galaxy_data *read_gal(const char *name, uint *ngal, int radec, format_t format, coord_t coord_input, 
                      int check_zero_weight, int ndim, const char **col_names, int Ncol, int quiet, error **err)
{
   FILE *F;
   uint i, j;
   size_t nmem;
   int n, idx;
   galaxy_data *gal;
   double eps[2], dummy[5], w;
   char str[1024], *res;

   testErrorRet(ndim != 2 && format != f_position, gc_dim, "Dimension!=2 only implemented for format=f_position", *err, __LINE__, NULL);


   if (format == f_fits) {

#ifdef _WITH_FITS
      gal = read_gal_fits(name, ngal, radec, coord_input, 
                          check_zero_weight, ndim, col_names, Ncol, quiet, err);
      forwardError(*err, __LINE__, NULL);

      post_process_galaxy(gal, *ngal, coord_input, radec, check_zero_weight, ndim, err);
      forwardError(*err, __LINE__, NULL);

      return gal;
#else
      *err = addError(gc_fits, "Re-compile 'athena' with fits support", *err, __LINE__);
      return gal;
#endif

   }

   F = fopen_err(name, "r", err); forwardError(*err, __LINE__, NULL);

   /* Read first line to see whether it contains number of objects */
   res = fgets(str, 1024, F);
   testErrorRetVA(res==NULL, gc_io, "Error while reading file %s", *err, __LINE__, NULL, name);
   n = sscanf(str, "%lf %lf\n", dummy, dummy+1);

   if (n>1) {
      /* More than one number in first line: no header */
      fclose(F);
      if (!quiet) {
         fprintf(stderr, "Counting lines..."); fflush(stderr);
      }
      *ngal = numberoflines(name, err);
      forwardError(*err, __LINE__, NULL);
      if (!quiet) fprintf(stderr, "done\n");
      F = fopen(name, "r");
   } else if (n==1) {
      /* First line is header */
      sscanf(str, "%u", ngal);
   } else {
      *err = addErrorVA(gc_io, "Error while reading first line of file %s", *err, __LINE__, name);
      return NULL;
   }

   if (!quiet) fprintf(stderr, "%u galaxies\n", *ngal);

   testErrorRet(*ngal<=0, gc_ngal, "Number of galaxies has to be positive", *err, __LINE__, NULL);
   nmem = sizeof(galaxy_data)*(*ngal);
   gal = malloc_err(nmem, err);   forwardError(*err, __LINE__, NULL);
   gal->ind_resample = NULL;

   for (i=0; i<(*ngal); i++) {
      testErrorRet(feof(F), gc_io, "Premature eof encountered", *err, __LINE__, NULL);

      /* Default values */
      w = 1.0;
      eps[0] = eps[1] = 0.0;
      idx = -1;

      switch (format) {
         case f_standard :
            /* Format: x y g1 g2 w */
            n = fscanf(F, "%lg %lg %lg %lg %lg\n", gal[i].pos, gal[i].pos+1, eps, eps+1, &w);
            testErrorRetVA(n != 5, gc_io, "%d elements read from line %d, expected %d", *err, __LINE__, NULL, n, i, 5);
            break;

         case f_hamana :
            /* Format: x y z k g1 g2 */
            n = fscanf(F, "%lg %lg %lg %lg %lg %lg\n", gal[i].pos, gal[i].pos+1, dummy, dummy+1, eps, eps+1);
            testErrorRetVA(n != 6, gc_io, "%d elements read from line %d, expected %d, for galaxy at pos=(%g,%g)",
                  *err, __LINE__, NULL, n, i, 6, gal[i].pos[0], gal[i].pos[1]);
            break;

         case f_position :
            /* Format: x0 x1 [x2 ...] w */
            for (j=0; j<ndim; j++) {
               n = fscanf(F, "%lg", gal[i].pos + j);
               testErrorRetVA(n != 1, gc_io, "%d elements read from line %d, expected %d", *err, __LINE__, NULL, n, i, 1);
            }
            n = fscanf(F, "%lg\n", &w);
            testErrorRetVA(n != 1, gc_io, "%d elements read from line %d, expected %d", *err, __LINE__, NULL, n, i, 1);
            break;

         case f_position_jack_num :
            /* Format: x y w num */
            n = fscanf(F, "%lg %lg %lg %d\n", gal[i].pos, gal[i].pos+1, &w, &idx);
            testErrorRetVA(n != 4, gc_io, "%d elements read from column %d, expected %d", *err, __LINE__, NULL, n, i, 4);
            testErrorRetVA(idx<0 || idx>=USHRT_MAX, gc_overflow, "Jackknife index %d out of range [0; %d]", *err, __LINE__, NULL,
                           idx, USHRT_MAX); 
            break;

         case f_lensing_jack_num :
            /* Format: x y g1 g2 w num */
            n = fscanf(F, "%lg %lg %lg %lg %lg %d\n", gal[i].pos, gal[i].pos+1, eps, eps+1, &w, &idx);
            testErrorRetVA(n != 6, gc_io, "%d elements read from line %d, expected %d", *err, __LINE__, NULL, n, i, 6);
            break;

         default :
            *err = addErrorVA(gc_type, "Wrong catalogue format type %d", *err, __LINE__, format);
            return NULL;
      }

      gal[i].weight = (fdreal)w;
      gal[i].idx    = (unsigned short int)idx;

      /* Complex ellipticity */
      gal[i].ell = FDComplex((fdreal)eps[0], (fdreal)eps[1]);

   }

   fileclose(F);

   post_process_galaxy(gal, *ngal, coord_input, radec, check_zero_weight, ndim, err);
   forwardError(*err, __LINE__, NULL);

   return gal;
}

/* ============================================================ *
 * Fills in galaxy catalogue with missing pieces after reading  *
 * from file.                                                   *
 * ============================================================ */
void post_process_galaxy(galaxy_data *gal, uint ngal, coord_t coord_input, int radec, int check_zero_weight, int ndim, error **err)
{
   uint i, j;

   for (i=0; i<ngal; i++) {
      /* Change coordinates to radians */
      for (j=0; j<ndim; j++) {
         coordinate_to_rad(gal[i].pos + j, coord_input, err);   forwardError(*err, __LINE__,);
      }

      if (radec == 1) {
         for (j=0; j<ndim; j++) {
            gal[i].cos[j] = cos(gal[i].pos[j]);
            gal[i].sin[j] = sin(gal[i].pos[j]);
         }
      }

      if (check_zero_weight) {
         /* MKDEBUG: Zero weights cause failure of the tree code */
         testErrorRetVA(gal[i].weight<1.0e-10, gc_cat_entry,
               "Dubious entry with very small or negative weight %g in galaxy catalog found!\n"
               "gal[%d]: pos=(%g, %g), eps=(%g,%g)",
               *err, __LINE__,, gal[i].weight, i, gal[i].pos[0], gal[i].pos[1], gal[i].ell.r, gal[i].ell.i);
      }

   }

}

/* ============================================================ *
 * Prints galaxy #i to output file.                             *
 * ============================================================ */
void out_galaxy(const void *data, uint i, int ndim, FILE *OUT)
{
   galaxy_data *galaxy;
   int j;

   galaxy = ((galaxy_data*)data)+i;

   if (OUT==NULL) OUT = stderr;
   fprintf(OUT, "gal = (");
   for (j=0; j<ndim; j++) {
      fprintf(OUT, "%f ", galaxy->pos[j]);
   }
   fprintf(OUT, ") (%g %g) %g\n", galaxy->ell.r, galaxy->ell.i, galaxy->weight);
}

/* ============================================================ *
 * Writes the galaxy catalogue with data 'gal' in Jackknife     *
 * sample list mode.						*
 * ============================================================ */
void out_cat_jack(const galaxy_data *gal, uint ngal, int NJACK, const char *name_in, coord_t coord_input, error **err)
{
   char name_out[1024];
   FILE *F;
   uint j;
   int b;
   double xy[2];

   sprintf(name_out, "%s.jack_list", name_in);
   F = fopen_err(name_out, "w", err);
   forwardError(*err, __LINE__,);

   for (j=0; j<ngal; j++) {
      for (b=0; b<2; b++) {
	 xy[b] = gal[j].pos[b];
	 rad_to_coordinate(xy+b, coord_input, err);
	 forwardError(*err, __LINE__,);
      }
      fprintf(F, "% .10g % .10g % g", xy[0], xy[1], gal[j].weight);
      for (b=0; b<NJACK; b++) {
         fprintf(F, " %.0f", gal[j].ind_resample[b]);
      }
      fprintf(F, "\n");
   }
   fclose(F);
}

void free_galaxy_data(galaxy_data *galaxy)
{
   if (galaxy->ind_resample) free(galaxy->ind_resample);
   free(galaxy);
}

double get_coordinate_unit(coord_t coord, error **err)
{
   switch (coord) {

      case coord_arcsec : return ARCSEC; break;
      case coord_arcmin : return ARCMIN; break;
      case coord_deg    : return pi/180.0; break;
      case coord_rad    : return 1.0;
      default           : *err = addErrorVA(gc_coord, "Unknown coordinate type (%d)", *err, __LINE__, coord);
     			              return -1.0;
   }
}

void coordinate_to_rad(double *pos, coord_t coord, error **err)
{
   *pos *= get_coordinate_unit(coord, err);
   forwardError(*err, __LINE__,);
}

void rad_to_coordinate(double *pos, coord_t coord, error **err)
{
   *pos /= get_coordinate_unit(coord, err);
   forwardError(*err, __LINE__,);
}

/* ============================================================= *
 * Fills information from tmp_data (pointer to galaxy_data) into *
 * node (used by grow_tree).					 *
 * ============================================================= */
void node_fill_galaxy_common(node *nd, const void *tmp_data, uint start, uint end, int ndim, int Nresample)
{
   const galaxy_data *data;
   uint i, k;

   data = (galaxy_data*)tmp_data;

   /* Extension of node. Works for RADEC if coordinates do not go through alpha=0 */
   for (k=0; k<ndim; k++) {
      nd->min[k] = 1.0e30;
      nd->max[k] = -1.0e30;
   }

   nd->well.r = nd->well.i = nd->weight = 0.0;
   for (i=start; i<end; i++) {

      for (k=0; k<ndim; k++) {
         if (data[i].pos[k]<nd->min[k]) nd->min[k] = data[i].pos[k];
         if (data[i].pos[k]>=nd->max[k]) nd->max[k] = data[i].pos[k];
      }

      /* Sum of weighted ellipticities */
      nd->well.r += data[i].ell.r*data[i].weight;
      nd->well.i += data[i].ell.i*data[i].weight;

      /* Sum of weights */
      nd->weight += data[i].weight;

   }

   assert(nd->ngal>=1);

   nd->ngal_resample = init_ngal_boots(data, Nresample, start, end);
   nd->well_resample = init_well_resample(data, Nresample, start, end, &nd->weight_resample);
}

/* See node_fill_galaxy */
void node_fill_galaxy_xy(node *nd, const void *tmp_data, uint start, uint end, int ndim, int Nboots)
{
   galaxy_data *data;
   uint i, k;
   double dsqr, rsqr;

   node_fill_galaxy_common(nd, tmp_data, start, end, ndim, Nboots);

   data = (galaxy_data*)tmp_data;

   /* Node barycenter, weighted arithmetic mean of weights */
   for (k=0; k<ndim; k++) nd->barycenter[k] = 0.0;
   for (i=start; i<end; i++) {
      for (k=0; k<ndim; k++) {
         nd->barycenter[k] += data[i].pos[k] * data[i].weight;
      }
   }
   for (k=0; k<ndim; k++) {
      nd->barycenter[k] /= nd->weight;
   }

   /* Node extension */
   for (i=start,rsqr=0.0; i<end; i++) {
      dsqr = absnsqr(nd->barycenter, data[i].pos, ndim);
      if (dsqr>rsqr) rsqr = dsqr;
   }
   nd->radius = sqrt(rsqr);
}

/* See node_fill_galaxy */
void node_fill_galaxy_radec(node *nd, const void *tmp_data, uint start, uint end, int ndim, int Nboots)
{
   galaxy_data *data;
   uint i, k;
   double dsqr, rsqr, rb[3], rbabs, tmp, alpha_mid;
   int n_shift;

   if (ndim > 2) {
      fprintf(stderr, "Cannot calculate barycenter on sphere of dimension ndim>2.");
      fprintf(stderr, "Mode 'radec' not possible with more than two coordinates.");
      assert(0);
   }

   node_fill_galaxy_common(nd, tmp_data, start, end, ndim, Nboots);

   data = (galaxy_data*)tmp_data;

   /* Node barycenter, weighted arithmetic mean of weights */

   /* Calculate R^3 barycenter, rescale to unit vector later */
   for (k=0; k<3; k++) rb[k] = 0.0;
   for (i=start; i<end; i++) {  
      rb[0] += data[i].cos[0]*data[i].cos[1] * data[i].weight;      /* cos alpha cos delta */
      rb[1] += data[i].sin[0]*data[i].cos[1] * data[i].weight;      /* sin alpha cos delta */
      rb[2] += data[i].sin[1] * data[i].weight;                     /* sin delta           */
   }

   for (k=0; k<3; k++) rb[k] /= nd->weight;

   /* Extend barycenter vector to lie on sphere */
   rbabs = sqrt(rb[0]*rb[0] +  rb[1]*rb[1] + rb[2]*rb[2]);
   for (k=0; k<3; k++) rb[k] /= rbabs;

   nd->sinbarycenter[1] = rb[2];                         /* sin delta */
   nd->barycenter[1]    = asin(nd->sinbarycenter[1]);    /* delta */
   nd->cosbarycenter[1] = cos(nd->barycenter[1]);        /* cos delta */
   nd->cosbarycenter[0] = rb[0]/nd->cosbarycenter[1];    /* cos alpha */
   nd->sinbarycenter[0] = rb[1]/nd->cosbarycenter[1];    /* sin alpha */

   alpha_mid = (nd->max[0]+nd->min[0])/2.0;              /* alpha_midpoint */
   n_shift = shift_to_branch(alpha_mid);                 /* Shift to [0;pi) branch */
   tmp = nd->cosbarycenter[0];
   if (n_shift%2!=0) tmp *= -1.0;                        /* cos(n_shift*pi) */
   nd->barycenter[0] = acos(tmp)-n_shift*pi;             /* alpha */

   /* Node extension */
   for (i=start,rsqr=0.0; i<end; i++) {
      dsqr = abs2sqr(nd->barycenter, data[i].pos);
      if (dsqr>rsqr) rsqr = dsqr;
   }
   nd->radius = sqrt(rsqr);
}

/* Returns the xy-component of position of i-th galaxy. tmp_data points to galaxy_data. *
 * Used in partition_data								*/
double get_pos_galaxy(void *tmp_data, uint i, int dim)
{
   double pos;
   galaxy_data *data;

   data = (galaxy_data*)tmp_data;
   pos = data[i].pos[dim];

   return pos;
}

/* Copies all entries from b to a */
void copy_galaxy_info(galaxy_data *a, const galaxy_data *b, int ndim)
{
   int k;

   for (k=0; k<ndim; k++) {
      a->pos[k]    = b->pos[k];
      a->cos[k]    = b->cos[k];
      a->sin[k]    = b->sin[k];
   }

   a->ell.r        = b->ell.r;
   a->ell.i        = b->ell.i;
   a->weight       = b->weight;
   a->idx          = b->idx;
   //a->z            = b->z;
   //a->dzl          = b->dzl;
   //a->dzu          = b->dzu;
   a->ind_resample = b->ind_resample;
}

/* Swaps the information of the two galaxies i and j (used in partition_data) */
void swap_galaxy(void *tmp_data, uint i, uint j, int ndim)
{
   galaxy_data *dummy, *a, *b, temp;

   dummy = (galaxy_data*)tmp_data;
   a = &dummy[i];
   b = &dummy[j];

   copy_galaxy_info(&temp, a, ndim);
   copy_galaxy_info(a, b, ndim);
   copy_galaxy_info(b, &temp, ndim);
}

/* Returns an array of the resample numbers for the node [start; end) */
double *init_ngal_boots(const galaxy_data *data, int Nresample, uint start, uint end)
{
   int i, j;
   double *ngal_resample;

   if (Nresample>0) {

      ngal_resample = calloc(Nresample, sizeof(double));
      assert(ngal_resample);
      for (i=0; i<Nresample; i++) {
         for (j=start; j<end; j++) {
            /* Adds multiplicity of galaxy number (or weight) */
            ngal_resample[i] += data[j].ind_resample[i];
         }
      }

   } else {

      ngal_resample = NULL;

   }

   return ngal_resample;
}

/* Returns an array of the resample weighted ellipticities for the node [start; end) */
dcomplex *init_well_resample(const galaxy_data *data, int Nresample, uint start, uint end, double **weight_resample)
{
   int i, j;
   dcomplex *well_resample;

   if (Nresample > 0) {

      *weight_resample = (double*)calloc(Nresample, sizeof(double));
      well_resample    = (dcomplex*)calloc(Nresample, sizeof(dcomplex));
      assert(well_resample);
      for (i=0; i<Nresample; i++) {
         for (j=start; j<end; j++) {
            well_resample[i].r += data[j].ell.r * data[j].weight * data[j].ind_resample[i];
            well_resample[i].i += data[j].ell.i * data[j].weight * data[j].ind_resample[i];
            (*weight_resample)[i] += data[j].weight * data[j].ind_resample[i];
         }
      }

   } else {

      *weight_resample = NULL;
      well_resample    = NULL;

   }

   return well_resample;
}

/* ============================================================ *
 * Returns n (positive or negative) such that			*
 * alpha + n*pi is in the 'good' branch of acos [0:pi)		*
 * ============================================================ */
int shift_to_branch(double alpha)
{
   int n = 0;

   if (alpha>=0 && alpha<pi) {
      return 0;
   }
   /* MKDEBUG: New v1.5 (hint from Shahab Joudaki) */
   if (alpha>=pi) {
      while (alpha+n*pi>pi) n--;
      return n;
   }
   if (alpha<0) {
      while (alpha+n*pi<0) n++;
      return n;
   }
   assert(0);
   return -1;
}

/* ============================================================ *
 * Calculates extend (min and max) for the first 2 coordinates  *
 * from catalogue gal. If gal2 != NULL, takes the max of both   *
 * catalogues.							                               *
 * ============================================================ */
#define BIG 1e30
void get_extend_2d(const galaxy_data *gal, uint ngal, const galaxy_data *gal2, uint ngal2, double *min, double *max)
{
   int xy;
   uint j;

   min[0] = min[1] = BIG;
   max[0] = max[1] = -BIG;

   for (j=0; j<ngal; j++) {
      for (xy=0; xy<2; xy++) {
         if (gal[j].pos[xy] < min[xy]) min[xy] = gal[j].pos[xy];
         if (gal[j].pos[xy] > max[xy]) max[xy] = gal[j].pos[xy];
      }
   }

   if (gal2 != NULL) {
      for (j=0; j<ngal2; j++) {
         for (xy=0; xy<2; xy++) {
            if (gal2[j].pos[xy] < min[xy]) min[xy] = gal2[j].pos[xy];
            if (gal2[j].pos[xy] > max[xy]) max[xy] = gal2[j].pos[xy];
         }
      }
   }
}
#undef BIG


/* ============================================================ *
 * From a column description string of format "key:val" get     *
 * column key and value.                                        *
 * ============================================================ */
void get_col(const char *col_name, col_descr *col, error **err)
{
   char *val, *skey, col_name_copy[128], *brkb;
   int j;

   strcpy(col_name_copy, col_name);

   skey  = strtok_r(col_name_copy, COL_NAME_SEP, &brkb); 
   STRING2ENUM(col->key, skey, col_names_t, scol_names_t, j, Ncol_names_t, err);

   val = strtok_r(NULL, COL_NAME_SEP, &brkb); 
   strcpy(col->val, val);
}

#ifdef _WITH_FITS
galaxy_data *read_gal_fits(const char *name, uint *ngal, int radec, coord_t coord_input, 
                           int check_zero_weight, int Ndim, const char **col_names, int Ncol, int quiet, error **err)
{
   fitsfile *fptr;
   int status = 0, res, fits_col_num, anynull, hdu, nhdu, hdutype;
   long n;
   uint i, j;
   size_t nmem;
   galaxy_data *gal;
   col_descr col;
   double *array;


   res = fits_open_image(&fptr, name, READONLY, &status);
   FITS_STATUS_ERROR(status, quiet, "Error (status=%d) while opening fits file %s", *err, __LINE__, NULL, status, name);
   testErrorRetVA(res != 0, gc_fits, "Error (return value=%d) while opening fits file %s", *err, __LINE__, NULL, res, name);

   fits_get_num_hdus(fptr, &nhdu, &status);
   FITS_STATUS_ERROR(status, quiet, "Error (status=%d) while reading number of HDUs", *err, __LINE__, NULL, status);

   for (hdu=1; hdu<=nhdu; hdu++) {
      fits_movabs_hdu(fptr, hdu, &hdutype, &status);
      FITS_STATUS_ERROR(status, quiet, "Error (status=%d) while moving to HDU %d", *err, __LINE__, NULL, status, hdu);
      if (hdutype == BINARY_TBL || hdutype == ASCII_TBL) break;
   }
   testErrorRetVA(!(hdutype == BINARY_TBL || hdutype == ASCII_TBL), gc_fits, "No Table in fits file found, #hdu=%d",
                  *err, __LINE__, NULL, nhdu);

   fits_get_num_rows(fptr, &n, &status);
   FITS_STATUS_ERROR(status, quiet, "Error (status=%d) while reading number of rows in fits file", *err, __LINE__, NULL, status);

   *ngal = (uint)n;

   nmem = sizeof(galaxy_data) * (*ngal);
   gal = malloc_err(nmem, err);   forwardError(*err, __LINE__, NULL);
   gal->ind_resample = NULL;
   array = malloc_err(sizeof(double) * (*ngal), err);   forwardError(*err, __LINE__, NULL); 

   /* Default values */
   for (i=0; i<(*ngal); i++) {
      gal[i].weight = 1.0;
      gal[i].ell.r = gal[i].ell.i = 0.0;
      gal[i].idx = -1;
   }


   /* Read all columns with names as indicated in the config file */
   for (j=0; j<Ncol; j++) {

      /* Find fits column number according to column name*/
      get_col(col_names[j], &col, err);
      forwardError(*err, __LINE__, NULL);
      fits_get_colnum(fptr, CASEINSEN, col.val, &fits_col_num, &status);
      FITS_STATUS_ERROR(status, quiet, "Error (status=%d) while getting fits column number with (s)key, val) = (%d)%s, %s\n",
                        *err, __LINE__, NULL, status, col.key, scol_names_t(col.key), col.val);

      /* Read column */
      fits_read_col(fptr, tcol_fits_types_t(col.key), fits_col_num, 1, 1, n, NULL, array, &anynull, &status);
      FITS_STATUS_ERROR(status, quiet, "Error (status=%d) while reading column with (s)key, val) = (%d)%s, %s\n",
                        *err, __LINE__, NULL, status, col.key, scol_names_t(col.key), col.val);

      for (i=0; i<(*ngal); i++) {

         switch (col.key) {
            case c_x   : gal[i].pos[0] = array[i];
                         break;
            case c_y   : gal[i].pos[1] = array[i];
                         break;
            case c_e1  : gal[i].ell.r = array[i];
                         break;
            case c_e2  : gal[i].ell.i = array[i];
                         break;
            case c_w   : gal[i].weight = array[i];
                         break;
            case c_njk : gal[i].idx = (int)array[i];
                         break;
            default :
               *err = addErrorVA(gc_fits, "Unkown column (s)key (%d)%s", *err, __LINE__, col.key, scol_names_t(col.key));
               return NULL;
         }
      }

   }

   free(array);

   return gal;
}


#endif

