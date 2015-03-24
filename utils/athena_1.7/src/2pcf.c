/* ============================================================ *
 * Martin Kilbinger 2006-2012				        *
 * 2pcf.c							*
 * Common 2pcf-related functions for brute-force and tree.      *
 * ============================================================ */


#include "2pcf.h"


/* ============================================================ *
 * Reads galaxy data from ncat files and sets the arrays gal1   *
 * and gal2.							*
 * File format: first line ngal, all other lines x y e1 e2 w.   *
 * x, y in arcsec/rad for radec5=0/>1.				*
 * ============================================================ */
void read_all_gal(const char *name1, const char *name2, uint *ngal1, uint *ngal2, int ncat, int radec, format_t format, coord_t coord_input,
		  galaxy_data **gal1, galaxy_data **gal2, int ndim, const char **col_names, int Ncol, int quiet, error **err)
{
   /* Read first catalogue */
   *gal1 = read_gal(name1, ngal1, radec, format, coord_input, 1, ndim, col_names, Ncol, quiet, err);
   forwardError(*err, __LINE__,);

   if (ncat==1) {

      /* Auto-correlation, copy catalogue */
      *ngal2 = *ngal1;

      /* NEW: don't copy catalogue to 2nd, only first cat is swapped during partition to build first tree */
      *gal2  = *gal1;

   } else if (ncat==2) {

      /* Cross-correlation, read second catalogue */
      *gal2 = read_gal(name2, ngal2, radec, format, coord_input, 1, ndim, col_names, Ncol, quiet, err);
      forwardError(*err, __LINE__,);

   } else {
      *err = addError(tpcf_ncat, "Number of catalogues has to be 1 or 2", *err, __LINE__);
      return;
   }

}

/* ============================================================ *
 * Initialises and returns a bin structure.			*
 * ============================================================ */
bin_data *init_bin(const double *thmin, const double *thmax, const uint *nth, int ndimbin, bin_type I_BINTYPE,
		   error **err)
{
   bin_data *bin;
   int k;

   bin = malloc_err(sizeof(bin_data), err);
   forwardError(*err, __LINE__, NULL);

   for (k=0; k<ndimbin; k++) {
      testErrorRet(thmin[k] < 0 || thmax[k] < 0 || thmax[k] <= thmin[k], tpcf_outofrange,
		"Angular scales (in config file) out of range", *err, __LINE__, NULL);

      bin->min[k]     = thmin[k];
      bin->max[k]     = thmax[k];
      bin->minsqr[k]  = DSQR(bin->min[k]);
      bin->maxsqr[k]  = DSQR(bin->max[k]);
      bin->logmin[k]  = log(bin->min[k]);
      bin->logmax[k]  = log(bin->max[k]);
      bin->diff[k]    = bin->max[k] - bin->min[k];
      bin->difflog[k] = bin->logmax[k] - bin->logmin[k];
      bin->N[k]       = nth[k];
      bin->linlog  = I_BINTYPE;
   }

   bin->ndim = ndimbin;

   testErrorRetVA(bin->ndim > NDIMBIN_MAX, tpcf_ndim, "Bin dimension %d out of range, > %d",
		  *err, __LINE__, NULL, bin->ndim, NDIMBIN_MAX);

   testErrorRet(bin->linlog!=LIN && bin->linlog!=LOG, tpcf_bintype,
		"BINTYPE (in config file) has to be 'LIN' or 'LOG'", *err, __LINE__, NULL);

   return bin;
}

/* Returns the total number of bins, which is the product of the  *
 * number of bins in each dimension.				  */
uint bin_N(const bin_data *bin)
{
   uint N, k;

   for (k=0,N=1; k<bin->ndim; k++) {
      N *= bin->N[k];
   }

   return N;
}

/* Returns a bin index corresponding to distance */
int distance2bin_index(double distance, double min, double diff, int N, int round)
{
   int i;
   /* N is correct, not N-1! */
   i = (int)(N * (distance - min) / diff);

   /* New in v1.61:                               *
    * Due to rounding errors it's possible that   *
    * distance < max, and i = N. In that case, if *
    * desired (round=1), i is set to N-1.         */
   if (i==N && round) i = N-1;

   return i;
}

/* ============================================================ *
 * Returns the index for a ndim-array with dimensions N,        *
 * encoded in an 1d-array.					*
 * Returns b0 + N0 * (b1 + N1 * (b2 + N2 * (b3 + ... )))	*
 * ============================================================ */
int get_bin_index(const int *b, const uint *N, int ndim)
{
   int k, b_index;

   for (k=ndim-1,b_index=0; k>=0; k--) {
      b_index = b_index * N[k-1] + b[k];
   }

   return b_index;
}

double bin_index_to_scale(const bin_data *bin, const double *wtheta, int b)
{
   double theta;

   if (wtheta == NULL || wtheta[b] < 0) {
      if (bin->linlog==LIN) {
         theta = (double)(b+0.5)/(double)bin->N[0] * bin->diff[0] + bin->min[0];
      } else {
         theta = exp((double)(b+0.5)/(double)bin->N[0] * bin->difflog[0] + bin->logmin[0]);
      }
   } else {
      theta = wtheta[b];
   }

   return theta;
}

double bin_index_to_scale_low(const bin_data *bin, const double *wtheta, int b)
{
  double theta1;
  
  if (wtheta == NULL || wtheta[b] < 0) {
     if (bin->linlog==LIN) {
        theta1 = (double)b/(double)bin->N[0] * bin->diff[0] + bin->min[0];
     } else {
        theta1 = exp((double)b/(double)bin->N[0] * bin->difflog[0] + bin->logmin[0]);
     }
  } else {
     theta1 = wtheta[b];
  }

  return theta1;
}

double bin_index_to_scale_high(const bin_data *bin, const double *wtheta, int b)
{
  double theta2;
  
  if (wtheta == NULL || wtheta[b] < 0) {
     if (bin->linlog==LIN) {
        theta2 = (double)(b+1)/(double)bin->N[0] * bin->diff[0] + bin->min[0];
     } else {
        theta2 = exp((double)(b+1)/(double)bin->N[0] * bin->difflog[0] + bin->logmin[0]);
     }
  } else {
     theta2 = wtheta[b];
  }

  return theta2;   
}

/* Shape noise. See Schneider, van Waerbeke, Kilbinger & Mellier (2001), eq. (23) for D */
double calculate_sqrt_D(double sig_eps4, double wwww, double ww, int ncat)
{
   double D;

   if (ww > 0) {
      D  = sig_eps4 * wwww / DSQR(ww);
   } else {
      D = 0;
   }
   if (D<0) D = 0;  // MKDEBUG: Does this case every occur?

   /* Correction for sum_i<j -> sum_ij over pairs*/
   if (ncat==1) {
      D /= 2;
   }

   return sqrt(D);
}

/* Approximate shape noise correction for open angle */
double calculate_sqrt_Dcor(double sqrtD, unsigned long long int npair, int Nnode)
{
   double Dcor;

   if (sqrtD>0) {
      Dcor = sqrtD * sqrtD / ((double)npair / (double)Nnode);
   } else {
      Dcor = 0.0;
   }

   return sqrt(Dcor);
}

/* ============================================================ *
 * Initialises and returns a 'xi' (correlation data) structure. *
 * ============================================================ */
xi_data *init_xi(const bin_data *bin, wcorr_t wcorr, int Nresample, int do_wtheta, error **err)
{
   int b;
   uint binNN;
   xi_data *xi;

   xi = malloc_err(sizeof(xi_data), err);
   forwardError(*err, __LINE__, NULL);

   binNN = bin_N(bin);
   xi->npair = calloc_err(binNN, sizeof(unsigned long long int), err);
   forwardError(*err, __LINE__, NULL);

   xi->p = xi->m = xi->x = xi->ww = xi->wn = xi->wwww = xi->wnwn = xi->gt = xi->gx = xi->g11 = xi->g12 = xi->g21 = xi->g22 = NULL;

   xi->npair_resample = xi->p_resample = xi->m_resample = xi->ww_resample = xi->wn_resample = xi->gt_resample = NULL;

   if (wcorr&gg) {
      xi->p    = vector(0, binNN - 1);
      xi->m    = vector(0, binNN - 1);
      xi->x    = vector(0, binNN - 1);
      xi->ww   = vector(0, binNN - 1);
      xi->wwww = vector(0, binNN - 1);
      xi->g11  = vector(0, binNN - 1);
      xi->g12  = vector(0, binNN - 1);
      xi->g21  = vector(0, binNN - 1);
      xi->g22  = vector(0, binNN  -1);
   }

   if (wcorr&gn) {
      xi->gt   = vector(0, binNN - 1);
      xi->gx   = vector(0, binNN - 1);
   }

   if (wcorr&gn || wcorr&wn) {
      xi->wn   = vector(0, binNN - 1);
      xi->wnwn = vector(0, binNN - 1);
   }
   if (wcorr&gg || wcorr&gn) {
      xi->Nnode = ivector(0, binNN - 1);
   }

   /* Resamples (bootstrap or jackknife) */
   if (wcorr&nn || wcorr&wn) {
      xi->npair_resample = calloc_array_double(binNN, Nresample, err);
      forwardError(*err, __LINE__, NULL);
   }
   if (wcorr&gg) {
      xi->p_resample = calloc_array_double(binNN, Nresample, err);
      forwardError(*err, __LINE__, NULL);
      xi->m_resample = calloc_array_double(binNN, Nresample, err);
      forwardError(*err, __LINE__, NULL);
      xi->ww_resample = calloc_array_double(binNN, Nresample, err);
      forwardError(*err, __LINE__, NULL);
   }
   if (wcorr&gn) {
      xi->gt_resample = calloc_array_double(binNN, Nresample, err);
      forwardError(*err, __LINE__, NULL);
   }
   if (wcorr&gn || wcorr&wn) {
      xi->wn_resample = calloc_array_double(binNN, Nresample, err);
      forwardError(*err, __LINE__, NULL);
   }

   for (b=0; b<binNN; b++) {

      if (wcorr&gg) {
         xi->p[b] = xi->m[b] = xi->x[b] = xi->ww[b] = xi->wwww[b] = xi->g11[b] = xi->g12[b] = xi->g21[b] = xi->g22[b] = 0.0;
      }
      if (wcorr&gn) {
         xi->gt[b] = xi->gx[b] = xi->wnwn[b] = 0.0;
      }
      if (wcorr&gn || wcorr&wn) {
         xi->wn[b] = 0.0;
      }
      if (wcorr&gg || wcorr&gn) {
         xi->Nnode[b] = 0;
      }

   }

   if (do_wtheta) {
      xi->wtheta = malloc_err(binNN * sizeof(double), err);
      forwardError(*err, __LINE__, NULL);
   } else {
      xi->wtheta = NULL;
   }

   return xi;
}

/* Normalizes xi by the galaxy pair weights for each bin */
void weigh(const bin_data *bin, xi_data *xi, wcorr_t wcorr, int nresample, error **err)
{
   int b, i;
   double w;

   for (b=0; b<bin_N(bin); b++) {

      if (wcorr&gg && xi->npair[b]) {

         testErrorRetVA(xi->ww[b]<=0, tpcf_npositive, "Weight %g is not positive for bin %d",
               *err, __LINE__,, xi->ww[b], b);

         xi->p[b]   /= xi->ww[b];
         xi->m[b]   /= xi->ww[b];
         xi->x[b]   /= xi->ww[b];
         xi->g11[b] /= xi->ww[b];
         xi->g12[b] /= xi->ww[b];
         xi->g21[b] /= xi->ww[b];
         xi->g22[b] /= xi->ww[b];

         if (xi->p_resample != NULL && xi->m_resample != NULL) {

            for (i=0; i<nresample; i++) {

               if (xi->ww_resample[b][i] >= 0) {
                  xi->p_resample[b][i] /= xi->ww_resample[b][i];
                  xi->m_resample[b][i] /= xi->ww_resample[b][i];
               }

            }

         }

      }

      if (wcorr&gn && xi->npair[b] ) {  /* Normalizes gammat and gammax */

         testErrorRetVA(xi->wn[b]<=0, tpcf_npositive, "Weight %g is not positive for bin %d",
               *err, __LINE__,, xi->wn[b], b);

         xi->gt[b] /= xi->wn[b];
         xi->gx[b] /= xi->wn[b];

         if (xi->gt_resample != NULL) {
            for (i=0; i<nresample; i++) {
               if (xi->wn_resample[b][i] >= 0) {
                  xi->gt_resample[b][i] /= xi->wn_resample[b][i];
               }
            }
         }

      }

      if (wcorr&wn && xi->npair[b]) {
         xi->wn[b] /= (double)xi->npair[b];
      }

      if (xi->wtheta != NULL) {
         // MKDEBUG: If gg and gn are run jointly, gg weight is used for weighed angular scale.
         if (wcorr & gg) {
            w = xi->ww[b];
         } else if (wcorr & gn) {
            w = xi->wn[b];
         } else {
            w = -1;
         }

         if (w > 0) {
            xi->wtheta[b] /= w;
         } else {
            xi->wtheta[b] = -1.0;
         }
      }

   }
}

double sigma_epsilon_sqr(const galaxy_data *gal, uint ngal)
{
   uint i;
   double sig_eps_sqr[2] = {0.0, 0.0};

   for (i=0; i<ngal; i++) {
      sig_eps_sqr[0] += DSQR(gal[i].ell.r);
      sig_eps_sqr[1] += DSQR(gal[i].ell.i);
   }
   sig_eps_sqr[0] /= (double)ngal;
   sig_eps_sqr[1] /= (double)ngal;

   return sig_eps_sqr[0] + sig_eps_sqr[1];
}

/* ============================================================ *
 * Returns the resample rms of x_resample[b] for bin b. If      *
 * mean_resample is not NULL, sets the mean.                    *
 * Valid for bootstrap and Jackknife.                           *
 * ============================================================ */
double get_resample_error(const double **x_resample, int b, int nresample, error_t ERROR, double *mean_resample, error **err)
{
   int i;
   double x, xsqr, rms, dnresample;

   x = xsqr = 0.0;

   for (i=0; i<nresample; i++) {
      x    += x_resample[b][i];
      xsqr += DSQR(x_resample[b][i]);
   }

   testErrorRetVA(nresample < 0, tpcf_npositive, "Nresample = %d has to be larger than zero", *err, __LINE__, -1.0, nresample);
   dnresample = (double)nresample;
   x /= dnresample;

   if (mean_resample != NULL) {
      *mean_resample = x;
   }

   rms = (xsqr - nresample * x * x) / (dnresample - 1.0);
   testErrorRetVA(rms < 0, tpcf_npositive, "Resample variance %g cannot be negative", *err, __LINE__, -1.0, rms);

   rms = sqrt(rms);

   if (ERROR == jackknife) {
      rms = rms / sqrt(dnresample) * (dnresample - 1.0);
   }

   return rms;
}

/* ============================================================ *
 * Sets cov_ppmmpm[0,1,2], the resample covariance for xi+,     *
 * xi-, and the mixed covariance, for bins b and c.             *
 * Valid for bootstrap and Jackknife.                           *
 * ============================================================ */
void get_resample_cov_xi(double *cov_ppmmpm, const xi_data *xi, int b, int c, int nresample, error_t ERROR, error **err)
{
   int i, j;
   double pmb[2], pmc[2], pmsqr[3], dnresample;

   for (j=0; j<2; j++) {
      pmb[j] = pmc[j]= pmsqr[j] = 0.0;
   }
   pmsqr[2] = 0.0;

   for (i=0; i<nresample; i++) {
      pmb[0]    += xi->p_resample[b][i];
      pmb[1]    += xi->m_resample[b][i];
      pmc[0]    += xi->p_resample[c][i];
      pmc[1]    += xi->m_resample[c][i];

      pmsqr[0] += xi->p_resample[b][i] * xi->p_resample[c][i];   /* ++ */
      pmsqr[1] += xi->m_resample[b][i] * xi->m_resample[c][i];   /* -- */
      pmsqr[2] += xi->p_resample[b][i] * xi->m_resample[c][i];   /* +- */
   }

   testErrorRetVA(nresample < 0, tpcf_npositive, "Nresample = %d has to be larger than zero", *err, __LINE__,, nresample);
   dnresample = (double)nresample;
   for (j=0; j<2; j++) {
      pmb[j] /= dnresample;
      pmc[j] /= dnresample;
   }

   cov_ppmmpm[0] = (pmsqr[0] - nresample * pmb[0] * pmc[0]) / (dnresample - 1.0);  /* ++ */
   cov_ppmmpm[1] = (pmsqr[1] - nresample * pmb[1] * pmc[1]) / (dnresample - 1.0);  /* -- */
   cov_ppmmpm[2] = (pmsqr[2] - nresample * pmb[0] * pmc[1]) / (dnresample - 1.0);  /* +- */

   if (ERROR == jackknife) {
      for (j=0; j<3; j++) {
         cov_ppmmpm[j] = cov_ppmmpm[j] / dnresample * (dnresample - 1.0) * (dnresample - 1.0);
      }
   }
}

/* ============================================================ *
 * Returns the resample covariance for g_t, for bins b and c.   *
 * Valid for bootstrap and Jackknife.                           *
 * ============================================================ */
double get_resample_cov_gl(const xi_data *xi, int b, int c, int nresample, error_t ERROR, error **err)
{
   double cov, gtb, gtc, gtsqr, dnresample;
   int i;

   gtb = gtc = gtsqr = 0.0;
   for (i=0; i<nresample; i++) {
      gtb   += xi->gt_resample[b][i];
      gtc   += xi->gt_resample[c][i];
      gtsqr += xi->gt_resample[b][i] * xi->gt_resample[c][i];
   }

   testErrorRetVA(nresample < 0, tpcf_npositive, "Nresample = %d has to be larger than zero",
                  *err, __LINE__, 0.0, nresample);
   dnresample = (double)nresample;

   gtb /= dnresample;
   gtc /= dnresample;

   cov = (gtsqr - nresample * gtb * gtc) / (dnresample - 1.0);

   if (ERROR == jackknife) {
      cov = cov / dnresample * (dnresample - 1.0) * (dnresample - 1.0);
   }

   return cov;

}

/* ============================================================ *
 * Returns the bootstrap error of the weighted angular          *
 * correlation function (wn). *Not* used for w(theta), see	    *
 * woftheta_xcorr.pl instead.				                         *
 * ============================================================ */
double get_boots_error(const xi_data *xi, int b, int nboots)
{
   int i;
   double xsqr, x, var;

   if (nboots<2) return 0.0;

   for (i=0,x=xsqr=0.0; i<nboots; i++) {
      x    += xi->npair_resample[b][i];
      xsqr += DSQR(xi->npair_resample[b][i]);
   }
   x  /= (double)nboots;
   var = (xsqr - nboots*x*x)/(double)(nboots-1.0);
   //printf("%d %g %g %d %g\n", b, xsqr, x, nboots, var);

   return sqrt(var);
}

/* Fills header string 'str_header' from list of strings 'header' with length n, *
 * using the correct format with field length 'signif'.                          */
void fill_header(char str_header[], const char *header[], int n, int signif, int nresample, error_t ERROR, coord_t coord, int n_coord)
{
   int b, len;
   char format[32];
   char scoord[32];


   /* First column: account for '#' character */
   if (n_coord == 0) {
      sprintf(format, "%%%ds", signif - 2);
      sprintf(str_header, format, header[0]);
   } else {
      sprintf(scoord, "%s", scoord_t(coord));
      len    = strlen(scoord);
      sprintf(format, "%%%ds[%%s]", signif - 2 - 2 - len);
      sprintf(str_header, format, header[0], scoord);
   }

   /* Subsequent columns */
   for (b=1; b<n; b++) {
      if (n_coord > b) {
         sprintf(format, "%%s%%%ds[%%s]", signif - 2 - len);
         sprintf(str_header, format, str_header, header[b], scoord);
      } else {
         sprintf(format, "%%s%%%ds", signif);
         sprintf(str_header, format, str_header, header[b]);
      }
   }

   if (nresample > 0) {
      sprintf(str_header, "%s  (nresample = %d, type = %s)", str_header, nresample, serror_t(ERROR));
   }
}

/* ============================================================ *
 * Writes the shear-correlation function to file, default "xi". *
 * ============================================================ */
void out_xi(const xi_data *xi, const char *xi_name, const char *xiref_name, const char *xi_name2, const bin_data *bin,
            double sig_eps4, int ncat, coord_t coord_output, error **err)
{
  
   FILE *F, *F2;
   FILE *F3;
   int b;
   double x, sqrtD, sqrtDcor;
   double th1, th2;
   char *header[NCOL_XI] = {COL_NAME_DISTANCE_2D, COL_NAME_XI_P, COL_NAME_XI_M, COL_NAME_XI_X,
      COL_NAME_W_TOT, COL_NAME_SQRT_D, COL_NAME_SQRT_D_COR, COL_NAME_N_PAIR};
   char str_header[128 * NCOL_XI] = "";

   fill_header(str_header, (const char **)header, NCOL_XI, 16, -1, error_none, coord_output, 1);

   F = fileopen(xi_name, "w");
   fprintf(F, "#%s\n", str_header);

   if (strcmp(xiref_name, "") != 0) {
      F2 = fileopen(xiref_name, "w");
      fprintf(F2, "# theta[%s]       g11             g22             g12             g21               w      sqrt(D)   sqrt(Dcor)        n_pair\n",
            scoord_t(coord_output));
   } else {
      F2 = NULL;
   }

   if (strcmp(xi_name2, "") != 0) {
      F3 = fileopen(xi_name2, "w");
      fprintf(F3, "# th_min[%s]   th_max            xi+             xi-             xix               w      sqrt(D)   sqrt(Dcor)        n_pair\n",
            scoord_t(coord_output));
   } else {
      F3 = NULL;
   }			    

   for (b=0; b<bin->N[0]; b++) {

      x   = bin_index_to_scale(bin, xi->wtheta, b);
      th1 = bin_index_to_scale_low(bin, xi->wtheta, b);
      th2 = bin_index_to_scale_high(bin, xi->wtheta, b);

      /* Shape noise terms */
      sqrtD    = calculate_sqrt_D(sig_eps4, xi->wwww[b], xi->ww[b], ncat);
      sqrtDcor = calculate_sqrt_Dcor(sqrtD, xi->npair[b], xi->Nnode[b]);

      rad_to_coordinate(&x, coord_output, err); forwardError(*err, __LINE__,);
      rad_to_coordinate(&th1, coord_output, err); forwardError(*err, __LINE__,);
      rad_to_coordinate(&th2, coord_output, err); forwardError(*err, __LINE__,);

      fprintf(F, "%15.5f % .8e % .8e % .8e % .8e %15.3e %15.3e %15llu\n",
            x, xi->p[b], xi->m[b], xi->x[b], xi->ww[b], sqrtD, sqrtDcor, xi->npair[b]);
      if (F2) {
         fprintf(F2, "%9.5f % .8e % .8e % .8e % .8e % .8e %3e %3e %12llu\n",
               x, xi->g11[b], xi->g22[b], xi->g12[b], xi->g21[b], xi->ww[b], sqrtD, sqrtDcor, xi->npair[b]);
      }
      if (F3) {
         fprintf(F3, "%9.5f %9.5f % .8e % .8e % .8e % .8e %3e %3e %12llu\n",
               th1 , th2 , xi->p[b], xi->m[b], xi->x[b], xi->ww[b], sqrtD, sqrtDcor, xi->npair[b]);
      }
   }
   fileclose(F);
   if (F2) fileclose(F2);
   if (F3) fileclose(F3);
}

/* ============================================================ *
 * Writes the resampling covariance of xi in 3-column format.   *
 * Use cov_ppmmpm_col2block.pl to transform to block format.    *
 * ============================================================ */
void out_cov_xi_resample(const xi_data *xi, const char *cov_name, const bin_data *bin, coord_t coord_output, int nresample,
			 error_t ERROR, error **err)
{
   FILE *F;
   int b, c, j;
   double xb, xc, cov_ppmmpm[3];
   char *header[NCOL_XI_COV_RES] = {COL_NAME_DISTANCE_2D_1, COL_NAME_DISTANCE_2D_2, COL_NAME_COV_PP, COL_NAME_COV_MM, COL_NAME_COV_PM}; 
   char str_header[NCOL_XI_COV_RES * 128] = "";


   fill_header(str_header, (const char **)header, NCOL_XI_RES, 18, nresample, ERROR, coord_output, 2);

   F = fileopen(cov_name, "w");
   fprintf(F, "#%s\n", str_header);

   for (b=0; b<bin->N[0]; b++) {

      xb = bin_index_to_scale(bin, xi->wtheta, b);
      rad_to_coordinate(&xb, coord_output, err); forwardError(*err, __LINE__,);

      for (c=0; c<bin->N[0]; c++) {

         xc = bin_index_to_scale(bin, xi->wtheta, c);
         rad_to_coordinate(&xc, coord_output, err);   forwardError(*err, __LINE__,);

         get_resample_cov_xi(cov_ppmmpm, xi, b, c, nresample, ERROR, err);   forwardError(*err, __LINE__,);

         fprintf(F, "%17.5f %17.5f", xb, xc);
         for (j=0; j<3; j++) fprintf(F, " % 17.5e", cov_ppmmpm[j]);
         fprintf(F, "\n");
      }
   }

   fileclose(F);
}

/* Writes the resampled mean and rms xi to an ascii file */
void out_xi_resample(const xi_data *xi, const char *xi_name, const char *out_ALL_xip_resample, const char *out_ALL_xim_resample,
      const bin_data *bin, coord_t coord_output, int nresample, error_t ERROR, error **err)
{
   FILE *F;
   int b, i;
   double x, mean_pm[2], rms_pm[2];
   char *header[NCOL_XI_RES] = {COL_NAME_DISTANCE_2D, COL_NAME_XI_P_RESAMPLE, COL_NAME_XI_M_RESAMPLE,
                                COL_NAME_RMS_P_RESAMPLE, COL_NAME_RMS_M_RESAMPLE};
   char str_header[128 * NCOL_XI_RES] = "";

   fill_header(str_header, (const char **)header, NCOL_XI_RES, 16, nresample, ERROR, coord_output, 1);

   F = fileopen(xi_name, "w");
   fprintf(F, "#%s\n", str_header);

   /* MKDEBUG TODO: use get_b_index */
   for (b=0; b<bin->N[0]; b++) {

      x = bin_index_to_scale(bin, xi->wtheta, b);   
      rad_to_coordinate(&x, coord_output, err); forwardError(*err, __LINE__,);

      rms_pm[0] = get_resample_error((const double**)xi->p_resample, b, nresample, ERROR, mean_pm, err);
      forwardError(*err, __LINE__,);
      rms_pm[1] = get_resample_error((const double**)xi->m_resample, b, nresample, ERROR, mean_pm+1, err);
      forwardError(*err, __LINE__,);

      fprintf(F, "%15.5f % .8e % .8e %15.3e %15.3e\n", x, mean_pm[0], mean_pm[1], rms_pm[0], rms_pm[1]);
   }

   fileclose(F);
   
   if (strcmp(out_ALL_xip_resample, "") != 0) {
      F = fileopen(out_ALL_xip_resample, "w");
      for (b=0; b<bin->N[0]; b++) {
         x = bin_index_to_scale(bin, xi->wtheta, b);   
         rad_to_coordinate(&x, coord_output, err); forwardError(*err, __LINE__,);
         fprintf(F, "%9.5f", x);
         for (i=0; i<nresample; i++) {
            fprintf(F, " % .8e", xi->p_resample[b][i]);
         }
         fprintf(F, "\n");
      }    
      fileclose(F);
   }

   if (strcmp(out_ALL_xim_resample, "") != 0) { 
      F = fileopen(out_ALL_xim_resample, "w");
      for (b=0; b<bin->N[0]; b++) {
         x = bin_index_to_scale(bin, xi->wtheta, b);   
         rad_to_coordinate(&x, coord_output, err); forwardError(*err, __LINE__,);
         fprintf(F, "%9.5f", x);
           for (i=0; i<nresample; i++) {
	           fprintf(F, " % .8e", xi->m_resample[b][i]);
           }
           fprintf(F, "\n");
       }
       fileclose(F);
   }
}

/* ============================================================ *
 * Writes the angular correlation function to file, default "w" *
 * ============================================================ */
void out_w(const xi_data *xi, const char *w_name, const bin_data *bin, int ngal1, int ngal2,
           coord_t coord_output, int nresample, error_t ERROR, error **err)
{
   FILE *F;
   int b, i;
   double x;
   char **header;
   char *str_header;

   /* Header string */
   str_header = malloc_err(sizeof(char) * 128 * (NCOL_W + nresample), err);
   forwardError(*err, __LINE__,);
   header = malloc_err(sizeof(char*) * (NCOL_W + nresample), err);
   forwardError(*err, __LINE__,);
   for (i=0; i<NCOL_W+nresample; i++) {
      header[i] = malloc_err(sizeof(char) * 128, err);
      forwardError(*err, __LINE__,);
   }

   sprintf(header[0], "%s", COL_NAME_DISTANCE_2D);
   sprintf(header[1], " %s", COL_NAME_N_PAIR);
   for (i=0; i<nresample; i++) {
      sprintf(header[i+2], " %s_%d", COL_NAME_N_PAIR_RESAMPLE, i);
   }
   fill_header(str_header, (const char **)header, NCOL_W + nresample, 14, nresample, ERROR, coord_output, 1);

   F = fileopen(w_name, w_name);
   fprintf(F, "# Number of galaxies for both catalogues = (%u, %u)\n", ngal1, ngal2);
   fprintf(F, "#%s\n", str_header); 

   for (b=0; b<bin->N[0]; b++) {

      x = bin_index_to_scale(bin, xi->wtheta, b);
      rad_to_coordinate(&x, coord_output, err); forwardError(*err, __LINE__,);

      fprintf(F, "%13.5g %13llu", x, xi->npair[b]);
      for (i=0; i<nresample; i++) {
          fprintf(F, " %13.5g", xi->npair_resample[b][i]);
      }
      fprintf(F, "\n");
   }
   fileclose(F);

   for (i=0; i<NCOL_W+nresample; i++) {
      free(header[i]);
   }
   free(header);
   free(str_header);

}

/* ============================================================ *
 * Writes the 2d (radial+line-of-sight) correlation function to *
 * file, default "w".					                            *
 * ============================================================ */
void out_wp_rp_pi(const xi_data *xi, const char *w_name, const bin_data *bin, int nresample, coord_t coord_output,
	   uint ngal1, uint ngal2, error **err)
{
   FILE *F;
   int b[2], b_index, i;
   double x[2];

   F = fileopen(w_name, w_name);
   fprintf(F, "# Number of galaxies for both catalogues = (%u, %u)\n", ngal1, ngal2);
   // TODO: units. Mpc?
   fprintf(F, "# rp[%s] pi[%s] n_pair n_pair_resample (%d realizations)\n", scoord_t(coord_output), scoord_t(coord_output), nresample);

   for (b[0]=0; b[0]<bin->N[0]; b[0]++) {
      for (b[1]=0; b[1]<bin->N[1]; b[1]++) {

         b_index = get_bin_index(b, bin->N, bin->ndim);

         for (i=0; i<2; i++) {
            x[i] = bin_index_to_scale(bin, xi->wtheta, b[i]);
            rad_to_coordinate(x+i, coord_output, err); forwardError(*err, __LINE__,); // TODO units
         }

         fprintf(F, "%9.5g %9.5g %15llu", x[0], x[1], xi->npair[b_index]);
         for (i=0; i<nresample; i++) {
            fprintf(F, " %15.5g", xi->npair_resample[b_index][i]);
         }
         //double err_bt; err_bt = get_boots_error(xi, b, nresample); fprintf(F, " %g", err_bt);
         fprintf(F, "\n");
      }
   }

   fileclose(F);
}


/* ============================================================== *
 * Writes the weighted angular correlation function to file,      *
 * default "ww".						                                 *
 * ============================================================== */
void out_ww(const xi_data *xi, const char *w_name, const bin_data *bin, int nboots, coord_t coord_output, double wbar, error **err)
{
   FILE *F;
   int b, i;
   double x, err_bt, sigma;

   F = fileopen(w_name, "w");
   fprintf(F, "# theta[%s] weight-wbar sigma(weight)\n", scoord_t(coord_output));
   fprintf(F, "# wbar = %.8g\n", wbar);

   for (b=0; b<bin->N[0]; b++) {

      x = bin_index_to_scale(bin, xi->wtheta, b);

      /* Width of distribution */
      // MKDEBUG: correct weights used?
      sigma  = 1.0/((double)xi->npair[b]-1.0)*(xi->wnwn[b] - (double)xi->npair[b]*DSQR(xi->wn[b]));
      sigma /= (double)xi->npair[b];  /* Error of mean */
      sigma  = sqrt(sigma);

      err_bt = get_boots_error(xi, b, nboots);
      // MKDEBUG: wbar not subtracted

      rad_to_coordinate(&x, coord_output, err); forwardError(*err, __LINE__,);
      fprintf(F, "%9.5g %14.10g %14.10g %15llu %14.10g", x, xi->wn[b], sigma, xi->npair[b], err_bt);
      for (i=0; i<nboots; i++) {
         //fprintf(F, " %15.5g", xi->npair_resample[b][i]);
      }
      fprintf(F, "\n");
   }
   fileclose(F);
}


/* ============================================================ *
 * Writes the galaxy-shear cross-correlation function to file,  *
 * default "wgl".						                               *
 * ============================================================ */
void out_gl(const xi_data *xi, const char *gl_name, const bin_data *bin, double sig_eps2, coord_t coord_output, error **err)
{
   FILE *F;
   int b;
   double x, sqrtD, sqrtDcor;
   char *header[NCOL_WGL] = {COL_NAME_DISTANCE_2D, COL_NAME_GT, COL_NAME_GX, COL_NAME_W_TOT, COL_NAME_SQRT_D,
                             COL_NAME_SQRT_D_COR, COL_NAME_N_PAIR};
   char str_header[128 * NCOL_WGL];

   fill_header(str_header, (const char **)header, NCOL_WGL, 16, -1, error_none, coord_output, 1);
  
   F = fileopen(gl_name, "w");
   fprintf(F, "#%s\n", str_header);

   for (b=0; b<bin->N[0]; b++) {
      x = bin_index_to_scale(bin, xi->wtheta, b);
      rad_to_coordinate(&x, coord_output, err); forwardError(*err, __LINE__,);

      sqrtD    = calculate_sqrt_D(sig_eps2, xi->wnwn[b], xi->wn[b], 2);
      sqrtDcor = calculate_sqrt_Dcor(sqrtD, xi->npair[b], xi->Nnode[b]);

      fprintf(F, "%15.5f % .8e % .8e % .8e % .8e % .8e %15llu \n", x, xi->gt[b], xi->gx[b], xi->wn[b], sqrtD, sqrtDcor, xi->npair[b]);
   }

   fileclose(F);
}

void out_gl_resample(const xi_data *xi, const char *gl_name, const bin_data *bin,const char *out_ALL_gt_resample ,coord_t coord_output,
                      int nresample, error_t ERROR, error **err)
{
   FILE *F;
   int b,i;
   double x, rms, mean;
   char *header[NCOL_WGL_RES] = {COL_NAME_DISTANCE_2D, COL_NAME_GT_RESAMPLE, COL_NAME_GT_RMS_RESAMPLE};
   char str_header[128 * NCOL_WGL_RES];

   fill_header(str_header, (const char **)header, NCOL_WGL_RES, 18, nresample, ERROR, coord_output, 1);

   F = fileopen(gl_name, "w");
   fprintf(F, "#%s\n", str_header);

   for (b=0; b<bin->N[0]; b++) {

      x = bin_index_to_scale(bin, xi->wtheta, b);
      rad_to_coordinate(&x, coord_output, err); forwardError(*err, __LINE__,);

      rms = get_resample_error((const double**)xi->gt_resample, b, nresample, ERROR, &mean, err);
      forwardError(*err, __LINE__,);

      fprintf(F, "%17.5f % 17.8e % 17.8e\n", x, mean, rms);

   }

   fileclose(F);

   if (strcmp(out_ALL_gt_resample, "") != 0) {
       F = fileopen(out_ALL_gt_resample, "w");
       for (b=0; b<bin->N[0]; b++) {
           x = bin_index_to_scale(bin, xi->wtheta, b);   
           rad_to_coordinate(&x, coord_output, err); forwardError(*err, __LINE__,);
           fprintf(F, "%9.5f", x);
           for (i=0; i<nresample; i++) {
	           fprintf(F, " % .8e", xi->gt_resample[b][i]);
           }
           fprintf(F, "\n");
       }    
   fileclose(F);
   }

}

/* ============================================================ *
 * Writes the resampling covariance of g_t in column format.    *
 * ============================================================ */
void out_cov_gl_resample(const xi_data *xi, const char *cov_name, const bin_data *bin, coord_t coord_output, int nresample,
          error_t ERROR, error **err)
{
   FILE *F;
   int b, c;
   double xb, xc, cov;
   char *header[NCOL_COV_COL] = {COL_NAME_DISTANCE_2D_1, COL_NAME_DISTANCE_2D_2, COL_NAME_GT_COV};
   char str_header[NCOL_COV_COL * 128] = "";

   fill_header(str_header, (const char **)header, NCOL_COV_COL, 18, nresample, ERROR, coord_output, 2);
   
   F = fileopen(cov_name, "w");
   fprintf(F, "#%s\n", str_header);

   for (b=0; b<bin->N[0]; b++) {

      xb = bin_index_to_scale(bin, xi->wtheta, b);
      rad_to_coordinate(&xb, coord_output, err); forwardError(*err, __LINE__,);

      for (c=0; c<bin->N[0]; c++) {

         xc = bin_index_to_scale(bin, xi->wtheta, c);
         rad_to_coordinate(&xc, coord_output, err);   forwardError(*err, __LINE__,);

         cov = get_resample_cov_gl(xi, b, c, nresample, ERROR, err);
         forwardError(*err, __LINE__,);

         fprintf(F, "%17.5f %17.5f %17.5e\n", xb, xc, cov);
      }
   }

   fileclose(F);
}

/* ============================================================ *
 * Returns 1 and sets error if input format is not in           *
 * with correlation type wcorr, i.e., if necessary input        *
 * columns are not given.                                       * 
 * ============================================================ */
int check_input_format(wcorr_t wcorr, format_t format, const char **col_names, int Ncol, error **err)
{
   int result;

   if (format == f_fits) {

#ifdef _WITH_FITS
      result = check_input_format_fits(wcorr, col_names, Ncol, err);
      forwardError(*err, __LINE__, 1);
      return result;
#else
      *err = addError(gc_fits, "Re-compile 'athena' with fits support", *err, __LINE__);
      return 1;
#endif

   } else {

      result = check_input_format_format(wcorr, format, err);
      forwardError(*err, __LINE__, 1);
      return result;

   }

}

/* Checks input catalogue format for explicit format (not fits) */
int check_input_format_format(wcorr_t wcorr, format_t format, error **err)
{
   if (wcorr & gg || wcorr & gn) {
      testErrorRet(format == f_position || format == f_position_jack_num, gc_format,
                   "Format 'position' or 'position_jack_num' does not imply ellipticities, necessary for lensing",
                   *err, __LINE__, 1);
   }

   return 0;
}

#ifdef _WITH_FITS
/* Writes fits file with xi, and, if ERROR != none, xi.resample, an xi.resample.cov. *
 * part is written as separate hdu binary table.                                     */
void out_xi_fits(const xi_data *xi, const char *xi_name, const bin_data *bin,
                 double sig_eps4, int ncat, coord_t coord_output, int nresample, error_t ERROR, int quiet, error **err)
{
   fitsfile *fptr;
   int status = 0, b, c, i, j;
   LONGLONG naxis2 = 0;
   char *ttype[NCOL_XI] = {COL_NAME_DISTANCE_2D, COL_NAME_XI_P, COL_NAME_XI_M, COL_NAME_XI_X,
                           COL_NAME_W_TOT, COL_NAME_SQRT_D, COL_NAME_SQRT_D_COR, COL_NAME_N_PAIR};
   char *tform[NCOL_XI] = {"D",  "D", "D", "D",  "D",  "D", "D",  "K"};
   int data_type[NCOL_XI] = {TDOUBLE, TDOUBLE, TDOUBLE, TDOUBLE, TDOUBLE, TDOUBLE, TDOUBLE, TULONG};
   void *ptr_data[NCOL_XI];
   char excl_xi_name[128];
   double *theta, *sqrtD, *sqrtDcor, *mean_p_res, *mean_m_res, *rms_p_res, *rms_m_res, *theta2, *cov[3], cov_ppmmpm[3];


   /* xi */

   /* Create fits table */
   sprintf(excl_xi_name, "!%s", xi_name);
   fits_create_file(&fptr, excl_xi_name, &status);
   FITS_STATUS_ERROR(status, quiet, "Error (status=%d) while creating fits file '%s'", *err, __LINE__,, status, xi_name);

   fits_create_tbl(fptr, BINARY_TBL, naxis2, NCOL_XI, ttype, tform, NULL, NULL, &status);
   FITS_STATUS_ERROR(status, quiet, "Error (status=%d) while creating fits table", *err, __LINE__,, status);

   write_header_keys(fptr, SHEAR_CORR_NAME, "Shear two-point correlation function", coord_output, -1, error_none, quiet, err);
   forwardError(*err, __LINE__,);
   FITS_STATUS_ERROR(status, quiet, "Error (status=%d) while writing OBJECT to fits header", *err, __LINE__,, status);


   theta    = malloc_err(sizeof(double) * bin->N[0], err); forwardError(*err, __LINE__,);
   sqrtD    = malloc_err(sizeof(double) * bin->N[0], err); forwardError(*err, __LINE__,);
   sqrtDcor = malloc_err(sizeof(double) * bin->N[0], err); forwardError(*err, __LINE__,);

   for (b=0; b<bin->N[0]; b++) {
      theta[b] = bin_index_to_scale(bin, xi->wtheta, b);
      rad_to_coordinate(theta + b, coord_output, err);
      forwardError(*err, __LINE__,);

      sqrtD[b]    = calculate_sqrt_D(sig_eps4, xi->wwww[b], xi->ww[b], ncat);
      sqrtDcor[b] = calculate_sqrt_Dcor(sqrtD[b], xi->npair[b], xi->Nnode[b]); 
   }

   ptr_data[0] = theta;
   ptr_data[1] = xi->p;
   ptr_data[2] = xi->m;
   ptr_data[3] = xi->x;
   ptr_data[4] = xi->ww;
   ptr_data[5] = sqrtD;
   ptr_data[6] = sqrtDcor;
   ptr_data[7] = xi->npair;


   write_columns(fptr, data_type, (void **)ptr_data, (const char **)ttype, bin->N[0], NCOL_XI, quiet, err);
   forwardError(*err, __LINE__,);

   free(theta); free(sqrtD); free(sqrtDcor);


   if (ERROR != error_none) {

      /* Resampling mean and rms */

      ttype[1] = COL_NAME_XI_P_RESAMPLE;
      ttype[2] = COL_NAME_XI_M_RESAMPLE;
      ttype[3] = COL_NAME_RMS_P_RESAMPLE;
      ttype[4] = COL_NAME_RMS_M_RESAMPLE;

      fits_create_tbl(fptr, BINARY_TBL, naxis2, NCOL_XI_RES, ttype, tform, NULL, NULL, &status);
      FITS_STATUS_ERROR(status, quiet, "Error (status=%d) while creating fits table", *err, __LINE__,, status);

      write_header_keys(fptr, SHEAR_RESAMPLE_CORR_NAME, "Resampled shear two-point correlation function mean and rms",
                        coord_output, nresample, ERROR, quiet, err);
      forwardError(*err, __LINE__,);

      mean_p_res = malloc_err(sizeof(double) * bin->N[0], err); forwardError(*err, __LINE__,);
      mean_m_res = malloc_err(sizeof(double) * bin->N[0], err); forwardError(*err, __LINE__,);
      rms_p_res  = malloc_err(sizeof(double) * bin->N[0], err); forwardError(*err, __LINE__,);
      rms_m_res  = malloc_err(sizeof(double) * bin->N[0], err); forwardError(*err, __LINE__,);

      for (b=0; b<bin->N[0]; b++) {
         rms_p_res[b] = get_resample_error((const double**)xi->p_resample, b, nresample, ERROR, mean_p_res+b, err);
         rms_m_res[b] = get_resample_error((const double**)xi->m_resample, b, nresample, ERROR, mean_m_res+b, err);
      }

      ptr_data[1] = mean_p_res;
      ptr_data[2] = mean_m_res;
      ptr_data[3] = rms_p_res;
      ptr_data[4] = rms_m_res;

      write_columns(fptr, data_type, (void**)ptr_data, (const char **)ttype, bin->N[0], NCOL_XI_RES, quiet, err);
      forwardError(*err, __LINE__,);

      free(mean_p_res);
      free(mean_m_res);
      free(rms_p_res);
      free(rms_m_res);


      /* Resampling covariance */
      
      ttype[0] = COL_NAME_DISTANCE_2D_1;
      ttype[1] = COL_NAME_DISTANCE_2D_2;
      ttype[2] = COL_NAME_COV_PP;
      ttype[3] = COL_NAME_COV_MM;
      ttype[4] = COL_NAME_COV_PM;

      fits_create_tbl(fptr, BINARY_TBL, naxis2, NCOL_XI_COV_RES, ttype, tform, NULL, NULL, &status);
      FITS_STATUS_ERROR(status, quiet, "Error (status=%d) while creating fits table", *err, __LINE__,, status);

      write_header_keys(fptr, SHEAR_RESAMPLE_COV_NAME, "Resampled shear two-point correlation function covariance",
                         coord_output, nresample, ERROR, quiet, err);
      forwardError(*err, __LINE__,);

      theta  = malloc_err(sizeof(double) * bin->N[0] * bin->N[0], err); forwardError(*err, __LINE__,);
      theta2 = malloc_err(sizeof(double) * bin->N[0] * bin->N[0], err); forwardError(*err, __LINE__,);
      for (j=0; j<3; j++) {
         cov[j] = malloc_err(sizeof(double) * bin->N[0] * bin->N[0], err);
         forwardError(*err, __LINE__,);
      }

      for (b=0,i=0; b<bin->N[0]; b++) {
         for (c=0; c<bin->N[0]; c++,i++) {
            theta[i] = bin_index_to_scale(bin, xi->wtheta, b);
            rad_to_coordinate(theta+i, coord_output, err); forwardError(*err, __LINE__,);
            theta2[i] = bin_index_to_scale(bin, xi->wtheta, c);
            rad_to_coordinate(theta2+i, coord_output, err); forwardError(*err, __LINE__,);

            get_resample_cov_xi(cov_ppmmpm, xi, b, c, nresample, ERROR, err);   forwardError(*err, __LINE__,);
            for (j=0; j<3; j++) {
               cov[j][i] = cov_ppmmpm[j];
            }
         }
      }

      ptr_data[0] = theta;
      ptr_data[1] = theta2;
      ptr_data[2] = cov[0];
      ptr_data[3] = cov[1];
      ptr_data[4] = cov[2];

      write_columns(fptr, data_type, (void **)ptr_data, (const char **)ttype, bin->N[0]*bin->N[0], NCOL_XI_COV_RES, quiet, err);
      forwardError(*err, __LINE__,);

      free(theta); free(theta2);
      for (j=0; j<3; j++) free(cov[j]);

   }

   fits_close_file(fptr, &status);
   FITS_STATUS_ERROR(status, quiet, "Error (status=%d) while closing fits file", *err, __LINE__,, status);
}

void out_gl_fits(const xi_data *xi, const char *gl_name, const bin_data *bin, double sig_eps2,
      coord_t coord_output, int nresample, error_t ERROR, int quiet, error **err)
{
   fitsfile *fptr;
   int status = 0, b, c, i;
   LONGLONG naxis2 = 0;
   char *ttype[NCOL_WGL] = {COL_NAME_DISTANCE_2D, COL_NAME_GT, COL_NAME_GX, COL_NAME_W_TOT, COL_NAME_SQRT_D, COL_NAME_SQRT_D_COR,
      COL_NAME_N_PAIR};
   char *tform[NCOL_WGL] = {"D",  "D", "D",  "D",  "D", "D",  "K"};
   int data_type[NCOL_WGL] = {TDOUBLE, TDOUBLE, TDOUBLE, TDOUBLE, TDOUBLE, TDOUBLE, TULONG};
   void *ptr_data[NCOL_WGL];
   char excl_gl_name[128];
   double *theta, *sqrtD, *sqrtDcor, *mean_res, *rms_res, *theta2, *cov;


   /* Create fits table */
   sprintf(excl_gl_name, "!%s", gl_name);
   fits_create_file(&fptr, excl_gl_name, &status);
   FITS_STATUS_ERROR(status, quiet, "Error (status=%d) while creating fits file '%s'", *err, __LINE__,, status, gl_name);

   fits_create_tbl(fptr, BINARY_TBL, naxis2, NCOL_WGL, ttype, tform, NULL, NULL, &status);
   FITS_STATUS_ERROR(status, quiet, "Error (status=%d) while creating fits table", *err, __LINE__,, status);

   write_header_keys(fptr, GAL_SHEAR_XCORR_NAME, "Galaxy-shear cross-correlation function", coord_output, -1, error_none, quiet, err);
   forwardError(*err, __LINE__,);

   theta    = malloc_err(sizeof(double) * bin->N[0], err); forwardError(*err, __LINE__,);
   sqrtD    = malloc_err(sizeof(double) * bin->N[0], err); forwardError(*err, __LINE__,);
   sqrtDcor = malloc_err(sizeof(double) * bin->N[0], err); forwardError(*err, __LINE__,);

   for (b=0; b<bin->N[0]; b++) {
      theta[b] = bin_index_to_scale(bin, xi->wtheta, b);
      rad_to_coordinate(theta + b, coord_output, err);
      forwardError(*err, __LINE__,);

      sqrtD[b]    = calculate_sqrt_D(sig_eps2, xi->wnwn[b], xi->wn[b], 2);
      sqrtDcor[b] = calculate_sqrt_Dcor(sqrtD[b], xi->npair[b], xi->Nnode[b]);
   }

   ptr_data[0] = theta;
   ptr_data[1] = xi->gt;
   ptr_data[2] = xi->gx;
   ptr_data[3] = xi->wn;
   ptr_data[4] = sqrtD;
   ptr_data[5] = sqrtDcor;
   ptr_data[6] = xi->npair;

   write_columns(fptr, data_type, (void **)ptr_data, (const char **)ttype, bin->N[0], NCOL_WGL, quiet, err);
   forwardError(*err, __LINE__,);

   free(sqrtD); free(sqrtDcor);


   if (ERROR != error_none) {

      /* Resampling mean and rms */

      ttype[1] = COL_NAME_GT_RESAMPLE;
      ttype[2] = COL_NAME_GT_RMS_RESAMPLE;

      fits_create_tbl(fptr, BINARY_TBL, naxis2, NCOL_WGL_RES, ttype, tform, NULL, NULL, &status);
      FITS_STATUS_ERROR(status, quiet, "Error (status=%d) while creating fits table", *err, __LINE__,, status);

      write_header_keys(fptr, GAL_SHEAR_RESAMPLE_XCORR_NAME, "Resampled gal-shear xcorr fct mean and rms",
                 coord_output, nresample, ERROR, quiet, err);
      forwardError(*err, __LINE__,);

      mean_res = malloc_err(sizeof(double) * bin->N[0], err); forwardError(*err, __LINE__,);
      rms_res  = malloc_err(sizeof(double) * bin->N[0], err); forwardError(*err, __LINE__,);

      for (b=0; b<bin->N[0]; b++) {
         rms_res[b] = get_resample_error((const double**)xi->gt_resample, b, nresample, ERROR, mean_res+b, err);
      }

      ptr_data[1] = mean_res;
      ptr_data[2] = rms_res;

      write_columns(fptr, data_type, (void **)ptr_data, (const char **)ttype, bin->N[0], NCOL_WGL_RES, quiet, err);
      forwardError(*err, __LINE__,);

      free(mean_res);
      free(rms_res);


      /* Resampling covariance */

      ttype[0] = COL_NAME_DISTANCE_2D_1;
      ttype[1] = COL_NAME_DISTANCE_2D_2;
      ttype[2] = COL_NAME_GT_COV;

      fits_create_tbl(fptr, BINARY_TBL, naxis2, NCOL_COV_COL, ttype, tform, NULL, NULL, &status);
      FITS_STATUS_ERROR(status, quiet, "Error (status=%d) while creating fits table", *err, __LINE__,, status);
      forwardError(*err, __LINE__,);

      write_header_keys(fptr, GAL_SHEAR_COV_XCORR_NAME, "Resampled gal-shear xcorr fct covariance",
                        coord_output, nresample, ERROR, quiet, err);
      forwardError(*err, __LINE__,);

      free(theta); /* Used in wgl and wgl.resample */

      theta  = malloc_err(sizeof(double) * bin->N[0] * bin->N[0], err); forwardError(*err, __LINE__,);
      theta2 = malloc_err(sizeof(double) * bin->N[0] * bin->N[0], err); forwardError(*err, __LINE__,);
      cov    = malloc_err(sizeof(double) * bin->N[0] * bin->N[0], err); forwardError(*err, __LINE__,);

      for (b=0,i=0; b<bin->N[0]; b++) {
         for (c=0; c<bin->N[0]; c++,i++) {
            theta[i] = bin_index_to_scale(bin, xi->wtheta, b);
            rad_to_coordinate(theta+i, coord_output, err); forwardError(*err, __LINE__,);
            theta2[i] = bin_index_to_scale(bin, xi->wtheta, c);
            rad_to_coordinate(theta2+i, coord_output, err); forwardError(*err, __LINE__,);

            cov[i] = get_resample_cov_gl(xi, b, c, nresample, ERROR, err); 
         }
      }

      ptr_data[0] = theta;
      ptr_data[1] = theta2;
      ptr_data[2] = cov;

      write_columns(fptr, data_type, (void **)ptr_data, (const char **)ttype, bin->N[0]*bin->N[0], NCOL_COV_COL, quiet, err);
      forwardError(*err, __LINE__,);

      free(theta);
      free(theta2);
      free(cov);

   }

   fits_close_file(fptr, &status);
   FITS_STATUS_ERROR(status, quiet, "Error (status=%d) while closing fits file", *err, __LINE__,, status);
}

/* Writes fits file with the 2D projected clustering correlation function w. *
 * Resampled number of pairs are written as extra columns in the same hdu.   */
void out_w_fits(const xi_data *xi, const char *w_name, const bin_data *bin,
                 int ngal1, int ngal2, coord_t coord_output, int nresample,
                 error_t ERROR, int quiet, error **err)
{  
   fitsfile *fptr;
   int status = 0, b, i;
   LONGLONG naxis2 = 0;
   char **ttype, **tform;
   int *data_type;
   void **ptr_data;
   char excl_w_name[128];
   double *theta, **npair_resample;


   data_type = malloc_err(sizeof(int) * (NCOL_W + nresample), err);
   forwardError(*err, __LINE__,);
   ttype    = malloc_err(sizeof(char*) * (NCOL_W + nresample), err);
   forwardError(*err, __LINE__,);
   tform    = malloc_err(sizeof(char*) * (NCOL_W + nresample), err);
   forwardError(*err, __LINE__,);
   for (i=0; i<NCOL_W+nresample; i++) {
      ttype[i] = malloc_err(sizeof(char) * 128, err);
      forwardError(*err, __LINE__,);
      tform[i] = malloc_err(sizeof(char) * 8, err);
      forwardError(*err, __LINE__,);
   }

   data_type[0] = TDOUBLE;
   data_type[1] = TULONG;
   sprintf(tform[0], "%s", "D");
   sprintf(tform[1], "%s", "K");
   sprintf(ttype[0], "%s", COL_NAME_DISTANCE_2D);
   sprintf(ttype[1], "%s", COL_NAME_N_PAIR);
   for (i=0; i<nresample; i++) {
      data_type[i+2] = TDOUBLE;
      sprintf(tform[i+2], "%s", "D");
      sprintf(ttype[i+2], "%s_%d", COL_NAME_N_PAIR_RESAMPLE, i); 
   }


   /* Create fits table */
   sprintf(excl_w_name, "!%s", w_name);
   fits_create_file(&fptr, excl_w_name, &status);
   FITS_STATUS_ERROR(status, quiet, "Error (status=%d) while creating fits file '%s'", *err, __LINE__,, status, w_name);

   fits_create_tbl(fptr, BINARY_TBL, naxis2, NCOL_W + nresample, ttype, tform, NULL, NULL, &status);
   FITS_STATUS_ERROR(status, quiet, "Error (status=%d) while creating fits table", *err, __LINE__,, status);

   write_header_keys(fptr, CORR_NAME, "Spatial correlation function", coord_output, nresample, ERROR, quiet, err);
   forwardError(*err, __LINE__,);

   fits_write_key(fptr, TINT, "NGAL1", &ngal1, "Number of galaxies in catalogue 1", &status);
   FITS_STATUS_ERROR(status, quiet, "Error (status=%d) while writing NGAL1 to fits header", *err, __LINE__,, status);
   fits_write_key(fptr, TINT, "NGAL2", &ngal2, "Number of galaxies in catalogue 2", &status);
   FITS_STATUS_ERROR(status, quiet, "Error (status=%d) while writing NGAL2 to fits header", *err, __LINE__,, status);

   theta          = malloc_err(sizeof(double) * bin->N[0], err);
   forwardError(*err, __LINE__,);
   npair_resample = calloc_array_double(nresample, bin->N[0], err);
   forwardError(*err, __LINE__,);
  
   for (b=0; b<bin->N[0]; b++) {
      theta[b] = bin_index_to_scale(bin, xi->wtheta, b);
      rad_to_coordinate(theta + b, coord_output, err);
      forwardError(*err, __LINE__,);
      for (i=0; i<nresample; i++) {
         npair_resample[i][b] = xi->npair_resample[b][i];
      }
   }

   ptr_data = malloc_err(sizeof(void*) * (NCOL_W + nresample), err);
   forwardError(*err, __LINE__,);

   ptr_data[0] = theta;
   ptr_data[1] = xi->npair;
   for (i=0; i<nresample; i++) {
      ptr_data[i+2] = npair_resample[i];
   }

   write_columns(fptr, data_type, (void **)ptr_data, (const char **)ttype, bin->N[0], NCOL_W + nresample, quiet, err);
   forwardError(*err, __LINE__,);

   for (i=0; i<NCOL_W+nresample; i++) {
      free(ttype[i]);
      free(tform[i]);
   }
   for (i=0; i<nresample; i++) {
      free(npair_resample[i]);
   }
   free(npair_resample);
   free(ttype);
   free(tform);
   free(ptr_data);
   free(data_type);
   free(theta);

   fits_close_file(fptr, &status);
   FITS_STATUS_ERROR(status, quiet, "Error (status=%d) while closing fits file", *err, __LINE__,, status);
}

/* Checks input catalogue format for implicit format (fits) */
int check_input_format_fits(wcorr_t wcorr, const char **col_names, int Ncol, error **err)
{
   int j;
   col_descr col;
   int n_pos, n_ell;

   n_pos = n_ell = 0;
   for (j=0; j<Ncol; j++) {
      get_col(col_names[j], &col, err);
      forwardError(*err, __LINE__, 1);

      switch (col.key) {
         case c_x : case c_y :
            n_pos ++;
            break;
         case c_e1 : case c_e2 :
            n_ell ++;
            break;
         default :
            break;
      }
   }

   testErrorRet(n_pos < 2, gc_format, "Less than two spatial coordinate columns in input catalogue",
         *err, __LINE__, 1);

   if ((wcorr & gg) || (wcorr & gn)) {
      testErrorRet(n_ell < 2, gc_format, "Less than two ellipticity columns in input catalogue",
            *err, __LINE__, 1);
   }

   return 0;
}


void write_header_keys(fitsfile *fptr, char *object, const char *object_comment,
                       coord_t coord, int nresample, error_t ERROR, int quiet, error **err)
{
   int status = 0;

   fits_write_key(fptr, TSTRING, "OBJECT", object, object_comment, &status);
   FITS_STATUS_ERROR(status, quiet, "Error (status=%d) while writing OBJECT to fits header", *err, __LINE__,, status);
   fits_write_key(fptr, TSTRING, "UNITS", scoord_t(coord), "Coordinate units", &status);
   FITS_STATUS_ERROR(status, quiet, "Error (status=%d) while writing OBJECT to fits header", *err, __LINE__,, status);

   if (ERROR != error_none) {
      fits_write_key(fptr, TSTRING, "TYPE", serror_t(ERROR), "Resample type", &status);
      FITS_STATUS_ERROR(status, quiet, "Error (status=%d) while writing TYPE to fits header", *err, __LINE__,, status);
      fits_write_key(fptr, TINT, "NRESAMPLE", &nresample, "Number of resamples", &status);
      FITS_STATUS_ERROR(status, quiet, "Error (status=%d) while writing NRESAMPLE to fits header", *err, __LINE__,, status);
   }
}

void write_columns(fitsfile *fptr, const int *data_type, void *ptr_data[], const char *ttype[],
                    int Nbin, int Ncol, int quiet, error **err)
{
   int colnum, status = 0;
   long firstrow=1, firstelem=1;

   for (colnum=0; colnum<Ncol; colnum++) {
      fits_write_col(fptr, data_type[colnum], colnum+1, firstrow, firstelem, Nbin, ptr_data[colnum], &status);
      FITS_STATUS_ERROR(status, quiet, "Error (status=%d) while writing column '%s' to fits file",
                        *err, __LINE__,, status, ttype[colnum]);
   }
   
}

#endif


void free_xi_data(xi_data *xi, uint N)
{
   int b;

   free(xi->npair);
   if (xi->p)     free_vector(xi->p , 0, N - 1);
   if (xi->m)     free_vector(xi->m , 0, N - 1);
   if (xi->x)     free_vector(xi->x , 0, N - 1);
   if (xi->g11)   free_vector(xi->g11 , 0, N - 1);
   if (xi->g12)   free_vector(xi->g12 , 0, N - 1);
   if (xi->g21)   free_vector(xi->g21 , 0, N - 1);
   if (xi->g22)   free_vector(xi->g22 , 0, N - 1);
   if (xi->ww)    free_vector(xi->ww , 0, N - 1);
   if (xi->wn)    free_vector(xi->wn, 0, N - 1);
   if (xi->wwww)  free_vector(xi->wwww, 0, N - 1);
   if (xi->wnwn)  free_vector(xi->wnwn, 0, N - 1);
   if (xi->Nnode) free_ivector(xi->Nnode, 0, N - 1);
   if (xi->gt)    free_vector(xi->gt, 0, N - 1);
   if (xi->gx)    free_vector(xi->gx, 0, N - 1);
   if (xi->npair_resample) {
      for (b=0; b<N; b++) free(xi->npair_resample[b]);
      free(xi->npair_resample);
   }
   if (xi->p_resample) {
      for (b=0; b<N; b++) free(xi->p_resample[b]);
      free(xi->p_resample);
   }
   if (xi->m_resample) {
      for (b=0; b<N; b++) free(xi->m_resample[b]);
      free(xi->m_resample);
   }
   if (xi->ww_resample) {
      for (b=0; b<N; b++) free(xi->ww_resample[b]);
      free(xi->ww_resample);
   }
   if (xi->wn_resample) {
      for (b=0; b<N; b++) free(xi->wn_resample[b]);
      free(xi->wn_resample);
   }
   if (xi->wtheta) free(xi->wtheta);
}
