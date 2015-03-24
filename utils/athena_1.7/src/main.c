/* ============================================================ *
 * main.c							                                  *
 * Tree-code 'athena'					                       	    *
 * Martin Kilbinger, Christopher Bonnett                        *
 * 2007-    						                                  *
 * ============================================================ */


#include "main.h"


unsigned long long npairtot;
long double p=0, f=0; /* global variable for progress indicator */
int quiet_gl;


/* ============================================================ *
 * Reads the configuration file with name `cname' and copies    *
 * the information to the structure `config'.			          *
 * ============================================================ */
void read_config_file(config_info *config, char cname[], int *ndim, int *ndimbin, error **err)
{
   FILE *F;
   config_element c = {0, 0.0, ""};
   int i, j, k;
   char s[512];

   F = fileopen(cname, "r");

   CONFIG_READ_S(config, GALCAT1, s, F, c, err); /* First catalogue  */
   CONFIG_READ_S(config, GALCAT2, s, F, c, err); /* Second catalogue */

   /* Check if there are one or two catalogues to be considered */
   if (strcmp(config->GALCAT2, "-")==0) strcpy(config->GALCAT2, config->GALCAT1);
   if (strcmp(config->GALCAT1, config->GALCAT2)==0) config->NCAT = 1;
   else config->NCAT = 2;

   CONFIG_READ(config, WCORR, i, F, c, err);

   if ((config->WCORR & nn) || (config->WCORR & wn)) {

      CONFIG_READ_S(config, SWCORR_SUBTYPE, s, F, c, err);
      STRING2ENUM(config->WCORR_SUBTYPE, config->SWCORR_SUBTYPE, wcorr_subtype_t, swcorr_subtype_t,
            j, Nwcorr_subtype, err);

      testErrorRetVA(config->WCORR == nn && (config->WCORR_SUBTYPE != nn_2d && config->WCORR_SUBTYPE != nn_3d &&
               config->WCORR_SUBTYPE != nn_rp_pi),
            tree_wcorr, "WCORR == %d (nn) not compatible with WCORR_SUBTYPE = %d (%s)", *err, __LINE,,
            (int)config->WCORR, (int)config->WCORR_SUBTYPE, swcorr_subtype_t(config->WCORR_SUBTYPE));

      if (config->WCORR_SUBTYPE==wn_pairwise) {
         CONFIG_READ(config, WN_ZGAP_MIN, d, F, c, err);
         CONFIG_READ(config, WN_ZGAP_MAX, d, F, c, err);
      }

   } else {

      config->WCORR_SUBTYPE = wcorr_subtype_undef;

   }


   /* Number of data dimensions. Default dimension is 2, for wcorr_subtype_undef. */
   if (config->WCORR_SUBTYPE == nn_3d || config->WCORR_SUBTYPE == nn_rp_pi) *ndim = 3;
   else *ndim = 2;

   /* Number of bin dimensions. Can be lower than data dimensions. Default is 1. */
   if (config->WCORR_SUBTYPE == nn_rp_pi) *ndimbin = 2;
   else *ndimbin = 1;

   if (*ndim > 2 && config->RADEC == 1) {
      fprintf(stderr, "Mode 'radec' not possible with more than two data dimensions.");
      assert(0);
   }


   /* Exit with an error if WCORR is out of range */
   testErrorRet(config->WCORR<1 || config->WCORR>gg+gn+nn+wn, tree_outofrange,
         "WCORR (in config file) out of range", *err, __LINE__,);
   //testErrorRet(config->WCORR==3, tree_outofrange,
         //"WCORR=1+2 not yet supported, first fix weight bug!!", *err, __LINE__,);

   CONFIG_READ_S(config, SFORMAT, s, F, c, err);
   STRING2ENUM(config->FORMAT, config->SFORMAT, format_t, sformat_t, j, Nformat_t, err);

   if (config->FORMAT == f_fits) {
      CONFIG_READ(config, NCOL, i, F, c, err);
      if (config->NCOL > 0) {
         config->COL_NAMES = malloc_err(sizeof(char*) * config->NCOL, err);
         forwardError(*err, __LINE__,);
         for (i=0; i<config->NCOL; i++) {
            config->COL_NAMES[i] = malloc_err(sizeof(char) * 128, err);
            forwardError(*err, __LINE__,);
         }
         CONFIG_READ_S_ARR(config, COL_NAMES, s, i, config->NCOL, F, c, err);
      }
   } else {
      config->NCOL = 0;
   }
   if (config->NCOL <= 0) {
      config->COL_NAMES = NULL;
   }

   CONFIG_READ_S(config, SCOORD_INPUT, s, F, c, err);
   STRING2ENUM(config->COORD_INPUT, config->SCOORD_INPUT,  coord_t, scoord_t, j, Ncoord_t, err);

   CONFIG_READ_S(config, SCOORD_OUTPUT, s, F, c, err);
   STRING2ENUM(config->COORD_OUTPUT, config->SCOORD_OUTPUT,  coord_t, scoord_t, j, Ncoord_t, err);

   for (k=0; k<*ndimbin; k++) {
      CONFIG_READ(config, THMIN, d, F, c, err);
      CONFIG_READ(config, THMAX, d, F, c, err);
      CONFIG_READ(config, NTH, i, F, c, err);

      config->thmin[k] = config->THMIN;
      config->thmax[k] = config->THMAX;
      config->nth[k]   = config->NTH;

      coordinate_to_rad(&config->thmin[k], config->COORD_OUTPUT, err);  forwardError(*err, __LINE__,);
      coordinate_to_rad(&config->thmax[k], config->COORD_OUTPUT, err);  forwardError(*err, __LINE__,);      
   }

   /* Dummy variables */
   config->THMIN = config->THMAX = -1.0;
   config->NTH = 0;

   CONFIG_READ_S(config, BINTYPE, s, F, c, err);
   CONFIG_READ(config, RADEC, i, F, c, err);
   CONFIG_READ(config, OATH, d, F, c, err);
   CONFIG_READ_S(config, SERROR, s, F, c, err);

   STRING2ENUM(config->ERROR, config->SERROR, error_t, serror_t, j, Nerror_t, err);
   if (config->ERROR != error_none) {
      // MKDEBUG NEW v1.55: If single catalogue, only read one NRESAMPLE value
      CONFIG_READ_ARR(config, NRESAMPLE, i, i, config->NCAT, s, F, c, err);
   } else {
      config->NRESAMPLE[0] = 0;
   }
   if (config->NCAT == 1) {
      config->NRESAMPLE[1] = config->NRESAMPLE[0];
   }

   /* Bug in v1.6, fixed in v1.61 (added 'if' condition) */
   if (config->ERROR != error_none) {
      testErrorRet(config->NRESAMPLE[0] == 0 || config->NRESAMPLE[1] == 0, tree_nboots,
            "The number of resamples has to be larger than zero for SERROR != 'none'", *err, __LINE__,);
   }

   testErrorRet(config->NRESAMPLE[0]!=config->NRESAMPLE[1] &&
         (config->NRESAMPLE[0]!=0 && config->NRESAMPLE[1]!=0), tree_nboots,
         "The resample sizes for both catalogues (NRESAMPLE in config file) can not be different "
         "if one of them is not zero", *err, __LINE__,);

   testErrorRet(config->NCAT==1 && config->NRESAMPLE[0]!=config->NRESAMPLE[1], tree_nboots,
         "With a single catalogue the resample sizes have to be equal", *err, __LINE__,);

   for (k=0; k<*ndimbin; k++) {
      testErrorRet(config->thmin[k] < 0 || config->thmax[k] < 0 ||
            config->thmax[k] <= config->thmin[k], tree_outofrange,
            "Angular scales (in config file) out of range", *err, __LINE__,);
   }
   testErrorRet(config->OATH < 0, tree_outofrange, "OATH (in config file) out of range", *err, __LINE__,);

   if (strcmp(config->BINTYPE, "LIN")==0) {
      config->I_BINTYPE = LIN;
   } else if (strcmp(config->BINTYPE, "LOG")==0) {
      config->I_BINTYPE = LOG;
   } else {
      *err = addError(tree_bintype, "BINTYPE (in config file) has to be 'LIN' or 'LOG'", *err, __LINE__);
      return;
   }

   testErrorRet(config->RADEC!=0 && config->RADEC!=1, tree_radec, "RADEC (in config file) has to be 0 or 1",
         *err, __LINE__,);

   fileclose(F);

   /* Items not read from the config file */
   config->do_wtheta = 0;
}


/* ============================================================ *
 * Prints the information read from the config file into the    *
 * file LOG.
 * ============================================================ */
void out_config(FILE *LOG, const char *cname, config_info config, int ndimbin)
{
   int k;

   fprintf(LOG, "read from config file (%s): galcat1=%s galcat2=%s (ncat=%d) wcorr=%d (s)format=(%s)%d ncol=%d\n",
         cname, config.GALCAT1, config.GALCAT2, config.NCAT, (int)config.WCORR, config.SFORMAT, (int)config.FORMAT,
         config.NCOL);

   if (config.NCOL >= 0) {
      fprintf(LOG, "COL_NAMES =");
      for (k=0; k<config.NCOL; k++) {
         fprintf(LOG, " %s", config.COL_NAMES[k]);
      }
      fprintf(LOG, "\n");
   }

   fprintf(LOG, "(s)coord_input=(%s)%d " "(s)coord_output=(%s)%d\n",
         config.SCOORD_INPUT, (int)config.COORD_INPUT, config.SCOORD_OUTPUT, (int)config.COORD_OUTPUT);


   for (k=0; k<ndimbin; k++) {
      fprintf(LOG, "thmin[%d]=%g (%g') thmax[%d]=%g (%g') Nth[%d]=%u\n",
            k, config.thmin[k], config.thmin[k]/ARCMIN, k, config.thmax[k], config.thmax[k]/ARCMIN, k, config.nth[k]);
   }

   fprintf(LOG, "type=%s(%d) oath=%.3f serror=%s(%d) nresample=(%d %d)\n",
         config.BINTYPE, (int)config.I_BINTYPE, config.OATH, serror_t(config.ERROR), (int)config.ERROR,
         config.NRESAMPLE[0], config.NRESAMPLE[1]);

   fflush(LOG);

   if (config.WCORR==nn || config.WCORR==wn) {
      fprintf(LOG, "  (s)wcorr_subtype=(%s)%d", swcorr_subtype_t(config.WCORR_SUBTYPE), (int)config.WCORR_SUBTYPE);
      if (config.WCORR_SUBTYPE==wn_pairwise) fprintf(LOG, " wn_zgap_min=%g wn_zgap_max=%g", config.WN_ZGAP_MIN,
            config.WN_ZGAP_MAX);
      fprintf(LOG, "\n");
   }
}

/* ============================================================ *
 * Prints message 'str' to stderr and to log file LOG.		*
 * ============================================================ */
void print_err_log_func(FILE *LOG, int my_quiet, const char *str)
{
   fprintf(LOG, "%s", str);
   fflush(LOG);

   if (!my_quiet) {
      fprintf(stderr, "%s", str);
      fflush(stderr);
   }
}

/* ============================================================ *
 * The following three functions are used when growing the      *
 * tree. They are specific to the type galaxy_data.		*
 * ============================================================ */
void init_tree(node **root1, node **root2, galaxy_data *gal1, galaxy_data *gal2,
      uint ngal1, uint ngal2, int ndim, int twins_merge, const config_info *config, int quiet, FILE *LOG)
{
   fill_function *node_fill_galaxy;

   if (config->RADEC==0) node_fill_galaxy = node_fill_galaxy_xy;
   else node_fill_galaxy = node_fill_galaxy_radec;
 

   /* Create tree for first catalogue */
   if (!quiet) fprintf(stderr, "Growing tree 1\n");
   *root1 = grow_tree(gal1, 0, ngal1, ndim, config->NRESAMPLE[0], twins_merge, quiet, node_fill_galaxy, swap_galaxy,
		      get_pos_galaxy, out_galaxy);
   some_stats(1, ndim, *root1, config->RADEC, LOG);

   /* If there is second catalogue create second tree */
   if (config->NCAT==2) {
      if (!quiet) fprintf(stderr, "Growing tree 2\n");
      *root2 = grow_tree(gal2, 0, ngal2, ndim, config->NRESAMPLE[1], twins_merge, quiet, node_fill_galaxy, swap_galaxy,
			 get_pos_galaxy, out_galaxy);
      some_stats(2, ndim, *root2, config->RADEC, LOG);
   } else {
      *root2 = *root1;
   }
}


/* ============================================================ *
 * Direct correlation of two nodes.				                   *
 * ============================================================ */

/* === For wn corrletaion === */
//const galaxy_data *gal1_TMP, *gal2_TMP; FILE *ZGAP_TMP;


void correl_two_nodes(const node *nd1, const node *nd2, const double distance, 
		      const bin_data *bin, const config_info *config, xi_data *xi)
{
   double x, y, dd, cos2phi=0.0, sin2phi=0.0, cos4phi=0.0, sin4phi=0.0,
          g11, g22, g12, g21,  ww=-1.0, logd, gt, gx;
   double sind=0.0, cosphi1=0.0, sinphi1=0.0, cosphi2=0.0, sinphi2=0.0, cos2phi1=0.0,
          sin2phi1=0.0, cos2phi2=0.0, sin2phi2=0.0, cosDeltaalpha=0.0, sinDeltaalpha=0.0;
   int b[NDIMBIN_MAX], b_index=-1, i;
   double nresample1, nresample2, cos2Deltaphi_p=0.0, cos2Deltaphi_m=0.0, sin2Deltaphi_p=0.0, sin2Deltaphi_m=0.0;



   /* Don't correlate twice if auto-correlation of single catalog */
   if (config->NCAT==1 && nd1->start < nd2->start) return;


   /* Range check for distance, and bin index calculation */
   if (bin->linlog == LIN) {

      if (config->WCORR_SUBTYPE != nn_rp_pi) {
         /* One bin dimension */
         if (distance < bin->min[0] || distance > bin->max[0]) {
            return;
         }
         b_index = distance2bin_index(distance, bin->min[0], bin->diff[0], bin->N[0], 1);
      } else {
         /* Two bin dimensions (first: xy, second: z) */
         assert(config->RADEC == 0);

         /* Distance in xy-plane */
         dd = sqrt( absnsqr(nd1->barycenter, nd2->barycenter, 2) );
         if (dd < bin->min[0] || dd > bin->max[0]) {
            return;
         }
         b[0] = distance2bin_index(dd, bin->min[0], bin->diff[0], bin->N[0], 1);

         /* Distance in z-direction */
         dd = fabs(nd1->barycenter[2] - nd2->barycenter[2]);
         if (dd < bin->min[1] || dd > bin->max[1]) {
            return;
         }
         b[1] = distance2bin_index(dd, bin->min[1], bin->diff[1], bin->N[1], 1);
      }

   } else if (bin->linlog == LOG) {

      if (config->WCORR_SUBTYPE != nn_rp_pi) {
         /* One bin dimension */
         logd = log(distance);
         if (logd < bin->logmin[0] || logd > bin->logmax[0]) {
            return;
         }
         b_index = distance2bin_index(logd, bin->logmin[0], bin->difflog[0], bin->N[0], 1);
      } else {
         /* Two bin dimensions (fist: xy, second z) */
         assert(config->RADEC == 0);

         /* Distance in xy-plane */
         logd = log(sqrt( absnsqr(nd1->barycenter, nd2->barycenter, 2) ));
         if (logd < bin->logmin[0] || logd > bin->logmax[0]) {
            return;
         }
         b[0] = distance2bin_index(logd, bin->logmin[0], bin->difflog[0], bin->N[0], 1);

         /* Distance in z-direction */
         logd = log(fabs(nd1->barycenter[2] - nd2->barycenter[2]));
         if (logd < bin->logmin[1] || logd > bin->logmax[1]) {
            return;
         }
         b[1] = distance2bin_index(logd, bin->logmin[1], bin->difflog[1], bin->N[1], 1);

      }

   } else {
      assert(0);
   }


   /* Consistency check, nodes should not be identical here. Has to be after distance range check. */
   assert(!(nd1==nd2 && nd1->start==nd2->start));


   if (config->WCORR_SUBTYPE == nn_rp_pi) {
      /* Multi-dimensional bin index for 1d array */
      b_index = get_bin_index(b, bin->N, bin->ndim);
   }

   assert(b_index >= 0);


   if (config->WCORR&wn && config->WCORR_SUBTYPE==wn_pairwise) {
      /* Return if no redshifts were assigned to galaxies in node */
      assert(0);
      //if (nd2->zmax<0 || nd1->zmin<0) return; // MKDBEUG commented

      /* Return if foreground node maximum redshift too close to *
       * background node minimum redshift.			 */
      //if (nd2->zmax + config->WN_ZGAP_MIN > nd1->zmin) return; // MKDEBUG commented

      /* Return if foreground max redshift too far away from bg min redshift */
      //if (config->WN_ZGAP_MAX>0 && gal2_TMP[nd2->start].z + config->WN_ZGAP_MAX < gal1_TMP[nd1->start].z) return;

      // For wn correlation
      //assert(nd2->start==nd2->end-1);
      //assert(nd1->start==nd1->end-1);
      //fprintf(ZGAP_TMP, "%.3f %.3f  %.3f %.3f\n", gal2_TMP[nd2->start].z, gal1_TMP[nd1->start].z,
      //      gal2_TMP[nd2->start].weight, gal1_TMP[nd1->start].weight);

   }



   /* Number of pairs is product of number of galaxies in node 1 and 2 */
   xi->npair[b_index] += (unsigned long long int)nd1->ngal*(unsigned long long int)nd2->ngal;


   /* Resampled number of objects (bootstrap; jackknife) */
   if (config->WCORR&nn || config->WCORR&wn) {
      for (i=0; i<IMAX(config->NRESAMPLE[0], config->NRESAMPLE[1]); i++) {
         if (config->NRESAMPLE[0]>0) nresample1 = nd1->ngal_resample[i];
         else nresample1 = (double)nd1->ngal;

         if (config->NRESAMPLE[1]>0) nresample2 = nd2->ngal_resample[i];
         else nresample2 = (double)nd2->ngal;

         xi->npair_resample[b_index][i] += nresample1*nresample2;
      }
   }


   /* Distances and angles */
   if (config->RADEC==0) {

      if (config->WCORR&gg || config->WCORR&gn) {
         dd = DSQR(distance);
         x  = nd1->barycenter[0] - nd2->barycenter[0];
         y  = nd1->barycenter[1] - nd2->barycenter[1];

         cos2phi = 2.0*x*x/dd - 1.0;
         sin2phi = 2.0*x*y/dd;     

         if (config->WCORR&gg) {
            cos4phi = 2.0*DSQR(cos2phi) - 1.0;
            sin4phi = 2.0*sin2phi*cos2phi;
         }
      }

   } else if (config->RADEC >= 1) {

      if (config->WCORR&gg ||config->WCORR&gn) {

         /* New in version 1.5: two course angles, to take into account    *
          * geodesic difference for projection at the two galaxy positions */

         sind    = sin(distance);
         cosDeltaalpha = nd1->cosbarycenter[0] * nd2->cosbarycenter[0] + nd1->sinbarycenter[0] * nd2->sinbarycenter[0];
         sinDeltaalpha = nd1->cosbarycenter[0] * nd2->sinbarycenter[0] - nd1->sinbarycenter[0] * nd2->cosbarycenter[0];

         cosphi1  = sinDeltaalpha * nd2->cosbarycenter[1] / sind;
         sinphi1  = (nd1->cosbarycenter[1] * nd2->sinbarycenter[1] - nd1->sinbarycenter[1] * nd2->cosbarycenter[1] * cosDeltaalpha) / sind;
         cos2phi1 = 2.0*DSQR(cosphi1) - 1.0;
         sin2phi1 = 2.0*cosphi1*sinphi1;

         cosphi2  = -sinDeltaalpha * nd1->cosbarycenter[1] / sind;
         sinphi2  = (nd2->cosbarycenter[1] * nd1->sinbarycenter[1] - nd2->sinbarycenter[1] * nd1->cosbarycenter[1] * cosDeltaalpha) / sind;
         cos2phi2 = 2.0*DSQR(cosphi2) - 1.0;
         sin2phi2 = 2.0*cosphi2*sinphi2;

      }

   } else {
      /* Wrong RADEC mode */
      assert(0);
   }

   /* Shear-shear correlation (cosmic shear) */
   if (config->WCORR&gg) {

      /* Shear two-point correlators g_i g_j */
      g11 = nd1->well.r * nd2->well.r;
      g22 = nd1->well.i * nd2->well.i;
      g12 = nd1->well.r * nd2->well.i;
      g21 = nd1->well.i * nd2->well.r;

      if (config->RADEC == 1) {

         /* For xi+: 2*Deltaphi = 2phi1 - 2phi2 */
         cos2Deltaphi_p  = cos2phi1 * cos2phi2 + sin2phi1 * sin2phi2;
         sin2Deltaphi_p  = sin2phi1 * cos2phi2 - cos2phi1 * sin2phi2;

         /* For xi- and xix: 2*Deltaphi = 2phi1 + 2phi2 (earlier: 4phi) */
         cos2Deltaphi_m  = cos2phi1 * cos2phi2 - sin2phi1 * sin2phi2;
         sin2Deltaphi_m  = sin2phi1 * cos2phi2 + cos2phi1 * sin2phi2;

         xi->p[b_index] += (g11 + g22) * cos2Deltaphi_p + (g12 - g21) * sin2Deltaphi_p;
         xi->m[b_index] += (g11 - g22) * cos2Deltaphi_m + (g12 + g21) * sin2Deltaphi_m;
         xi->x[b_index] += 0.5 * ( (-g11 + g22) * sin2Deltaphi_m + (g12 + g21) * cos2Deltaphi_m );

      } else {

         /* xi_+ = <gamma_t gamma_t> + <gamma_x gamma_x> */
         xi->p[b_index] += g11 + g22;

         /* xi_- = <gamma_t gamma_t> - <gamma_x gamma_x> */
         xi->m[b_index] += (g11 - g22) * cos4phi + (g12 + g21) * sin4phi;

         /* xi_x = <gamma_t gamma_x> */
         xi->x[b_index] += 0.5*( (-g11 + g22) * sin4phi + (g12 + g21) * cos4phi );

      }

      /* Reference-frame-dependent quantities (for systematics check) */
      xi->g11[b_index] += g11;
      xi->g22[b_index] += g22;
      xi->g12[b_index] += g12;
      xi->g21[b_index] += g21;

      /* Sum of weights, for Poisson error */
      ww = nd1->weight * nd2->weight;              /* ww = w_i * w_j */
      xi->ww[b_index]   += ww;                       /* Sum_ij[ w_i w_j ] = Npair */
      xi->wwww[b_index] += ww*ww;                    /* Sum_ij[ (w_i w_j)^2 ] */
   }

   if (config->WCORR&gg || config->WCORR&gn) {
      xi->Nnode[b_index] ++;
   }

   /* Shear-position correlation (galaxy-galaxy lensing) */
   if (config->WCORR&gn) {

      /* Bug fix in v.1.53: Added missing weighting with number of foreground objects, *
       * this caused underestimation of gt for oath>0.                                 */

      /* nd1 is background, nd2 foreground catalogue */
      if (config->RADEC == 0) {
          gt = - ( nd1->well.r * cos2phi + nd1->well.i * sin2phi) * nd2->ngal;
          gx = - (-nd1->well.r * sin2phi + nd1->well.i * cos2phi) * nd2->ngal;
          xi->gt[b_index]   += gt;
          xi->gx[b_index]   += gx;
      } else {
          gt = - ( nd1->well.r * cos2phi2 + nd1->well.i * sin2phi2) * nd2->ngal;
          gx = - (-nd1->well.r * sin2phi2 + nd1->well.i * cos2phi2) * nd2->ngal;
          xi->gt[b_index]   += gt;
          xi->gx[b_index]   += gx;  
      }
   }

   if (config->WCORR & gn || config->WCORR & wn) {
      xi->wn[b_index]   += nd1->weight * nd2->ngal;
      xi->wnwn[b_index] += DSQR(nd1->weight * nd2->ngal);
   }


   /* Weighted scale */
   if (config->do_wtheta) {

      if (config->WCORR & gg) {
         /* ww = w^2 already calculated */
      } else if (config->WCORR & gn || config->WCORR & wn) {
         ww = nd1->weight * nd2->ngal;
      } else {
         ww = 1.0;         /* No weighting for nn */
      }

      assert(ww > 0); 

      if (config->WCORR_SUBTYPE != nn_rp_pi) {
         xi->wtheta[b_index] += distance * ww;
         //if (b_index == 3) fprintf(stderr, "%g %g %g\n", distance/ARCMIN, ww, ww * distance/ARCMIN);
      } else {
         xi->wtheta[b[0]] += sqrt( absnsqr(nd1->barycenter, nd2->barycenter, 2) ) * ww;
         xi->wtheta[b[1]] += fabs(nd1->barycenter[2] - nd2->barycenter[2]) * ww;
      }
   }


   /* Resampled shear */
   // MKDEBUG: Put resampled numbers here from above
   if (config->WCORR & gg) {

      for (i=0; i<IMAX(config->NRESAMPLE[0], config->NRESAMPLE[1]); i++) {

         /* MKDEBUG TODO test: probably, both nresample values have to be equal */
         g11 = nd1->well_resample[i].r * nd2->well_resample[i].r;
         g22 = nd1->well_resample[i].i * nd2->well_resample[i].i;
         g12 = nd1->well_resample[i].r * nd2->well_resample[i].i;
         g21 = nd1->well_resample[i].i * nd2->well_resample[i].r;

         if (config->RADEC == 0) {
            xi->p_resample[b_index][i] += g11 + g22;
            xi->m_resample[b_index][i] += (g11 - g22) * cos4phi + (g12 + g21) * sin4phi;
         } else {
            xi->p_resample[b_index][i] += (g11 + g22) * cos2Deltaphi_p + (g12 - g21) * sin2Deltaphi_p;
            xi->m_resample[b_index][i] += (g11 - g22) * cos2Deltaphi_m + (g12 + g21) * sin2Deltaphi_m;
         }

         /* Angles and distances corresponds to original barycenters, not exactly = resampled ones */
         xi->ww_resample[b_index][i] += nd1->weight_resample[i] * nd2->weight_resample[i];

      }
   } else if (config->WCORR & gn) {

      for (i=0; i<IMAX(config->NRESAMPLE[0], config->NRESAMPLE[1]); i++) {

         if (config->RADEC == 0) {
            gt = - (nd1->well_resample[i].r * cos2phi + nd1->well_resample[i].i * sin2phi) * nd2->ngal_resample[i];
            xi->gt_resample[b_index][i] += gt;
            xi->wn_resample[b_index][i] += nd1->weight_resample[i] * nd2->ngal_resample[i];
         } else {
            gt = - (nd1->well_resample[i].r * cos2phi2 + nd1->well_resample[i].i * sin2phi2) * nd2->ngal_resample[i];
            xi->gt_resample[b_index][i] += gt;
            xi->wn_resample[b_index][i] += nd1->weight_resample[i] * nd2->ngal_resample[i];
         }

      }

   }

}


/* ============================================================ *
 * Called recursively. According to an open-angle criterium,    *
 * splits or correlates nodes.					                   *
 * No error messaging because of recursive calls.		          *
 * ============================================================ */
#define EPSILON 1.0e-6

void twopcf(const galaxy_data *gal1, const galaxy_data *gal2, uint ngal1, uint ngal2,
	    int ndim, node *nd1, node *nd2, const bin_data *bin, 
	    const config_info *config, xi_data *xi)
{
   double oa1=0.0, oa2=-1.0, distance=-1.0, OATHd=0.0;
   int split1, split2;

   split1 = split2 = 0;

   assert(nd1);
   assert(nd2);

   if (nd1->ngal>1) {

      /* Calculates open angle between nodes */
      oa1 = open_angle(nd2, nd1, &distance, ndim, config->RADEC);
 
  } else {

      /* Distance between 2 centers of the nodes */
      if (config->RADEC==0) {

	 distance = sqrt( absnsqr(nd1->barycenter, nd2->barycenter, ndim) );

      } else {

	 /* Distance in spherical coordinates, 2d only.                   *
	  * cos(distance) =  cos(a1-a2) cos(d1) cos(d2) + sin(s1) sin(d2) */
	 distance = (nd1->cosbarycenter[0]*nd2->cosbarycenter[0] +
		     nd1->sinbarycenter[0]*nd2->sinbarycenter[0])
	   * nd1->cosbarycenter[1]*nd2->cosbarycenter[1]
	   + nd1->sinbarycenter[1]*nd2->sinbarycenter[1];

	 if (fabs(distance)>1.0+EPSILON) {
	    fprintf(stderr, "Trying to take acos of %.10f!\n", distance);
	    assert(0);
	 } else if (distance>1.0) {
	    distance = 0.0;   /* acos(1) */
	 } else if (distance<-1.0) {
	    distance = pi;    /* acos(-1) */
	 } else {
	    distance = acos(distance);
	 }
      }

   }

   OATHd = config->OATH;

   if (nd1->left!=NULL && nd1->right!=NULL && oa1>OATHd) split1 = 1;

   if (nd2->left!=NULL && nd2->right!=NULL) {
      if (nd2->ngal > 1) {
	 oa2 = open_angle(nd1, nd2, &distance, ndim, config->RADEC);
      }
      assert(oa2 >= 0.0); // Nnew
      if (oa2>OATHd) split2 = 1;
   }


   switch (2 * split1 + split2) {

      case 3 :
	 /* 2 + 1: Split 1 and 2 */
	 twopcf(gal1, gal2, ngal1, ngal2, ndim, nd1->left, nd2->left, bin, config, xi);
	 twopcf(gal1, gal2, ngal1, ngal2, ndim, nd1->left, nd2->right, bin, config, xi);
	 twopcf(gal1, gal2, ngal1, ngal2, ndim, nd1->right, nd2->left, bin, config, xi);
	 twopcf(gal1, gal2, ngal1, ngal2, ndim, nd1->right, nd2->right, bin, config, xi);
	 return;
      case 2 :
	 /* 2 + 0: Split 1 */
	 twopcf(gal1, gal2, ngal1, ngal2, ndim, nd1->left, nd2, bin, config, xi);
	 twopcf(gal1, gal2, ngal1, ngal2, ndim, nd1->right, nd2, bin, config, xi);
	 return;
      case 1 :
	 /* 0 + 1: Split 2 */
	 twopcf(gal1, gal2, ngal1, ngal2, ndim, nd1, nd2->left, bin, config, xi);
	 twopcf(gal1, gal2, ngal1, ngal2, ndim, nd1, nd2->right, bin, config, xi);
	 return;
      case 0 :
	 /* 0 + 0: Don't split */

	 /* Progress indicator */
	 if (! quiet_gl) {
	    f += (long double)nd1->ngal * (long double)nd2->ngal/(long double)npairtot;
	    if (f - p > 0.0001) {
	       fprintf(stderr, "~ %5.1f%%                          \r", (double)f*100.0); fflush(stderr);
	       p = f;
	    }
	 }
	 break;
      default :
	 assert(0);
   }


   /* Direct correlation */
   correl_two_nodes(nd1, nd2, distance, bin, config, xi);
}

#undef EPSILON

void fill_bootstrap_sample(int NRESAMPLE, wcorr_t WCORR, galaxy_data *gal, uint ngal, int quiet, error **err)
{
   int j, total, num;

   if (!quiet) fprintf(stderr, "Creating bootstrap sample");

   allocate_jackboot_sample(NRESAMPLE, gal, ngal, err);
   forwardError(*err, __LINE__,);

   for (j=0; j<NRESAMPLE; j++) {
      if (!quiet) { fprintf(stderr,"."); fflush(stderr); }
      total = 0;
      do {
         /* Draw random galaxy numbers, increase corresponding slot by 1 (or by weight for wn) */
         num = (int)(rand()/((double)RAND_MAX+1.0)*ngal);
         if (WCORR&wn) gal[num].ind_resample[j] += gal[num].weight;
         else          gal[num].ind_resample[j]++;
      } while (++total<ngal);
   }

   if (!quiet) { fprintf(stderr, "done\n"); fflush(stderr); }
}

int update_nresample(int *NRESAMPLE, int i, int njack, int quiet)
{
   int nboots0_prev;

   if (!quiet) {
      if (NRESAMPLE[i] != njack) {
         fprintf(stderr, "Njackknife for cat %d updated from %d to %d\n", i+1, NRESAMPLE[i], njack);
      } else {
         fprintf(stderr, "Njackknife for cat %d left to be %d\n", i+1, njack);
      }
   }
   nboots0_prev = NRESAMPLE[i];
   NRESAMPLE[i] = njack;

   return nboots0_prev;
}

/* Bootstrap index */
void init_jackboot(config_info *config, galaxy_data *gal1, galaxy_data *gal2, uint ngal1, uint ngal2,
                   int quiet, error **err)
{
   int j, njack, nboots0_prev=-1;
   double min[2], max[2];

   if (config->NRESAMPLE[0]>0) {

      switch (config->ERROR) {
         case bootstrap :

            fill_bootstrap_sample(config->NRESAMPLE[0], config->WCORR, gal1, ngal1, quiet, err);
            forwardError(*err, __LINE__,);
            break;

         case jackknife :

            switch (config->FORMAT) {
               case f_position_jack_num : case f_lensing_jack_num :
                  njack = fill_jackknife_sample_idx(config->NRESAMPLE[0], gal1, ngal1, err);
                  forwardError(*err, __LINE__,);
                  break;
               default :
                  get_extend_2d(gal1, ngal1, gal2, ngal2, min, max);
                  njack = fill_jackknife_sample(config->NRESAMPLE[0], gal1, ngal1, min, max, err);
                  forwardError(*err, __LINE__,);
                  break;
            }

            nboots0_prev = update_nresample(config->NRESAMPLE, 0, njack, quiet);
            if (config->NCAT == 1) config->NRESAMPLE[1] = njack; /* For consistency */

            break;

         default:

            *err = addError(tree_error_t, "Wrong error type", *err, __LINE__);
            return;

      }

   } else {

      for (j=0; j<ngal1; j++) {
         gal1[j].ind_resample = NULL;
      }

   }

   if (config->NCAT == 2) {

      //assert(nboots0_prev > 0); // New
      testErrorRetVA(config->ERROR == jackknife && nboots0_prev != config->NRESAMPLE[1], tree_jack_NJ,
            "Both catalogues have to have the same Jackknife sample size (found: %d, %d)", *err, __LINE__,,
            config->NRESAMPLE[0], config->NRESAMPLE[1]);

      if (config->NRESAMPLE[1]>0) {

         switch (config->ERROR) {

            case bootstrap :
               fill_bootstrap_sample(config->NRESAMPLE[1], config->WCORR, gal2, ngal2, quiet, err);
               forwardError(*err, __LINE__,);
               break;

            case jackknife :

               switch (config->FORMAT) {
                  case f_position_jack_num : case f_lensing_jack_num :
                     /* v1.54: Bug fixed, gal1 instead of gal2 was used, led to segmenation fault */
                     njack = fill_jackknife_sample_idx(config->NRESAMPLE[0], gal2, ngal2, err);
                     forwardError(*err, __LINE__,);
                     break;
                  default :
                     /* NEW: Use min, max from catalogue #1 to get consistent jackknife regions */
                     njack = fill_jackknife_sample(config->NRESAMPLE[1], gal2, ngal2, min, max, err);
                     forwardError(*err, __LINE__,);
                     break;
               }

               update_nresample(config->NRESAMPLE, 1, njack, quiet);
               break;

            default:

               *err = addError(tree_error_t, "Wrong error type", *err, __LINE__);
               return;

         }

      } else {

         /* NEW: 2 catalogues, NRESAMPLE[1] is zero */
         for (j=0; j<ngal2; j++) {
            gal2[j].ind_resample = NULL;
         }

      }

   } else {
      /* Do nothing: gal2 points to gal1 in this case */
   }

}

/* Allocates NRESAMPLE slots for each of the ngal galaxies */
void allocate_jackboot_sample(int NRESAMPLE, galaxy_data *gal, uint ngal, error **err)
{
   int j;

   for (j=0; j<ngal; j++) {
      gal[j].ind_resample = calloc_err(NRESAMPLE, sizeof(double), err);
      forwardError(*err, __LINE__,);
   }
}

/* Calculates integers Nx, Ny with Nx * Ny = N, and Ny/Nx = y/x = r */
void Nxy_jackknife(uint N, double r, uint *Nx, uint *Ny)
{
   double dNx, dNy;

   dNx = sqrt((double)N / r);
   dNy = sqrt((double)N * r);

   *Nx  = (uint)round(dNx);
   *Ny  = (uint)round(dNy);

   if (*Nx == 0) (*Nx) ++;
   if (*Ny == 0) (*Ny) ++;

   //printf("%g %g %u %u  %u\n", dNx, dNy, *Nx, *Ny, N);
}

/* ============================================================ *
 * Sets the jackknife index array to values (1 or 0), corres-   *
 * ponding to the Jackknife sample number.			*
 * ============================================================ */
int fill_jackknife_sample(int NJACK, galaxy_data *gal, uint ngal, const double *min, const double *max, error **err)
{
   double r, Dx, Dy, eps;
   uint Nx, Ny, j;
   int b, jx, jy, jind;

   /* Ratio of y to x extend. TODO: more than 2 dimensions. */
   eps = 1.0e-6;
   Dy = (max[1] - min[1]) * (1.0 + eps);
   Dx = (max[0] - min[0]) * (1.0 + eps);
   r =  Dy / Dx; 

   Nxy_jackknife(NJACK, r, &Nx, &Ny);

   allocate_jackboot_sample(Nx * Ny, gal, ngal, err);
   forwardError(*err, __LINE__, -1);

   for (j=0; j<ngal; j++) {

      jx   = (int)((gal[j].pos[0] - min[0]) / Dx * (double)Nx);
      jy   = (int)((gal[j].pos[1] - min[1]) / Dy * (double)Ny);
      jind = jx + jy * Nx;

      testErrorRetVA(jind < 0, tree_jack_NJ,
		     "Jackknife sample index is -1 for gal #%d.", *err, __LINE__, -1, j);
      testErrorRetVA(jind >= Nx * Ny, tree_jack_NJ,
		     "Jackknife sample index %d (%d, %d) >= NJACK=%d from config file, for gal #%d",
		     *err, __LINE__, -1, jind, jx, jy, Nx * Ny, j);

      for (b=0; b<Nx * Ny; b++) {
         if (b != jind) gal[j].ind_resample[b] = 1.0;
      }
   }

   return Nx * Ny;
}

/* ============================================================ *
 * Fills Jackknife index array according to Jackknife sample    *
 * number gal->idx. NJACK has to be larger or equal to the      *
 * number of distinct idx values.				                   *
 * ============================================================ */
int fill_jackknife_sample_idx(int NJACK, galaxy_data *gal, uint ngal, error **err)
{
   uint j;
   int b, NJmax;
   unsigned short int jind;

   allocate_jackboot_sample(NJACK, gal, ngal, err);
   forwardError(*err, __LINE__, -1);

   for (j=0,NJmax=-1; j<ngal; j++) {

      jind = gal[j].idx; 
      if (jind > NJmax) NJmax = jind;

      testErrorRetVA(NJmax >= NJACK, tree_jack_NJ,
		     "Jackknife sample index %d larger than NJACK=%d from config file, for gal #%d",
		     *err, __LINE__, -1, jind, NJACK, j);

      for (b=0; b<NJACK; b++) {
         if (b != jind) gal[j].ind_resample[b] = 1.0;
      }
   }

   /* Return one plus the largest index */
   return NJmax + 1;
}

/* ============================================================ *
 * Calls the sub-routines which write the correlation functions *
 * to files.							*
 * ============================================================ */
void write_all_correlations(const xi_data *xi, const char *names[], const bin_data *bin, const config_info *config,
			    double sig_eps2, double sig_eps4, double wbar, uint ngal1, uint ngal2, error **err)
{
   int imax;
   char name_resample[128];

   imax = IMAX(config->NRESAMPLE[0], config->NRESAMPLE[1]);

   if (config->WCORR&gg) {
      out_xi(xi, names[ID_SHEAR_CORR], names[ID_SHEAR_CORR_REF], names[ID_SHEAR_CORR2], bin, sig_eps4,
	          config->NCAT, config->COORD_OUTPUT, err);
      forwardError(*err, __LINE__,);

      if (config->ERROR != error_none) {
         out_xi_resample(xi, names[ID_SHEAR_RESAMPLE_CORR], names[ID_SHEAR_ALL_RESAMPLE_CORR_XIP],
		                   names[ID_SHEAR_ALL_RESAMPLE_CORR_XIM], bin, config->COORD_OUTPUT,
			                imax, config->ERROR, err);
         forwardError(*err, __LINE__,);

         sprintf(name_resample, "%s.cov", names[ID_SHEAR_RESAMPLE_CORR]);
         out_cov_xi_resample(xi, name_resample, bin, config->COORD_OUTPUT, imax, config->ERROR, err);
         forwardError(*err, __LINE__,);
      }

   }

   if (config->WCORR&nn) {
      switch (config->WCORR_SUBTYPE) {
         case nn_2d : case nn_3d :
            out_w(xi, names[ID_CORR], bin, ngal1, ngal2, config->COORD_OUTPUT, imax, config->ERROR, err);
            forwardError(*err, __LINE__,);
            break;
         case nn_rp_pi :
            out_wp_rp_pi(xi, names[ID_CORR], bin, imax, config->COORD_OUTPUT, ngal1, ngal2, err);
            forwardError(*err, __LINE__,);
            break;
         default :
            assert(0);
      }
   }

   if (config->WCORR&gn) {
      out_gl(xi, names[ID_GAL_SHEAR_XCORR], bin, sig_eps2, config->COORD_OUTPUT, err);
      forwardError(*err, __LINE__,);

      if (config->ERROR != error_none) {
         out_gl_resample(xi, names[ID_GAL_SHEAR_RESAMPLE_XCORR], bin,names[ID_GAL_SHEAR_ALL_RESAMPLE_CORR_GT],
                config->COORD_OUTPUT, imax, config->ERROR, err);
         forwardError(*err, __LINE__,);

         sprintf(name_resample, "%s.cov", names[ID_GAL_SHEAR_RESAMPLE_XCORR]);
         out_cov_gl_resample(xi, name_resample, bin, config->COORD_OUTPUT, imax, config->ERROR, err);
         forwardError(*err, __LINE__,);
      }
   }

   if (config->WCORR&wn) {
      out_ww(xi, names[4], bin, imax, config->COORD_OUTPUT, wbar, err);
      forwardError(*err, __LINE__,);
   }
}

#ifdef _WITH_FITS
void write_all_correlations_fits(const xi_data *xi, const char *names[], const bin_data *bin, const config_info *config,
      double sig_eps2, double sig_eps4, double wbar, uint ngal1, uint ngal2, int quiet, error **err)
{
   int imax;

   imax = IMAX(config->NRESAMPLE[0], config->NRESAMPLE[1]);

   if (config->WCORR & gg) {
      out_xi_fits(xi, names[ID_SHEAR_CORR], bin, sig_eps4, config->NCAT, config->COORD_OUTPUT, imax, config->ERROR, quiet, err);
      forwardError(*err, __LINE__,);
   }

   if (config->WCORR & gn) {
      out_gl_fits(xi, names[ID_GAL_SHEAR_XCORR], bin, sig_eps2, config->COORD_OUTPUT, imax, config->ERROR, quiet, err);
      forwardError(*err, __LINE__,);
   }

   if (config->WCORR & nn) {

      switch (config->WCORR_SUBTYPE) {
         case nn_2d : case nn_3d :
            out_w_fits(xi, names[ID_CORR], bin, ngal1, ngal2, config->COORD_OUTPUT, imax, config->ERROR, quiet, err);
            forwardError(*err, __LINE__,);
            break;
         default :
            assert(0);
      }

   }

}
#endif


/* ============================================================ *
 * If `str' is NULL, copies the content of 'name_to_set'.        *
 * ============================================================ */
void set_char_if_null(char **str, const char *name_to_set, error **err)
{
   if (*str==NULL) {
      *str = malloc_err(512*sizeof(char), err);
      forwardError(*err, __LINE__,);
      strcpy(*str, name_to_set);
    }
}


/* ============================================================ *
 * Print some catalogue stats to the log file.			*
 * ============================================================ */
void some_stats(int i, int ndim, const node *root, int radec, FILE *LOG)
{
   double d;

   if (radec == 0) {
      d = abs22(root->min, root->max);
   } else {
      d = cos(root->max[0] - root->min[0]) * cos(root->min[1]) * cos(root->max[1]) + sin(root->min[1]) * sin(root->max[1]);
      d = acos(d);
   }

   fprintf(LOG, "Galaxy catalogue #%d:\n", i);
   fprintf(LOG, "  radius = %g rad (%g arcmin)\n", root->radius, root->radius / ARCMIN);
   fprintf(LOG, "  largest separation = %g rad (%g arcmin)\n", d, d / ARCMIN);
   fprintf(LOG, "  area   = %g sq rad (%g sq arcmin)\n", volume(root, ndim), volume(root, ndim) / ARCMIN / ARCMIN);
   fprintf(LOG, "  number density = %g/[sq rad] (%g/[sq arcmin])\n", numberdensity(root, ndim),
	   numberdensity(root, ndim) * ARCMIN * ARCMIN);

   if (ndim > 2) {
      fprintf(LOG, "The above area and number density are in ndim=%d\n", ndim);
      fprintf(LOG, "Considering only the first two coordinates:\n");
      fprintf(LOG, "  area (2d) = %g sq rad (%g sq arcmin)\n", volume(root, 2), volume(root, 2) / ARCMIN / ARCMIN);
      fprintf(LOG, "  number density (2d) = %g/[sq rad] (%g/[sq arcmin])\n", numberdensity(root, 2),
	      numberdensity(root, 2) * ARCMIN * ARCMIN);
   }
}

void usage(int ex)
{
   fprintf(stderr, "Usage: athena [OPTIONS]\n");
   fprintf(stderr, "OPTIONS:\n");
   fprintf(stderr, "  -c, --config CONFIG        Configuration file (default: config_pmc)\n");
   fprintf(stderr, "  -t, --twins_nomerge        Do not merge identical twin objects into one node. If there are twin\n");
   fprintf(stderr, "                              objects in the input catalogue, 'athena -t' exits with an error meseage\n");
   fprintf(stderr, "  -w                         Weighted angular scales at output\n");
   fprintf(stderr, "  -s                         Show angular scales and exit\n");
   fprintf(stderr, "  --out_xi  NAME             Output file name NAME for shear-shear correlation\n");
   fprintf(stderr, "                              (default '%s')\n", SHEAR_CORR_NAME);
   fprintf(stderr, "  --out_gl NAME              Output file name NAME for shear-position correlation (default '%s')\n",
	   GAL_SHEAR_XCORR_NAME);
   fprintf(stderr, "  --out_w NAME               Output file name NAME for position correlation (default '%s')\n", CORR_NAME);
   fprintf(stderr, "  --out_xi_resample NAME     Output file name NAME for resampled shear-shear correlation (default '%s')\n",
	   SHEAR_RESAMPLE_CORR_NAME);
   fprintf(stderr, "  --out_gl_resample NAME     Output file name NAME for resampled shear-position correlation (default '%s')\n",
	   GAL_SHEAR_RESAMPLE_XCORR_NAME);
   fprintf(stderr, "  --out_xi2 NAME             Output file name NAME for shear-shear correlation with min and max of bins\n");
   fprintf(stderr, "                              (default '%s')\n", SHEAR_CORR_NAME2);
   fprintf(stderr, "  --out_xiref NAME           Output file name NAME for shear reference-frame-dependent correlation\n");
   fprintf(stderr, "                              (default: file is not written)\n");
   fprintf(stderr, "  --out_ww NAME              Output file name NAME for weighted position correlation (default '%s')\n",
	   WEIGHTED_ANGULAR_CORR_NAME);
   fprintf(stderr, "  --out_ALL_xip_resample NAME  Output file name NAME for All values for resamples xi_plus values \n");
   fprintf(stderr, "                              (default '%s')\n", SHEAR_ALL_RESAMPLE_CORR_NAME_XIP);
   fprintf(stderr, "  --out_ALL_xim_resample NAME  Output file name NAME for All values for resamples xi_minus values \n");
   fprintf(stderr, "                              (default '%s')\n", SHEAR_ALL_RESAMPLE_CORR_NAME_XIM);
   fprintf(stderr, "  --out_ALL_gt_resample NAME  Output file name NAME for All values for resamples gamma-t values \n");
   fprintf(stderr, "                              (default '%s')\n", ID_GAL_SHEAR_ALL_RESAMPLE_CORR_NAME_GT);
   fprintf(stderr, "  --out_log NAME             Log file name NAME (default '%s')\n", LOG_NAME);
   fprintf(stderr, "  --out_suf SUFFIX           SUFFIX is appended to all output file names\n");
   fprintf(stderr, "  -q, --quiet                Quiet mode\n");
   fprintf(stderr, "  -v, --version              Print version and exit\n");
   fprintf(stderr, "  -h, --help                 This message\n");

   if (ex>=0) exit(ex);
}

/* ============================================================ *
 * Main program.						*
 * ============================================================ */
int main(int argc, char *argv[])
{
   galaxy_data *gal1=NULL, *gal2=NULL;
   node *root1, *root2;
   bin_data *bin;
   xi_data *xi;
   time_t t_start;
   long double fnpairtot;
   FILE *LOG;
   config_info config;
   uint ngal1, ngal2;
   int b, c, twins_merge, ndim, ndimbin, do_wtheta, show_theta_and_exit, version_and_exit;
   char *cname, *names[NNAMES], str_argv[4096], *p_str;
   double sig_eps2, sig_eps4, wbar;
   error *myerr = NULL, **err;
   int quiet;

   err = &myerr;


   /* Copy command line arguments to string for log file */
   p_str = str_argv;
   for (c=0; c<argc; c++) {
      p_str = stpcpy(p_str, argv[c]);
      p_str = stpcpy(p_str, " ");

      // Alternatively (E. Jullo):
      //p_str = strcpy(p_str, argv[c]);
      //p_str = strcat(p_str, " ");
   }




   /* Command line arguments */
   cname = NULL;
   twins_merge = 1;
   do_wtheta = 0;
   quiet = 0;
   show_theta_and_exit = 0;
   version_and_exit = 0;
   for (b=0; b<NNAMES; b++) names[b] = NULL;
   while (1)
   {
      static struct option long_options[] =
      {
         {"config", required_argument, 0, 'c'},
         {"twins_merge", no_argument, 0, 't'},
         {"out_xi", required_argument, 0, 0},
         {"out_xiref", required_argument, 0, 1},
         {"out_w", required_argument, 0, 2},
         {"out_gl", required_argument, 0, 3},
         {"out_ww", required_argument, 0, 4},
         {"out_log", required_argument, 0, 5},
         {"out_xi2", required_argument, 0, 6},
         {"out_xi_resample", required_argument, 0, 7},
         {"out_ALL_xip_resample", required_argument, 0, 8},
         {"out_ALL_xim_resample", required_argument, 0, 9},
	      {"out_gl_resample", required_argument, 0, 10},
	      {"out_ALL_gt_resample", required_argument, 0, 11},
         {"out_suf", required_argument, 0, 12},
         {"version", no_argument, 0, 'v'},
         {"quiet", no_argument, 0, 'q'},
         {"help", no_argument, 0, 'h'},
         {0, 0, 0, 0}
      };

      int option_index = 0;

      c = getopt_long(argc, argv, "c:tswvqh", long_options, &option_index);

      switch (c) {
         case 'c' :
            cname = optarg;
            break;
         case 0 : case 1 : case 2 : case 3 : case 4 : case 5 : case 6 : case 7 : case 8 : case 9 : case 10 : case 11 :
         case 12 :
            names[c] = optarg;
            break;
         case 't' :
            twins_merge = 0;
            break;
         case 'w' :
            do_wtheta = 1;
            break;
         case 's' :
            show_theta_and_exit = 1;
            break;
         case 'v' :
            version_and_exit = 1;
            break;
         case 'q' :
            quiet = 1;
            break;
         case 'h' :
            usage(0);
         case -1  :
            goto end_arg;
         default :
            usage(2);
      }

   }
end_arg:

   set_char_if_null(&cname,  CONFIG_DEF_NAME, err);  exitOnError(*err,stderr);
   set_char_if_null(names + ID_SHEAR_CORR, SHEAR_CORR_NAME, err);  exitOnError(*err,stderr);
   set_char_if_null(names + ID_SHEAR_CORR_REF, SHEAR_CORR_REF_NAME, err);  exitOnError(*err,stderr);
   set_char_if_null(names + ID_CORR, CORR_NAME, err);  exitOnError(*err,stderr);
   set_char_if_null(names + ID_GAL_SHEAR_XCORR, GAL_SHEAR_XCORR_NAME, err);  exitOnError(*err,stderr);
   set_char_if_null(names + ID_WEIGHTED_ANGULAR_CORR, WEIGHTED_ANGULAR_CORR_NAME, err);  exitOnError(*err,stderr);
   set_char_if_null(names + ID_LOG, LOG_NAME, err);  exitOnError(*err,stderr);
   set_char_if_null(names + ID_SHEAR_CORR2, SHEAR_CORR_NAME2, err);  exitOnError(*err,stderr);
   set_char_if_null(names + ID_SHEAR_RESAMPLE_CORR, SHEAR_RESAMPLE_CORR_NAME, err);  exitOnError(*err,stderr);
   set_char_if_null(names + ID_GAL_SHEAR_RESAMPLE_XCORR, GAL_SHEAR_RESAMPLE_XCORR_NAME, err);  exitOnError(*err,stderr);
   set_char_if_null(names + ID_SHEAR_ALL_RESAMPLE_CORR_XIP, SHEAR_ALL_RESAMPLE_CORR_NAME_XIP, err);  exitOnError(*err,stderr);
   set_char_if_null(names + ID_SHEAR_ALL_RESAMPLE_CORR_XIM, SHEAR_ALL_RESAMPLE_CORR_NAME_XIM, err);  exitOnError(*err,stderr);
   set_char_if_null(names + ID_GAL_SHEAR_ALL_RESAMPLE_CORR_GT, ID_GAL_SHEAR_ALL_RESAMPLE_CORR_NAME_GT, err);  exitOnError(*err,stderr);

   /* Add suffix to names if required */
   if (names[ID_SUFF] != NULL) {
      for (b=0; b<NNAMES-1; b++) {
         if (strcmp(names[b], "")!=0) strcat(names[b], names[ID_SUFF]);
      }
   }

   if (version_and_exit) {
      printf("%g\n", VERSION);
      return 0;
   }

   LOG = fileopen(names[ID_LOG], "w");
   fprintf(LOG, "athena (%s) (v%g) compiled on %s %s\n", __FILE__, VERSION, __DATE__, __TIME__);
   fprintf(LOG, "Command line:");
   fprintf(LOG, " %s", str_argv);
   fprintf(LOG, "\n");
   t_start = start_time(LOG);
   if (!quiet) fprintf(stderr, "athena v%g started\n", VERSION);


   /* Read config file */
   if (!quiet) fprintf(stderr, "Reading config file %s...\n", cname);
   read_config_file(&config, cname, &ndim, &ndimbin, err);
   exitOnError(*err, stderr);
   out_config(LOG, cname, config, ndimbin);
   fprintf(LOG, "ndim    = %d\nndimbin = %d\n", ndim, ndimbin);


   /* Add file type suffix to names */
   for (b=0; b<NNAMES-1; b++) {
      if (config.FORMAT == f_fits) {
         if (strcmp(names[b], "") != 0) { strcat(names[b], SUFFIX_FITS); }
      } else {
         if (strcmp(names[b], "") != 0) strcat(names[b], SUFFIX_ASCII);
      }
   }


   /* Set config entries from command line */
   config.do_wtheta = do_wtheta;

   /* MKDEBUG: Bug with weighted logarithmic scales. *
    * This error message in athena v1.55.            */
   if (config.do_wtheta == 1 && config.I_BINTYPE != LIN) {
      fprintf(stderr, "Weighted angular scales (option '-w') and logarithmic bins (config_tree:BINTYPE=LOG) not yet implemented\n");
      return 2;
   }

   /* Initialize angular bins */
   bin = init_bin(config.thmin, config.thmax, config.nth, ndimbin, config.I_BINTYPE, err);
   exitOnError(*err, stderr);

   if (show_theta_and_exit) {
      double th,th1,th2;
      /* TODO: dims>1 */
      printf("# index theta[rad] theta[%s] bin_min[%s] bin_max[%s]\n",
	     scoord_t(config.COORD_OUTPUT),scoord_t(config.COORD_OUTPUT),scoord_t(config.COORD_OUTPUT));
      for (b=0; b<bin->N[0]; b++) {
         th = bin_index_to_scale(bin, NULL, b);
         th1 = bin_index_to_scale_low(bin, NULL, b);
         th2 = bin_index_to_scale_high(bin, NULL, b);
         printf("%d %g  ", b, th);
         rad_to_coordinate(&th, config.COORD_OUTPUT, err);
         rad_to_coordinate(&th1, config.COORD_OUTPUT, err);
         rad_to_coordinate(&th2, config.COORD_OUTPUT, err);
         exitOnError(*err, stderr);
         printf("%g %13.10f %13.10f \n", th, th1 ,th2); 
      }
      return 0;
   }

   if (config.ERROR == bootstrap) {
      random_init("init");
   }

   /* Read galaxy catalogue(s) */
   read_all_gal(config.GALCAT1, config.GALCAT2, &ngal1, &ngal2, config.NCAT, config.RADEC,
		config.FORMAT, config.COORD_INPUT, &gal1, &gal2, ndim, (const char**)config.COL_NAMES, config.NCOL, quiet, err);
   exitOnError(*err, stderr);

   check_input_format(config.WCORR, config.FORMAT, (const char**)config.COL_NAMES, config.NCOL, err);
   exitOnError(*err, stderr);

   /* Some debug output */
   fprintf(LOG, "First and last galaxy (for debug info):\n");
   fprintf(LOG, "Galaxy #%d (%s), cat #1:\n", 0, "first"); out_galaxy(gal1, 0, ndim, LOG);
   fprintf(LOG, "Galaxy #%u (%s), cat #1:\n", ngal1-1, "last"); out_galaxy(gal1, ngal1-1, ndim, LOG);
   if (config.NCAT==2) {
      fprintf(LOG, "Galaxy #%d (%s), cat #2:\n", 0, "first"); out_galaxy(gal2, 0, ndim, LOG);
      fprintf(LOG, "Galaxy #%u (%s), cat #2:\n", ngal2-1, "last"); out_galaxy(gal2, ngal2-1, ndim, LOG);
   }


   /* === Initialise stuff === */

   init_jackboot(&config, gal1, gal2, ngal1, ngal2, quiet, err);
   exitOnError(*err, stderr);
   if (config.ERROR == jackknife) {
      fprintf(LOG, "After resample initialisation: NRESAMPLE = (%d, %d)\n", config.NRESAMPLE[0], config.NRESAMPLE[1]);
      out_cat_jack(gal1, ngal1, config.NRESAMPLE[0], config.GALCAT1, config.COORD_INPUT, err);
      if (config.NCAT == 2) out_cat_jack(gal2, ngal2, config.NRESAMPLE[1], config.GALCAT2, config.COORD_INPUT, err);
   }

   init_tree(&root1, &root2, gal1, gal2, ngal1, ngal2, ndim, twins_merge, &config, quiet, LOG);

   xi = init_xi(bin, config.WCORR, IMAX(config.NRESAMPLE[0], config.NRESAMPLE[1]), config.do_wtheta, err);
   exitOnError(*err, stderr);


#ifdef ULLONG_MAX
   fnpairtot = (long double)root1->ngal*(long double)root2->ngal;
   if (fnpairtot>ULLONG_MAX) {
      print_err_log_func(LOG, 0, "*** Warning: nb of pairs larger than ULLONG_MAX\n");
   }
#endif
   npairtot = (unsigned long long)root1->ngal*(unsigned long long)root2->ngal;

   /* === Calculate correlations === */


   print_err_log(LOG, quiet, "Calculating twopcf, npairtot=%llu (%.2e)\n", npairtot, (double)npairtot);

   quiet_gl = quiet;
   twopcf((const galaxy_data*)gal1, (const galaxy_data*)gal2, ngal1, ngal2, ndim, root1, root2, bin, &config, xi);

   //if (config.WCORR&wn) fclose(ZGAP_TMP);


   weigh(bin, xi, config.WCORR, IMAX(config.NRESAMPLE[0], config.NRESAMPLE[1]), err);
   exitOnError(*err, stderr);


   /* MKDEBUG TODO: Use real number of pairs for each bootstrap sample? */
   if (config.WCORR&wn) {
      for (b=0; b<bin_N(bin); b++){
	 for (c=0; c<IMAX(config.NRESAMPLE[0], config.NRESAMPLE[1]); c++) {
	    xi->npair_resample[b][c] /= (double)xi->npair[b];
	 }
      }
   }

   sig_eps2 = sig_eps4 = 0.0;
   if (config.WCORR & gg || config.WCORR & gn) {
       sig_eps2 = sigma_epsilon_sqr(gal1, ngal1);
   }
   if (config.WCORR & gg) {
      sig_eps4 = sig_eps2 * sigma_epsilon_sqr(gal2, ngal2);
      fprintf(LOG, "cat1*cat2: sig_eps4 = %.3e, sig_eps2 = %.3e, sig_eps = %.3f\n", sig_eps4, sqrt(sig_eps4), pow(sig_eps4, 0.25));
   }
   if (config.WCORR & gn) {
      fprintf(LOG, "cat1: sig_eps2 = %.3e, sig_eps = %.3e\n", sig_eps2, sqrt(sig_eps2));
   }


   /* === Write correlations to files === */

   /* Mean weight (used for wn) */ 
   if (config.WCORR==wn) {
      if (config.WCORR_SUBTYPE==wn_all) {
         /* Mean weight over whole bg catalogue */
         wbar = root1->weight/ngal1;
      } else {
         /* Pairwise */
         for (b=0,wbar=0.0; b<bin->N[0]; b++) {
            wbar += xi->wn[b];
            print_err_log(LOG, quiet, "mean magnitude of bin %d = %g\n", b, xi->wn[b]);
         }
         wbar /= (double)bin->N[0];
         print_err_log(LOG, quiet, "mean magnitude averaged over all bins = %g\n", wbar);
      }
   } else {
      wbar = 0.0;
   }

   if (config.FORMAT == f_fits) {

#ifdef _WITH_FITS
      write_all_correlations_fits(xi, (const char**)names, bin, &config, sig_eps2, sig_eps4, wbar, ngal1, ngal2, quiet, err);
      exitOnError(*err, stderr);
#else
      *err = addError(gc_fits, "Re-compile 'athena' with fits support", *err, __LINE__);
      exitOnError(*err, stderr);
#endif
   } else {

      write_all_correlations(xi, (const char**)names, bin, &config, sig_eps2, sig_eps4, wbar, ngal1, ngal2, err);
      exitOnError(*err, stderr);

   }


   fprintf(LOG, "\n"); end_time(t_start, LOG);

   fileclose(LOG);


   /* === Clean up === */

   free_xi_data(xi, bin_N(bin));
   free(bin);
   free_galaxy_data(gal1);
   free(root1);
   if (config.NCAT==2) {
      free_galaxy_data(gal2);
      free(root2);
   }

   if (!quiet) fprintf(stderr, "athena finished successfully\n");

   return 0;
}
