#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <getopt.h>
#include <gsl/gsl_rng.h>

#include "mathstuff.h"
#include "nrcomplex.h"
#include "config.h"
#include "errorlist.h"
#include "2pcf.h"
#include "kdtree.h"
#include "gal_cat.h"
#include "io.h"
#include "stdnames.h"



void add_if_filled(int **number_good, int **number_mask, int i, int j, int N[2], int *n_good, int *n_mask)
{
   if (i>0 && j>0 && i<N[0] && j<N[1]) {
      if (number_good[i][j] > number_mask[i][j]) {
	 (*n_good) ++;
      }
      if (number_good[i][j] < number_mask[i][j]) {
	 (*n_mask) ++;
      }
   }
}

void scoord_to_coord(coord_t *coord, const char scoord[], error **err)
{
   int j;
   STRING2ENUM(*coord, scoord, coord_t, scoord_t, j, Ncoord_t, err);
   forwardError(*err, __LINE__,);
}

/* ============================================================ *
 * Returns GSL random number generator. Initialised with seed.  *
 * If seed==-1, initialised with current time.			*
 * From exec_helper.c						*
 * ============================================================ */
gsl_rng *init_random(int seed, FILE *FLOG)
{
   gsl_rng *rng;
   struct tms Time;
   clock_t cl;
   unsigned u;

   if (seed==-1) {
      cl = times(&Time);
      u  = (unsigned)cl;
   } else {
      u  = (unsigned)seed;
   }

   rng = gsl_rng_alloc(gsl_rng_default);
   gsl_rng_set(rng, u);

   if (FLOG) {
      fprintf(FLOG, "Seed for random generator = %u\n", u);
   }

   return rng;
}


/* ============================================================ *
 * Returns sum of weights in box [x; x+dx]^2.		        *
 * ============================================================ */
void count_in_box(const double x[2], const double dx[2], const node *nd, int *number_good, int *number_mask, const galaxy_data *gal)
{
   uint m;

   /* Node completely in box? */
   if (nd->min[0] >= x[0] && nd->max[0] <= x[0] + dx[0] &&
       nd->min[1] >= x[1] && nd->max[1] <= x[1] + dx[1])
   {
      for (m=nd->start; m<nd->end; m++) {
	 if (gal[m].weight > 0.99) (*number_good) ++;
	 if (gal[m].weight < 0.01) (*number_mask) ++;
      }
      return;
   }

   /* Node completely outside of box? */
   if (nd->max[0] <= x[0] ||            /* Node to the bottom of box */
       nd->max[1] <= x[1] ||            /* Node to the left of box */
       nd->min[0] >= x[0] + dx[0] ||    /* Node to the top of box */
       nd->min[1] >= x[1] + dx[1]) {    /* Node to the right of box */
      return;
   }
   
   /* Call recursively with child nodes */
   count_in_box(x, dx, nd->left, number_good, number_mask, gal);
   count_in_box(x, dx, nd->right, number_good, number_mask, gal);
}

/* ============================================================ *
 * Print some catalogue stats to the log file.			*
 * ============================================================ */
void some_stats(int i, const node *root, int radec, FILE *LOG)
{
   double f;

   f = radec==0 ? 60.0 : 1.0;

   fprintf(LOG, "Galaxy catalogue #%d:\n", i);
   fprintf(LOG, "  radius = %g %s\n", root->radius/f, radec==0 ? "arcmin" : "rad");
   fprintf(LOG, "  area   = %g sq %s\n", area(root)/f/f, radec==0 ? "arcmin" : "rad");
   fprintf(LOG, "  number density = %g/(sq %s)\n", numberdensity(root)*f*f,
	   radec==0 ? "arcmin" : "rad");
}

#define usage_VA(ex, str, ...) char message_VA_usage[TXT_SZ]; sprintf(message_VA_usage, str, __VA_ARGS__); usage(ex, str);
void usage(int ex, const char* str)
{
   if (str!=NULL) {
      fprintf(stderr, "%s\n", str);
   }

   fprintf(stderr, "Usage: bck_eng_random.c OPTIONS\n");
   fprintf(stderr, "OPTIONS:\n");
   fprintf(stderr, "  -i INPUT         Input catalogue INPUT, in 'position' format\n");
   fprintf(stderr, "  -o OUTPUT        Output catalogue OUTPUT base (default 'INPUT.out'.\n");
   fprintf(stderr, "                    The files 'OUTPUT.pix' (pixel mask) and 'OUTPUT.cat' random\n");
   fprintf(stderr, "                    catalogue) will be created\n");
   fprintf(stderr, "  -r               Spherical coordinates (radec=1)\n");
   fprintf(stderr, "  -c COORD         Catalogue coordinates are 'COORD' (default 'arcsec')\n");
   fprintf(stderr, "  -N N             NxN pixels (default 100)\n");
   fprintf(stderr, "  -n R             R random objects (default: number of input objects\n");
   fprintf(stderr, "  -h               This message\n");
   
   if (ex>=0) exit(ex);
}

int main(int argc, char *argv[])
{
   galaxy_data *gal=NULL;
   uint ngal;
   node *root;
   char *cname, *outname;
   char rname[512];
   char *scoord;
   coord_t coord;
   int c, i, j, radec, quiet=0, N[2], Nrand, n, n_good, n_mask;
   uint ir[2];
   double num_avg, r[2];
   int **number_good, **number_mask, **number_second, num_min;
   double x[2], xc[2], dx[2];
   fill_function *node_fill_galaxy[2] = {node_fill_galaxy_xy, node_fill_galaxy_radec};
   gsl_rng *rng;
   extern char *optarg;
   extern int optind, optopt;
   error *myerr = NULL, **err;
   FILE *F;


   err = &myerr;


   fprintf(stderr, "bckw_eng_random (%s) compiled on %s %s\n", __FILE__, __DATE__, __TIME__);


   /* Command line arguments */

   cname   = NULL;
   outname = NULL;
   scoord  = NULL;
   coord   = coord_arcsec;
   radec   = 0;
   Nrand   = 0;
   N[0] = N[1] = 100;
   while ((c = getopt(argc, argv, ":i:o:rc:N:n:h")) != -1) {

      switch (c) {
	 case 'i' :
	    cname = optarg;
	    break;
	 case 'o' :
	    outname = optarg;
	 case 'r' :
	    radec = 1;
	    break;
	 case 'c' :
	    scoord = optarg;
	    scoord_to_coord(&coord, scoord, err);
	    quitOnError(*err, __LINE__, stderr);
	    break;
	 case 'N':
	    N[0] = atoi(optarg);
	    N[1] = N[0];
	    break;
	 case 'n' :
	    Nrand = atoi(optarg);
	    break;
	 case ':' :
	    usage_VA(1, "Argument for -%c option missing\n", optopt);
	 case '?' :
	    usage_VA(2, "Unknown option -%c\n", optopt);
	 case 'h' :
	    usage(3);
      }

   }

   if (cname==NULL) usage(1, NULL);
   if (outname==NULL) {
      outname = malloc_err(sizeof(char)*512, err);
      quitOnError(*err, __LINE__, stderr);
      sprintf(outname, "%s.out", cname);
   }


   /* Read input catalogue and create tree */

   gal = read_gal(cname, &ngal, radec, f_position, coord, 0.0, 0, err);
   quitOnError(*err, __LINE__, stderr);

   if (Nrand==0) Nrand = ngal;

   root = grow_tree(gal, 0, ngal, 0, 0, quiet, node_fill_galaxy[radec], swap_galaxy, get_pos_galaxy, out_galaxy);
   some_stats(1, root, radec, stderr);


   number_mask   = imatrix(0, N[0]-1, 0, N[1]-1);
   number_good   = imatrix(0, N[0]-1, 0, N[1]-1);
   number_second = imatrix(0, N[0]-1, 0, N[1]-1);
   num_min = ngal;
   num_avg = 0.0;

   for (j=0; j<2; j++) {
      dx[j] = (root->max[j] - root->min[j])/((double)N[j]);
   }


   /* Create pixel mask */

   for (i=0,x[0]=root->min[0]; i<N[0]; i++,x[0]+=dx[0]) {

      for (j=0,x[1]=root->min[1]; j<N[1]; j++,x[1]+=dx[1]) {

	 number_good[i][j]  = number_mask[i][j] = 0;
	 count_in_box(x, dx, root, &(number_good[i][j]), &(number_mask[i][j]), gal);
	 n = number_good[i][j] + number_mask[i][j];
	 if (n < num_min) num_min = n;
	 num_avg += n;

      }

   }
   num_avg = num_avg / (double)N[0] / (double)N[1];
   fprintf(stderr, "# Minimum number of galaxies in a box = %d\n", num_min);
   fprintf(stderr, "# Average number of galaxies per box  = %g\n", num_avg);


   /* Fill empty pixels (containing no galaxies), according to neighbouring pixels */
   for (i=0; i<N[0]; i++) {
      for (j=0; j<N[0]; j++) {

	 number_second[i][j] = 0;
	 if (number_good[i][j] == 0 && number_mask[i][j] == 0) {
	    n_good = n_mask = 0;
	    add_if_filled(number_good, number_mask, i-1, j, N, &n_good, &n_mask);
	    add_if_filled(number_good, number_mask, i+1, j, N, &n_good, &n_mask);
	    add_if_filled(number_good, number_mask, i, j-1, N, &n_good, &n_mask);
	    add_if_filled(number_good, number_mask, i, j+1, N, &n_good, &n_mask);
	    if (n_good > n_mask) {
	       //number_good[i][j] = 1;
	       number_second[i][j] = 1;
	    } else {
	       //number_mask[i][j] = 1;
	       number_second[i][j] = -1;
	    }
	 }

      }
   }

   /* Copy data from second run */
   for (i=0,n=0; i<N[0]; i++) {
      for (j=0; j<N[0]; j++) {
	 if (number_good[i][j] == 0 && number_mask[i][j] == 0) {
	    if (number_second[i][j] == 1) number_good[i][j] = 1;
	    else if (number_second[i][j] == -1) number_mask[i][j] = 1;
	    n++;
	 }
      }
   }
   fprintf(stderr, "In second run: n = %d (%.1f%%) empty pixels filled\n", n, 100.0*(double)n/(double)N[0]/(double)N[1]);

   /* Output pixel mask */

   sprintf(rname, "%s.pix", outname);
   F = fileopen(rname, "w");

   for (i=0,x[0]=root->min[0]; i<N[0]; i++,x[0]+=dx[0]) {

      xc[0] = x[0];
      rad_to_coordinate(xc, coord, err);
      quitOnError(*err, __LINE__, stderr);

      for (j=0,x[1]=root->min[1]; j<N[1]; j++,x[1]+=dx[1]) {

	 xc[1] = x[1];
	 rad_to_coordinate(xc+1, coord, err);
	 quitOnError(*err, __LINE__, stderr);

	 fprintf(F, "% 10.3f % 10.3f  %d %d\n", xc[0], xc[1], number_good[i][j], number_mask[i][j]);

      }

   }
   fileclose(F);


   /* Create random catalogue */

   sprintf(rname, "%s.cat", outname);
   F = fileopen(rname, "w");
   rng = init_random(-1, stderr);
   i = 0;
   while (i<Nrand) {
      
      for (j=0; j<2; j++) {
	 /* Draw random integer number from [ 0 ; N[j] ) */
	 ir[j] = gsl_rng_uniform_int(rng, N[j]);
      }

      if (number_good[ir[0]][ir[1]] > number_mask[ir[0]][ir[1]]) {
	 for (j=0; j<2; j++) {
	    /* Draw random double number from [ 0 ; 1 ) */
	    r[j]  = gsl_rng_uniform(rng);
	    x[j]  = root->min[j] + ((double)ir[j] + r[j])*dx[j];
	    xc[j] = x[j];
	    rad_to_coordinate(xc+j, coord, err);
	    quitOnError(*err, __LINE__, stderr);
	 }
	 fprintf(F, "%20.10g %20.10g 1.0\n", xc[0], xc[1]);
	 i++;
      }

   }
   fileclose(F);


   /* Clean up */

   free_imatrix(number_good, 0, N[0]-1, 0, N[1]-1);
   free_imatrix(number_mask, 0, N[0]-1, 0, N[1]-1);
   free_imatrix(number_second, 0, N[0]-1, 0, N[1]-1);
   free_galaxy_data(gal);
   free(root);
   gsl_rng_free(rng);

   return 0;
}
