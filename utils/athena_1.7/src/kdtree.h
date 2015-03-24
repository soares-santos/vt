/* ============================================================ *
 * kdtree.h							*
 * Martin Kilbinger, Melody Wolk				*
 * 01/2006							*
 * 09/2011 Add w_p(r_p), k dimensions				*
 * ============================================================ */


#ifndef __KDTREE_H
#define __KDTREE_H

#include <sys/times.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <sys/types.h>

#include "mathstuff.h"
#include "nrcomplex.h"

typedef unsigned int uint;

/* Dimensionality of data. The 'k' in kd-tree */
#define KTREE_MAX 3

/* Dimensionality of bins. */
#define NDIMBIN_MAX 2


typedef struct dummy_node {

  double min[KTREE_MAX], max[KTREE_MAX];     /* Rectangle corners          */

  dcomplex well;                    /* Summed weighted ellipticity         */
  double weight;                    /* Mean weight		           */
  double barycenter[KTREE_MAX];     /* Mean weighted position              */
  double cosbarycenter[KTREE_MAX],
    sinbarycenter[KTREE_MAX];       /* Used for spherical coordinates only */
  double radius;                    /* Distance from bc to farest galaxy   */
  double ellsqr;                    /* Ellipticity modulus square          */
   //double zmin, zmax;                /* Minimum and maximum redshift, used for WCORR=wn */

  uint ngal, start, end;

  double *ngal_resample, *weight_resample;
  dcomplex *well_resample;

  struct dummy_node *left, *right;
} node;

/* Link in a chain of nodes */
typedef struct dummy_link {

  struct dummy_link *next;
  node *nd;

} node_link;

typedef void fill_function(node *nd, const void *tmp_data, uint start, uint end, int ndim, int Nboots);


/* ============================================================ *
 * The following three functions have to be provided for a      *
 * given type of extra information in each node. Void-pointer	*
 * are used for general types.					*
 * ============================================================ */
typedef void node_fill_func(node *, const void*, uint, uint, int, int);
typedef void swap_func(void *, uint, uint, int);
typedef double get_pos_func(void *, uint, int);
typedef void out_func(const void *, uint, int, FILE *F);


uint partition_data(void *data, uint start, uint end, int ndim, int split_dim, 
		    double middle, swap_func *swap, get_pos_func *get_pos);
int check_intersect(const double pos[2], double th, double thsqr, const node *nd);
int proximity(double pos[2], const node *nd, double th);
node_link *list_init();
void list_append(node_link *link, node *nd);
void list_unchain(node_link *link);
int get_split_dim(const double *min, const double *max, int ndim);
node *grow_tree(void *data, uint start, uint end, int ndim, int Nboots, int twins_merge, int quiet,
		node_fill_func *node_fill, swap_func *swap, get_pos_func *get_pos, out_func *out);
void free_tree(node *nd);
double open_angle(const node *nbase, const node *ntan, double *distance, int ndim, int radec);
void minmax(node *n1, node *n2, double *mm);
void pair_list(node *n1, node *n2, double dmin, double dmax, 
	       uint *gal_num, uint *cur, uint ngal);
double numberdensity(const node *nd, int ndim);
double volume(const node *nd, int ndim);
#endif
