/* ============================================================ *
 * kdtree.c							*
 * Martin Kilbinger 2006 - 201.					*
 * ============================================================ */


#include "kdtree.h"


/* ============================================================ *
 * Partition algorithm. Used in grow_tree.			*
 * ============================================================ */
uint partition_data(void *data, uint start, uint end, int ndim, int split_dim,
		    double middle, swap_func *swap, get_pos_func *get_pos)
{
   int s, e;

   s = start;
   e = end - 1;

   /* MKDEBUG: CHECK. Is the following condition ever met? It returns 0/1 ! */
   if (e==s) return get_pos(data, s, split_dim)>=middle;   /* Only one element */

   while (1) {

      while (get_pos(data, s, split_dim)<middle) {
	 s++;
	 if (s==end) return s;     /* All galaxy positions >= middle */
      }
      while (get_pos(data, e, split_dim)>=middle && e>s) {
	 e--;
      }

      if (e<=s) return s;

      swap(data, s, e, ndim);
      s++; e--;

      if (e<s) return s;

   }

}

/* ============================================================ *
 * Returns 1 if aperture with radius th around pos is not       *
 * completely in node nd. Only valid in 2d			*
 * ============================================================ */
int proximity(double pos[2], const node *nd, double th)
{
   return (pos[0]-nd->min[0]<th || nd->max[0]-pos[0]<th ||
           pos[1]-nd->min[1]<th || nd->max[1]-pos[1]<th);
}

node_link *list_init()
{
   node_link *link = malloc(sizeof(node_link)); assert(link);
   link->nd   = 0;
   link->next = 0;
   return link;
}

void list_append(node_link *link, node *nd)
{
   node_link *new = malloc(sizeof(node_link)); assert(new);
   link->next = new;
   new->next = 0;
   new->nd = nd;
}

void list_unchain(node_link *link)
{
   node_link *tmp, *tmp2;

   tmp = link;
   while (tmp) {
      tmp2 = tmp->next;
      // ??
      if (tmp->nd->ngal_resample!=NULL) free(tmp->nd->ngal_resample);
      if (tmp->nd->weight_resample!=NULL) free(tmp->nd->weight_resample);
      //if (tmp->nd->well_resample!=NULL) free(tmp->nd->well_resample); // MKDEBUG complex
      free(tmp->nd);
      free(tmp);
      tmp = tmp2;
   }
}

/* ============================================================ *
 * Returns 1 if the circle with center pos and radius th inter- *
 * sects the node nd. Only valid in 2d.				*
 * MKDEBUG TODO: Has to be revised for RADEC.			*
 * ============================================================ */
int check_intersect(const double pos[2], double th, double thsqr, const node *nd)
{
   double min[2], max[2];
   int k;

   /* Translate coordinates to shift circle centre to origin  */
   for (k=0; k<2; k++) {
      min[k] = nd->min[k]-pos[k];
      max[k] = nd->max[k]-pos[k];
   }

   if (max[0]<0)              /* R to left of circle center   */
     if (max[1]<0)            /* R in lower left corner       */
       return ((Dsqr(max[0]) + Dsqr(max[1]))<thsqr);
     else if (min[1]>0)       /* R in upper left corner       */
       return ((Dsqr(max[0]) + Dsqr(min[1]))<thsqr);
     else                     /* R due West of circle         */
       return (-max[0]<th);
   else if (min[0]>0)         /* R to right of circle center  */
     if (max[1]<0)            /* R in lower right corner      */
       return ((Dsqr(min[0]) + Dsqr(max[1]))<thsqr);
     else if (min[1]>0)       /* R in upper right corner      */
       return ((Dsqr(min[0]) + Dsqr(min[1]))<thsqr);
     else                     /* R due East of circle         */
       return (min[0]<th);
   else			      /* R on y-axis                  */
     if (max[1]<0)            /* R due South of circle        */
       return (-max[1]<th);
     else if (min[1]>0)       /* R due North of circle        */
       return (min[1]<th);
     else                     /* R contains circle center     */
       return 1;
} 

/* ============================================================ *
 * Returns the dimension of largest node extension. Used in     *
 * grow_tree to define coordinate to be split into sub-nodes.	*
 * ============================================================ */
int get_split_dim(const double *min, const double *max, int ndim)
{
   int split_dim = -1;
   double extension = -1.0;
   int k;

   for (k=0; k<ndim; k++) {
      if (extension < max[k] - min[k]) {
	 extension = max[k] - min[k];
	 split_dim = k;
      }
   }

   return split_dim;
}

/* ============================================================ *
 * Handles the galaxies \in [start; end[, creates a new node    *
 * comprising of these galaxies and takes further action if     *
 * necessary, i.e. splits the node by recursive calls.		*
 * ============================================================ */
node *grow_tree(void *data, uint start, uint end, int ndim, int Nboots, int twins_merge, int quiet,
		node_fill_func *node_fill, swap_func *swap, get_pos_func *get_pos, out_func *out)
{
   node *new;
   int split_dim, mid, i;
   double middle;

   new        = malloc(sizeof(node)); assert(new);
   new->start = start;
   new->end   = end;
   new->ngal  = end-start;

   /* This assertion fails if there are multiple identical positions in the catalog *
    * or if the catalogue contains non-ASCII characters.			    */
   if (end<=start) {
      fprintf(stderr, "Error: end<=start (%ud<=%ud)\n", end, start);
      out(data, start, ndim, stderr);
      assert(0);
   }

   /* Fill information into new node (from data) */
   node_fill(new, data, start, end, ndim, Nboots);

   if (new->ngal!=1) {    /* Not a leaf node */

      /* Use minimum size condition or split down to one galaxy ? */

      // MKDEBUG TODO! non-Cartesian coordinates?
      /* Returns the dimension of largest node extension. *
       * Earlier, only 2d was defined, x = 0, y = 1.      */
      split_dim = get_split_dim(new->min, new->max, ndim);
      middle    = (new->min[split_dim] + new->max[split_dim]) / 2.0;

      mid = partition_data(data, start, end, ndim, split_dim, middle, swap, get_pos);

      if ((start==mid || mid==end) && twins_merge==1) {

	 if (!quiet) {
	    fprintf(stderr, "Merging twins %u ... %u\n", start, end-1);
	    for (i=start; i<end; i++) {
	       fprintf(stderr, "%d ", i);
	       out(data, i, ndim, stderr);
	    }
	 }

	 new->left = new->right = NULL;

      } else {

	 new->left  = grow_tree(data, start, mid, ndim, Nboots, twins_merge, quiet, node_fill, swap, get_pos, out);
	 new->right = grow_tree(data, mid, end, ndim, Nboots, twins_merge, quiet, node_fill, swap, get_pos, out);

      }

   } else {
      new->left = new->right = NULL;
   }

   return new;
}

void free_tree(node *nd)
{
   if (nd->left!=NULL) free_tree(nd->left);
   if (nd->right!=NULL) free_tree(nd->right);

   if (nd->ngal_resample!=NULL) free(nd->ngal_resample);
   if (nd->weight_resample!=NULL) free(nd->weight_resample);
   // MKDEBUG well_resample
   free(nd);
}

/* ============================================================ *
 * Calculates the minimum (mm[0]) and maximum (mm[1]) distance  * 
 * between the two nodes n1 and n2. Only applicable in 		*
 * Cartesian coordinates. Only valid in 2d.			*
 * Only needed for close_pairs. Todo: Maximum without corner    *
 * variable, which was removed from node (09/2011).		*
 * ============================================================ */
void minmax(node *n1, node *n2, double *mm)
{
   /* Minimum distance */
   if (n1->max[0]<n2->min[0])            /* n1 left of n2 */
     if (n1->max[1]<n2->min[1])          /* n1 lower left of n2 */
       mm[0] = sqrt(Dsqr(n1->max[0]-n2->min[0]) + Dsqr(n1->max[1]-n2->min[1]));
     else if (n1->min[1]>n2->max[1])     /* n1 upper left of n2 */
       mm[0] = sqrt(Dsqr(n1->max[0]-n2->min[0]) + Dsqr(n1->min[1]-n2->max[1]));
     else			         /* n1 due left of n2 */
       mm[0] = n2->min[0]-n1->max[0];
   else if (n1->min[0]>n2->max[0])       /* n1 right of n2 */
     if (n1->max[1]<n2->min[1])          /* n1 lower right of n2 */
       mm[0] = sqrt(Dsqr(n1->min[0]-n2->max[0]) + Dsqr(n1->max[1]-n2->min[1]));
     else if (n1->min[1]>n2->max[1])     /* n1 upper right of n2 */
       mm[0] = sqrt(Dsqr(n1->min[0]-n2->max[0]) + Dsqr(n1->min[1]-n2->max[1]));
     else				 /* n1 due right of n2 */
       mm[0] = n1->min[0]-n2->max[0];
   else					 /* n1 and n2 share same y coord */
     if (n1->max[1]<n2->min[1])          /* n1 below n2 */
       mm[0] = n2->min[1]-n1->max[1];
     else if (n1->min[1]>n2->max[1])     /* n1 above n2 */
       mm[0] = n1->min[1]-n2->max[1];
     else                                /* n1 and n2 intersect */
       mm[0] = 0.0;

   /* Maximum distance */
   assert(0);

   mm[1] = sqrt(mm[1]);

}

/* ============================================================ *
 * Returns the open angle (in the small-angle approximation,    *
 * tan(x)=x), seen from the base node barycenter to the tangent *
 * node. Works also for RADEC except that the result is		*
 * probably not an "angle" any more.				*
 * ============================================================ */
#define EPSILON  1.0e-6
#define EPSILON2 1.0e-15
double open_angle(const node *nbase, const node *ntan, double *distance, int ndim, int radec)
{
   if (*distance<0) {

      if (radec==0) {

	 *distance = sqrt( absnsqr(ntan->barycenter, nbase->barycenter, ndim) );

      } else {

	 *distance = (ntan->cosbarycenter[0]*nbase->cosbarycenter[0] +
		      ntan->sinbarycenter[0]*nbase->sinbarycenter[0])
	   * ntan->cosbarycenter[1]*nbase->cosbarycenter[1]
	   + ntan->sinbarycenter[1]*nbase->sinbarycenter[1];
	 if (fabs(*distance)>1.0+EPSILON) {
	    fprintf(stderr, "Trying to take acos of %.10f!\n", *distance);
	    assert(0);
	 } else if (*distance>1.0) {
	    *distance = 0.0;   /* acos(1) */
	 } else if (*distance<-1.0) {
	    *distance = pi;    /* acos(-1) */
	 } else {
	    *distance = acos(*distance);
	 }
	 //assert(fabs(*distance)<1.0);
	 //*distance = acos(*distance);

      }

   }
   //return atan2(ntan->radius, *distance);
   /* TODO: ok for radec?? */

   if (*distance < EPSILON2) {
      return 1.0e30;
   } else {
      return ntan->radius/(*distance);
   }
}
#undef EPSILON
#undef EPSILON2

/* ============================================================ *
 * Determines close galaxy pairs. Only valid in 2d.		*
 * ============================================================ */
void pair_list(node *n1, node *n2, double dmin, double dmax, 
	       uint *gal_num, uint *cur, uint ngal)
{
   real mm[2];
   uint i;
   node *nmax, *nmin;

   minmax(n1, n2, mm);

   /* nodes further apart than largest requested dist. */
   if (mm[0]>dmax) return;

   /* nodes nearer than smallest requested distance */
   if (mm[1]<dmin) return;

   if (n1->end-n1->start>n2->end-n2->start) {
      nmax = n1;
      nmin = n2;
   } else {
      nmax = n2;
      nmin = n1;
   }

   /* all galaxies within both nodes match request */
   if (dmin<=mm[0] && mm[1]<=dmax) {

      /* add close pairs if more than one galaxy in node */
      if (nmax->end-nmax->start>1)
	for (i=nmax->start; i<nmax->end; i++) gal_num[(*cur)++] = i;

      if (nmin->end-nmin->start>1 &&
	  /* nmin not descendant of nmax */
	  (nmin->start<=nmax->end || nmin->end<=nmax->start))
	for (i=nmin->start; i<nmin->end; i++) gal_num[(*cur)++] = i;

      assert(*cur<ngal*ngal);
      return;

   } else {

      assert(nmax->end-nmax->start>1);
      if (nmax->end-nmax->start==1) {
	 minmax(n1, n2, mm);
      }

      /* dist = abs22(nmax->corner, nmin->corner); */
      /* if (dmin<=dist && dist<=dmax) gal+num[]; */
      /* } */
      assert(0);

      pair_list(nmax->left, nmin, dmin, dmax, gal_num, cur, ngal);
      pair_list(nmax->right, nmin, dmin, dmax, gal_num, cur, ngal);

   }

   return;
}

double numberdensity(const node *nd, int ndim)
{
   return nd->ngal/volume(nd, ndim);
}

double volume(const node *nd, int ndim)
{
   double V;
   int k;

   for (k=0,V=1.0; k<ndim; k++) {
      V *= nd->max[k] - nd->min[k];
   }

   return V;
}
