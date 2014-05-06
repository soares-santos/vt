#include "voronoi.h"

void usage(void) {
  fprintf(stderr, "\nUsage:\n"                                              );
  fprintf(stderr, " vt [-avxPhTR] [-A amplitude] [-g power] [-r min_cl] "   );
  fprintf(stderr, "[-s scl]\n    [-C galaxycatalogcols] "                   );
  fprintf(stderr, "[-B maskfilecols] [-F framecoords] "                     );
  fprintf(stderr, "[-b maskfile]\n    [-i idprefix] galaxy_catalog\n"       );
  fprintf(stderr, "\n Detect galaxy clusters from the galaxy_catalog."      );
  fprintf(stderr, " Assuming the angular 2pt\n correlation function"        );
  fprintf(stderr, " of the background galaxy distribution to be\n"          );
  fprintf(stderr, " w(theta) = A * theta^(1-g).\n"                          );
  fprintf(stderr, "\nOptions:\n"                                            );
  fprintf(stderr, " -a: input files are ascii files, not fits        \n"    );
  fprintf(stderr, " -v: verbose                                      \n"    );
  fprintf(stderr, " -x: input pix coordinates (x,y), not (RA,DEC)    \n"    );
  fprintf(stderr, " -P: Poisson background                           \n"    ); 
  fprintf(stderr, " -h: print this help and exit                     \n"    );
  fprintf(stderr, " -T: do the Tessellation only, do not find clusters\n"   );
  fprintf(stderr, " -R: raw finder: dump candidates list and exit     \n"   );
  fprintf(stderr, " -A: w(theta) amplitude parameter       default: 0.005\n");
  fprintf(stderr, " -g: w(theta) power law parameter       default: -1.70\n");
  fprintf(stderr, " -r: min conf. level above background   default: 0.999\n");
  fprintf(stderr, " -s: fixed scl, not Ebeling's method      \n"            );
  fprintf(stderr, " -C: colnames in galaxy catalog         "                );
  fprintf(stderr, "default: id ra dec z\n"                                  );
  fprintf(stderr, " -B: colnames in mask file              "                );
  fprintf(stderr, "default: ramin ramax decmin decmax\n"                    );
  fprintf(stderr, " -F: frame limits for valid data        "                );
  fprintf(stderr, "default: 0 360 -90 90\n"                                 );
  fprintf(stderr, " -b: bad areas file name                \n"              );
  fprintf(stderr, " -i: zeropt for cluster id numbers      default: 0\n"    );
  fprintf(stderr, "\nNotes:\n"                                              );
  fprintf(stderr, " -- w(theta) parameters must obey 0 <= A  "              );
  fprintf(stderr, "(A = 0 -> Poisson) and\n    1 <= g < 2.\n"               );
  fprintf(stderr, " -- A and g can be omitted in the Poisson case.\n"       );
  fprintf(stderr, " -- Ebeling's  method is default.\n"                     );
  fprintf(stderr, " -- galaxycatalog must contain the cols: id ra dec z\n"  );
  fprintf(stderr, " -- maskfile must have the cols: "                       );
  fprintf(stderr, "ramin ramax decmin decmax\n"                             );
  fprintf(stderr, " -- If the '-a' option is set, the input files must "    );
  fprintf(stderr, "contain only the columns\n    listed (plus optional "    );
  fprintf(stderr, "header/comment lines) and in the exact order above.\n"   );
  fprintf(stderr, " -- In case of fits files, the order is not important "  );
  fprintf(stderr, "and there may be extra\n    "                            );
  fprintf(stderr, "cols. They will be selected by column name (case "       );
  fprintf(stderr, "insensitive). "                                          );
  fprintf(stderr, "Column\n    names can be changed via commandline opts.\n"); 
  fprintf(stderr, " -- Output .VTclusters.fit(.cat) will have the cols:\n  ");
  fprintf(stderr, "  id ra dec z_median z_mean z_rms contrast sig nvt nbg\n");
  fprintf(stderr, "    ncorr size a b pa area min_dens ra_c dec_c z_c"      );
  fprintf(stderr, " x_0 y_0\n"                                              );
  fprintf(stderr, " -- Output .VTgalaxies.fit(.cat) will have the cols:\n"  );
  fprintf(stderr, "    id ra dec z host_id local_dens central x_0 y_0    \n");
  fprintf(stderr, " -- A list of cell areas is saved in the .areas file. \n");
  fprintf(stderr, " -- The empirical and theoretical cell density "         );
  fprintf(stderr, "distributions are saved in\n    the "                    );
  fprintf(stderr, ".distributions file.\n"                                  );
  fprintf(stderr, " -- The tessellation diagrams are in files .voronoi "    );
  fprintf(stderr, "and .delaunay\n"                                         );
  fprintf(stderr, " -- The .raw file holds a preliminary candidates "       );
  fprintf(stderr, "list. Cols are:\n    number min_dens ncorr nvor nback\n" );
  fprintf(stderr, " -- A simple memebers list containing only glx_id and "  );
  fprintf(stderr, "host_id is saved into a\n    VTmembers.list file\n "     );
  fprintf(stderr, "\n"                                                      );
  exit(8);
} /* usage() */

int main(int argc, char *argv[]) {
  register int i, j;
  int **Prox, **Edge, *nprox, *p_prox, *p_edge, *mask2;
  int n_tmp,ia,ib,n1,file_output=0,i_gr;
  REAL  *f, tmp;
  Connect *connection_1, *p_conn;

  if(argc > 1) set_defaults();

  while(--argc > 0 && (*++argv)[0] == '-') while(argv[0][1] != '\0')
    switch(*++argv[0]) { 
    case 'a': Fits=0;    break;
    case 'v': Verbose=1; break;
    case 'x': XYcoord=1; break;
    case 'P': Poisson=1; break;
    case 'h': usage();   break;
    case 'T': NotFind=1; break;
    case 'R': RawFind=1; break;
    case 'A': 
      if(argv[0][1] == '\0') argc--, w_amp=atof((++argv)[0]);
      else w_amp=atof(++argv[0]);
      while(argv[0][1] != '\0') *++argv[0];
      if(w_amp < 0) usage();
      if(w_amp >0.1) w_amp=0.1;
      break;
    case 'g': 
      if(argv[0][1] == '\0') argc--, w_pow=atof((++argv)[0]);
      else w_pow=atof(++argv[0]);
      while(argv[0][1] != '\0') *++argv[0];
      if(w_pow < 1 || w_pow > 2) usage(); 
      break;
    case 'r': 
      if(argv[0][1] == '\0') argc--, reject_lev=atof((++argv)[0]);
      else reject_lev=atof(++argv[0]);
      while(argv[0][1] != '\0') *++argv[0];
      if(reject_lev < 0 || conf_lev > 1) usage(); 
      break;
    case 's':
      Ebeling=0;
      if(argv[0][1] == '\0') argc--, conf_lev=atof((++argv)[0]);
      else conf_lev=atof(++argv[0]);
      while(argv[0][1] != '\0') *++argv[0];
      if(conf_lev < 0 || conf_lev > 1) usage(); 
      break;
    case 'C':
      for(i=0;i<4;++i){
	if(argv[0][1] == '\0') argc--, ttype[i]=(++argv)[0];
	else ttype[i]=++argv[0];
	while(argv[0][1] != '\0') *++argv[0];
	if(argc<2 || ttype[i][0] == '-' ) usage();
      }
      break;
    case 'B':
      for(i=0;i<4;++i){
	if(argv[0][1] == '\0') argc--, ttypeH[i]=(++argv)[0];
	else ttypeH[i]=++argv[0];
	while(argv[0][1] != '\0') *++argv[0];
	if(argc<2 || ttypeH[i][0] == '-' ) usage();
      }
      break;
    case 'F':
      for(i=0;i<4;++i){
	if(argv[0][1] == '\0') argc--, frame[i]=atof((++argv)[0]);
	else frame[i]=atof(++argv[0]);
	while(argv[0][1] != '\0') *++argv[0];
	if(argc<2 || argv[0][0] == '-') usage();
      }
      break;
    case 'b': 
      Holes=1;
      if(argv[0][1] == '\0') argc--, file_holes=(++argv)[0];
      else file_holes= ++argv[0];
      while(argv[0][1] != '\0') *++argv[0];
      break;
    case 'i':
      if(argv[0][1] == '\0') argc--, FieldID=atoll((++argv)[0]);
      else FieldID=atoll(++argv[0]);
      while(argv[0][1] != '\0') *++argv[0];
      if(FieldID < 0) FieldID=0;
      break;
    default: usage();
    } /* switch loop */ 

  if (argc != 1) usage(); 

  if(Verbose) printf("\nVT :: 2D :: INITIALIZATION\n");

  file_name=argv[0];
  print_init();
  set_TF_params();
  read_data();
  if(Holes) read_holes();
  adjust_frame();

  if(Verbose) printf("\nVT :: 2D :: TESSELLATION\n");

  /* Triangle returns triangulation in `mid' and a voronoi diagram in `vor'*/
  mid.neighborlist = (int *) NULL;     /* Needed only if -n switch used. */
  mid.edgelist = (int *) NULL;         /* Needed only if -e switch used. */
  vor.pointlist = (REAL *) NULL;       /* Needed only if -v switch used. */
  vor.edgelist = (int *) NULL;
  vor.normlist = (REAL *) NULL;
  /* Triangulate the points.  Switches are chosen to preserve the convex */
  /* hull (c), number everything from zero (z), produce an edge list (e),*/
  /* a Voronoi diagram (v), and a triangle neighbor list (n).            */
  /* B: no boundary markers; N: no mid.output; E: no ele ...  ; Q:Silent */
  triangulate("QBENczevn", &in, &mid, &vor);
  if(Verbose) printf("Voronoi diagram succesfully built.\n");
  free(mid.neighborlist);
  free(vor.normlist);
  /* Prepare to assign neighbors and voronoi edges to each point */
  ARRAY2(int,Prox,ndata,20);
  ARRAY2(int,Edge,ndata,20);
  nprox=(int *) calloc(ndata, sizeof(int)); /* Initialize to zeros */
  for(i=0;i<mid.numberofedges;i++) {
    ia=mid.edgelist[2*i]; ib=mid.edgelist[2*i+1];
    Prox[ia][nprox[ia]]=ib; Prox[ib][nprox[ib]]=ia;
    Edge[ia][nprox[ia]]=Edge[ib][nprox[ib]]=i;
    nprox[ia]++; nprox[ib]++;
  }
  for(i=0;i<ndata;i++) if(nprox[i]==0)
    printf("Warning: point %d has no neighbors; probably double point\n", i);
  /* Assigning neighbors and voronoi edges to each point */
  n_tmp=ndata+2*mid.numberofedges+1;
  ARRAY1(int,prox,n_tmp);       /* the neighbors */
  ARRAY1(int,edge,n_tmp);       /* the edges of Voronoi polygons */
  prox[0]=edge[0]=ndata+1;
  p_prox=prox+ndata; p_edge=edge+ndata;
  for(i=0; i<ndata; i++) {
    n_tmp=nprox[i];
    prox[i+1]=edge[i+1]=prox[i]+n_tmp;
    for(j=0; j<n_tmp; j++) {
      *(++p_prox)=Prox[i][j];
      *(++p_edge)=Edge[i][j];
    }
  }
  free(nprox); free2(Prox); free2(Edge);
  if(Verbose) printf("Cell edges and neighbors computed.\n");
  /* Compute areas */
  ARRAY1(REAL, area, ndata);
  ARRAY1(REAL, f, ndata);
  ARRAY1(short, voronoi_inside, vor.numberofedges);
  for(i=0;i<vor.numberofedges;i++) voronoi_inside[i]=0;
  for(i=0;i<ndata;i++) {
    area[i]=tmp=poly_area(i); /* this is where the holes area considered! */
    if(tmp>0.) f[ndata_ok++]=1./tmp; /* Select internal points */
    tot_area_covered += tmp;
  }
  WriteAreasFile();
  WriteVoronoiCatalog(); 
  WriteDelaunayCatalog(); 
  area_mean_ok=tot_area_covered/ndata_ok; /* mean area of ALL tessels */
  free(edge); free(HolesCoord);
  ndataused=ndata_ok;
  qsort((REAL *) f, ndata_ok, sizeof(REAL), real_compare);
  Nbg=Background(f,tot_area_covered,ndata_ok);
  area_mean=tot_area_covered/Nbg;  
  background_density = Nbg/tot_area_covered; 
  if(Verbose){ 
    printf("Mean area: %14.8f\n", area_mean);
    printf("Mean density: %14.8f\n", background_density); 
  }
  if(NotFind){
    if(Verbose) printf("Not finding clusters. Ending VT.\n");
    return 0;
  }

  if(Verbose) printf("\nVT :: 2D :: DETECTION\n");

  for(i=0;i<ndata_ok;i++) f[i]/=background_density;
  //qsort((REAL *) f, ndata_ok, sizeof(REAL), real_compare);
  threshold=ComputeThreshold(f); 
  WriteDistributionsFile(f);
  free(f);   
  mask1=(int *) calloc(ndata,sizeof(int));
  mask2=(int *) calloc(ndata,sizeof(int));
  for(i=0;i<ndata;i++) {
    if((tmp=area[i])>0. && area_mean/tmp >= threshold) 
      mask1[i]=mask2[i]=1;   /* mask: high-density points=1, 0 otherwise. */
  }
  if(Verbose) printf("Selected points above the threshold.\n");
  WriteCandidatesFile(prox,mask2);
  free(mask2);
  if(RawFind){
    if(Verbose) printf("Raw cluster finder. Ending VT.\n");
    return 0;
  }
  connection_1=percola(prox, mask1);  /* percolation applying rcl */
  nclusters=n1=CountItem(connection_1);
  free(prox);

  if(Verbose) printf("\nVT :: 2D :: CHARACTERIZATION\n");

  /* Set flags and rewrite mask (0=background,1=cluster) */ 
  central=(int *) calloc(ndata,sizeof(int));
  edge1=(int *) calloc(ndata,sizeof(int));
  for(i=0;i<ndata;i++) {
    mask1[i]=0; 
    central[i]=0; 
    edge1[i]=0;
    if(in.pointlist[2*i]<frame[0] || 
       in.pointlist[2*i]>frame[1])     edge1[i]=1;
    if(in.pointlist[2*i+1]<frame[2] || 
       in.pointlist[2*i+1]>frame[3])   edge1[i]=1; 
  }
  p_conn=connection_1; 
  while(p_conn != NULL) { 
    for(i=0;i<p_conn->n;i++) mask1[p_conn->list[i]]=1; 
    p_conn=p_conn->link; 
  } 

  ARRAY1(Footprint,footprint,n1);
  SetFootprints(connection_1);
  host_id_list_1=(long long *) calloc(ndata,sizeof(long long));
  AssignMembers(connection_1);
  ComputeRedshifts(connection_1);

  if(Verbose) printf("\nVT :: 2D :: OUTPUT\n");

  WriteClusterCatalog(); 
  WriteGalaxyCatalog();

  free(voronoi_inside); 
  free(vor.edgelist); 
  free(vor.pointlist); 
  free(mid.edgelist); 
  free(in.pointlist); 
  free(mask1);
  free(footprint);
  free(connection_1); 
  free(central);
  free(host_id_list_1);
  //  free(ttype);
  //  free(ttypeH);

  return 0;
}

