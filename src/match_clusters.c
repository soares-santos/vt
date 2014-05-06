#include "voronoi.h"

long long *id_cls_1,*id_glx_1,*host_id_1,*id_match_1; // 1st set of input catalogs
long long *id_cls_2,*id_glx_2,*host_id_2,*id_match_2; // 2nd set of input catalogs
long long *id_cls_3,*id_cnd_3,*id_match_3; // output table
long ncls_1,nglx_1,ncls_2,nglx_2,ncls_3;
int *ngals_1,*ngals_2,*ngals_3,*nmatches_1,*nmatches_2,*nmatches_3;
int MoreVerbose,*primary_1,*primary_2,*primary_3;
char *fname[4];
double *mf1,*mf2,*mf3;

void match(void);
void write_match_table(void); 
void alloc_memory(void);
void free_memory(void);
void usage(void) {
  fprintf(stderr, "\nUsage: match_clusters [-vVh] [-C IDcol HOSTIDcol]");
  fprintf(stderr, " cluster_cat1 galaxy_cat1 cluster_cat2 galaxy_cat2\n\n");  exit(8);
} 
  
int main(int argc, char *argv[]) { 
  register int i, j;    
  if(argc > 1) set_defaults();
  while(--argc > 0 && (*++argv)[0] == '-') while(argv[0][1] != '\0')
    switch(*++argv[0]) { 
    case 'v': Verbose=1; break;
    case 'V': MoreVerbose=Verbose=1; break;
    case 'h': usage();   break;
    case 'C':
      for(i=0;i<2;++i){
        if(argv[0][1] == '\0') argc--, ttype[i]=(++argv)[0];
        else ttype[i]=++argv[0];
        while(argv[0][1] != '\0') *++argv[0];
        if(ttype[i][0] == '-' ) usage();
      }
      break;
    default: usage();
    }  
  if (argc != 4) usage(); 
  for(i=0;i<4;i++) fname[i]=argv[i];
  if(Verbose) printf("VT :: match clusters from 1st and 2nd run\n");
  read_data();
  alloc_memory();
  match();
  write_match_table();
  free_memory();
  return 0;
}

void match(void){
  register int i,j,k,l,m,n,p,a,b;
  long long *memb1=NULL,*memb2=NULL,id;
  int *matched=NULL;

  if(Verbose) printf("matching catalogs...\n");

  matched = (int *) realloc(matched,ncls_2*sizeof(int));
  for (k=0;k<ncls_2;k++) matched[k]=0;

  /* for each cluster in catalog 1 ... */
  for(i=0;i<ncls_1;i++){
    /* ... set default values */
    ngals_1[i]     =  0;    m=0;
    nmatches_1[i]  =  0;    n=i;
    primary_1[i]   = -1;    
    id_match_1[i]  = -1;
    mf1[i]         = 0.;
    id_cls_3[n]    = id_cls_1[i]*10+1;
    id_cnd_3[n]    = id_cls_1[i];
    if(MoreVerbose) printf("\ncluster id : %lld   \t",id_cls_1[i]);
    /* ... find members */
    for(j=0;j<nglx_1;j++)
      if(host_id_1[j]==id_cls_1[i]){
	++m;
	ngals_1[i]=m;
	memb1 = (long long *) realloc(memb1,m*sizeof(long long));
	memb1[m-1]=id_glx_1[j];
      }
    ngals_3[n]=ngals_1[i];
    if(MoreVerbose) printf("nmembers : %d\n",ngals_1[i]);
    /* ... try matching to clusters in catalog 2 ... */
    for(k=0;k<ncls_2;k++) if(matched[k]==0) {
      /* ... set default values */
      ngals_2[k]     =  0;    m=0;
      nmatches_2[k]  =  0;    p=ncls_1+k;
      primary_2[k]   = -1;
      id_match_2[k]  = -1;
      mf2[k]         = 0.;
      id_cls_3[p]    = id_cls_2[k]*10+2;
      id_cnd_3[p]    = id_cls_2[k];
      /* ... find members */
      if(MoreVerbose) printf("\n trying match to cluster %lld ...\t",id_cls_2[k]);
      for(j=0;j<nglx_2;j++)
	if(host_id_2[j]==id_cls_2[k]){
	  ++m;
	  ngals_2[k]=m;
	  memb2 = (long long *) realloc(memb2,m*sizeof(long long));
	  memb2[m-1]=id_glx_2[j];
	}
      ngals_3[p]=ngals_2[k];
      if(MoreVerbose) printf("nmembers : %d ...\t",ngals_2[k]);
      /* ... count shared galaxies */     
      a=b=0;
      for(j=0;j<ngals_1[i];j++)
	for(l=0;l<ngals_2[k];l++)
	  if(memb1[j]==memb2[l]) {a++; b++;}  
      mf1[i]=(double)a / (double)ngals_1[i];
      mf2[k]=(double)b / (double)ngals_2[k];
      if(MoreVerbose) printf("match fraction 1 : %f ...\t\t match fraction 2: %f\n",mf1[i],mf2[k]);
      /* ... check for matching */
      if(mf1[i]>=0.5||mf2[k]>=0.5){
	/* flag the primary */
	if(mf1[i]<=mf2[k]) {primary_1[i]=1; primary_2[k]=0;} 
	else { 
	  primary_2[k]=(primary_1[i]==0)?0:1; /* previous candidate was flagged as primary ?? */  
	  primary_1[i]=0;
	}
	++nmatches_1[i];
	++nmatches_2[k];
	id_match_1[i]=id_cls_2[k];
	id_match_2[k]=id_cls_1[i];
	matched[k]=1;
	mf3[n]=mf1[i];
	mf3[p]=mf2[k];

      } 
      id_match_3[n]=id_match_1[i];
      nmatches_3[n]=nmatches_1[i];
      primary_3[n]=primary_1[i];
      id_match_3[p]=id_match_2[k];
      nmatches_3[p]=nmatches_2[k];
      primary_3[p]=primary_2[k];
    } /* for(k=0;k<ncls_2;k++)*/   
  }
  free(memb1);
  free(memb2);
  if(Verbose) printf("...matching complete.\n");
}

void set_defaults(void){
  int i;
  Verbose=0;
  MoreVerbose=0;
  for(i=0;i<4;i++){
    ttype[i] = (char *) malloc(FLEN_VALUE);
    fname[i] = (char *) malloc(FLEN_VALUE);
  }
  ttype[0]="id";         fname[0]="clusters_1.fit"; 
  ttype[1]="host_id";    fname[1]="members_1.fit";
  ttype[2]="id_run1";    fname[2]="clusters_2.fit";
  ttype[3]="id";         fname[3]="members_2.fit";
}

void write_match_table(void){
  register int i,j;
  fitsfile *fptr;
  int status=0, hdutype;
  long firstrow, firstelem;
  int tfields   = 8;       
  char extname[] = "Matched_Clusters";
  char *ttype[] = { "id", "id_cand", "id_match", "mfrac", "nmatches", "primary",  "ngals", "run"};
  char *tform[] = { "1K", "1K"     , "1K"      , "1D"   , "1J"      , "1J"     ,  "1J"   , "1J" };
  char *tunit[] = { "\0", "\0"     , "\0"      , "\0"   , "\0"      , "\0"     ,  "\0"   , "\0" };
  REAL **Rrow;
  long long naxis2=0,frow=1,felem=1,nelem=1, colnum=1, **llrow;
  int **irow;
  ARRAY2(REAL,Rrow,1,1); // 2d array: ARRAY2(type,name,ncols,nrows)
  ARRAY2(long long,llrow,3,1);
  ARRAY2(int,irow,4,1);
  /* open file */
  char * pch;
  char outnameF[1000]="";
  pch=strrchr(fname[0],'.');
  strncat (outnameF, fname[0], pch-fname[0]);
  strcat(outnameF,"_match_table.fit");
  fits_create_file(&fptr,outnameF,&status);
  fits_create_tbl(fptr,BINARY_TBL,naxis2,tfields,ttype,tform,tunit,extname,&status);
  /* write file */
  for(i=0;i<ncls_3;i++){
    llrow[0][0]=id_cls_3[i];
    llrow[1][0]=id_cnd_3[i];
    llrow[2][0]=id_match_3[i];
    Rrow[0][0]=mf3[i];
    irow[0][0]=nmatches_3[i];
    irow[1][0]=primary_3[i];
    irow[2][0]=ngals_3[i];
    irow[3][0]=(i<ncls_1)?1:2;
    for(j=0;j<3;j++)
      fits_write_col_lnglng(fptr, colnum++,frow,felem,nelem,llrow[j],&status);
    fits_write_col_dbl(fptr,colnum++,frow,felem,nelem,Rrow[0],&status);
    for(j=0;j<4;j++)
      fits_write_col_int(fptr,colnum++,frow,felem,nelem,irow[j],&status);
    cfitsio_error(status);
    colnum=1;
    frow++;
  } 
  fits_close_file(fptr, &status);
  free2(Rrow);
  free2(llrow);
  free2(irow);
  if(Verbose) printf("match table saved.\n");
} /* write_match_table() */

void read_data(void) {
  register int i;
  fitsfile *fptr;
  int status=0,hdunum=2,hdutype,frow=1,felem=1,colnum,anynull;
  long nrows;
  long long longnull=0;

  if(Verbose) printf("reading catalogs...\n");

    /* read run1 cluster catalog */
    fits_open_file(&fptr,fname[0],READONLY,&status), cfitsio_error(status);
    fits_movabs_hdu(fptr,hdunum,&hdutype,&status);
    if(hdutype != BINARY_TBL)
      fprintf(stderr,"Error: Binary table expected in HDU #2.\n"),exit(1);
    fits_get_num_rows(fptr, &nrows, &status);
    ndata = (int) nrows;
    if(ndata == 0) fprintf(stderr,"Error: No data points in file.\n"),exit(1);
    ncls_1=ndata;
    id_cls_1 = (long long *) malloc(ndata*sizeof(long long));
    fits_get_colnum(fptr,CASEINSEN,ttype[0],&colnum,&status); //ttype[0]="id"
    fits_read_col(fptr,TLONGLONG,colnum,frow,felem,ndata,&longnull,id_cls_1,
		  &anynull,&status), cfitsio_error(status);
    fits_close_file(fptr,&status);
    if(Verbose) printf(" nclusters_1 = %d\n",ndata);

    /* read run1 members list */
    fits_open_file(&fptr,fname[1],READONLY,&status), cfitsio_error(status);
    fits_movabs_hdu(fptr,hdunum,&hdutype,&status);
    if(hdutype != BINARY_TBL)
      fprintf(stderr,"Error: Binary table expected in HDU #2.\n"),exit(1);
    fits_get_num_rows(fptr, &nrows, &status);
    ndata = (int) nrows;
    if(ndata == 0) fprintf(stderr,"Error: No data points in file.\n"),exit(1);
    nglx_1=ndata;
    id_glx_1 = (long long *) malloc(ndata*sizeof(long long));
    host_id_1 = (long long *) malloc(ndata*sizeof(long long));
    fits_get_colnum(fptr,CASEINSEN,ttype[0],&colnum,&status); //ttype[0]="id" 
    fits_read_col(fptr,TLONGLONG,colnum,frow,felem,ndata,&longnull,id_glx_1,
		  &anynull,&status), cfitsio_error(status);
    fits_get_colnum(fptr,CASEINSEN,ttype[1],&colnum,&status); //ttype[1]="host_id"
    fits_read_col(fptr,TLONGLONG,colnum,frow,felem,ndata,&longnull,host_id_1,
		  &anynull,&status), cfitsio_error(status);
    fits_close_file(fptr,&status);
    if(Verbose) printf(" ngalaxies_1 = %d\n",ndata);

    /* read run2 cluster catalog */
    fits_open_file(&fptr,fname[2],READONLY,&status), cfitsio_error(status);
    fits_movabs_hdu(fptr,hdunum,&hdutype,&status);
    if(hdutype != BINARY_TBL)
      fprintf(stderr,"Error: Binary table expected in HDU #2.\n"),exit(1);
    fits_get_num_rows(fptr, &nrows, &status);
    ndata = (int) nrows;
    if(ndata == 0) fprintf(stderr,"Error: No data points in file.\n"),exit(1);
    ncls_2=ndata;
    id_cls_2 = (long long *) malloc(ndata*sizeof(long long));
    fits_get_colnum(fptr,CASEINSEN,ttype[0],&colnum,&status); //ttype[0]="id"
    fits_read_col(fptr,TLONGLONG,colnum,frow,felem,ndata,&longnull,id_cls_2,
		  &anynull,&status), cfitsio_error(status);
    fits_close_file(fptr,&status);
    if(Verbose) printf(" nclusters_2 = %d\n",ndata);

    /* read run2 members list */
    fits_open_file(&fptr,fname[3],READONLY,&status), cfitsio_error(status);
    fits_movabs_hdu(fptr,hdunum,&hdutype,&status);
    if(hdutype != BINARY_TBL)
      fprintf(stderr,"Error: Binary table expected in HDU #2.\n"),exit(1);
    fits_get_num_rows(fptr, &nrows, &status);
    ndata = (int) nrows;
    if(ndata == 0) fprintf(stderr,"Error: No data points in file.\n"),exit(1);
    nglx_2=ndata;
    id_glx_2 = (long long *) malloc(ndata*sizeof(long long));
    host_id_2 = (long long *) malloc(ndata*sizeof(long long));
    fits_get_colnum(fptr,CASEINSEN,ttype[0],&colnum,&status); //ttype[0]="id" 
    fits_read_col(fptr,TLONGLONG,colnum,frow,felem,ndata,&longnull,id_glx_2,
		  &anynull,&status), cfitsio_error(status);
    fits_get_colnum(fptr,CASEINSEN,ttype[1],&colnum,&status); //ttype[1]="host_id"
    fits_read_col(fptr,TLONGLONG,colnum,frow,felem,ndata,&longnull,host_id_2,
		  &anynull,&status), cfitsio_error(status);
    fits_close_file(fptr,&status);
    if(Verbose) printf(" ngalaxies_2 = %d \n",ndata);

    if (Verbose) printf("...reading complete.\n");

} /* read_data() */

void cfitsio_error(int status){
  if(status) fits_report_error(stderr,status), exit(status);    
  return;
} 

void alloc_memory(void){
  ncls_3=ncls_1+ncls_2;
  id_cls_3 = (long long *) malloc(ncls_3*sizeof(long long));
  id_cnd_3 = (long long *) malloc(ncls_3*sizeof(long long));
  ngals_1 = (int *) malloc(ncls_1*sizeof(int));
  ngals_2 = (int *) malloc(ncls_2*sizeof(int));
  ngals_3 = (int *) malloc(ncls_3*sizeof(int));
  mf1 = (double *) malloc(ncls_1*sizeof(double));
  mf2 = (double *) malloc(ncls_2*sizeof(double));
  mf3 = (double *) malloc(ncls_3*sizeof(double));
  nmatches_1 = (int *) malloc(ncls_1*sizeof(int));
  nmatches_2 = (int *) malloc(ncls_2*sizeof(int));
  nmatches_3 = (int *) malloc(ncls_3*sizeof(int));
  primary_1 = (int *) malloc(ncls_1*sizeof(int));
  primary_2 = (int *) malloc(ncls_2*sizeof(int));
  primary_3 = (int *) malloc(ncls_3*sizeof(int));
  id_match_1 = (long long *) malloc(ncls_1*sizeof(long long));
  id_match_2 = (long long *) malloc(ncls_2*sizeof(long long));
  id_match_3 = (long long *) malloc(ncls_3*sizeof(long long));
  if(Verbose) printf("memory allocated.\n");
  return;
}

void free_memory(){
  free(ngals_1);    free(mf1);       free(nmatches_1); free(primary_1);
  free(ngals_2);    free(mf2);       free(nmatches_2); free(primary_2);
  free(ngals_3);    free(mf3);       free(nmatches_3); free(primary_3);
  free(id_match_1); free(id_cls_1);  free(id_glx_1);   free(host_id_1);
  free(id_match_2); free(id_cls_2);  free(id_glx_2);   free(host_id_2);
  free(id_match_3); free(id_cls_3);  free(id_cnd_3); 
  if(Verbose) printf("memory freed.\n");
  return;
}
