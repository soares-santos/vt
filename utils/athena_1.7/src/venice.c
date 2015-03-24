#include "venice.h"

/*
------------------------
Jean Coupon - Jan 2009
------------------------

Programs that reads a mask file (DS9 type) and a catalogue 
of objects and computes one of the following tasks:
1. Creates a pixelized mask. 
   venice -m mask.reg [OPTIONS]
2. Finds objects inside/outside a mask. 
   venice -m mask.reg [OPTIONS]
3. Generates a random catalogue of objects inside/outside a mask. 
   venice -m mask.reg [OPTIONS]

Modifications:

v 2.01 - March 5 2009
- bug corrected in flagcat (wrong ymax)
v 2.0 - Feb 16th 2009
- All parameters put in a single structure "Config".
- File reading improved and simplified.
- Now integrates random catalogue on a sphere.
v 1.02 - Jan 17th 2009
- the option format [-f] can be use in task 2.
- the limits options can be used in task 2.
- README updated.
V 1.01 - April 2nd 2008
- The pixelized mask limits can be definied by the user as well
as the random catalogue limits.
- More explanations added in the README file.

V 1.00 - March 26th 2008
First version.
*/

int main(int argc, char **argv)
{
  //Initialization
  srand((unsigned int)time(NULL));
  EPS = determineMachineEpsilon();
  IDERR = determineSize_tError();
  Config para;
 
  //Tasks
  switch (readParameters(argc,argv,&para)){
  case 1:
    mask2d(&para);
    break;
  case 2:
    flagCat(&para);
    break;
  case 3:
    randomCat(&para);
    break;
  }
  
  return EXIT_SUCCESS;
}

/*----------------------------------------------------------------*
 *Main routines                                                   *
 *----------------------------------------------------------------*/

gsl_histogram2d *mask2d(const Config *para){
  /*Return the mask in gsl histogram format and writes the mask in fileOut. 
    The limits are the extrema of the extreme polygons in fileRegIn. 
    The pixel is set to 0 when inside the mask and 1 otherwise. 
  */  
  
  double x0,y0,x,y,xmin,xmax,ymin,ymax;
  int Npolys,poly_id;
  size_t i,j,count,total;

  
  FILE *fileRegIn = fopenAndCheck(para->fileRegInName,"r");
  FILE *fileOut = fopen(para->fileOutName,"w");
  
  //Reads the mask file (.reg DS9 file)
  //and construct the array of polygons.
  Polygon *polys=readPolygonFile(fileRegIn,&Npolys);
  
  //Sets the limit of the mask
  xmin = polys[0].xmin;
  xmax = polys[0].xmax;
  ymin = polys[0].ymin;
  ymax = polys[0].ymax;
  for(i=1;i<Npolys;i++){
    xmin = MIN(xmin,polys[i].xmin);
    xmax = MAX(xmax,polys[i].xmax);
    ymin = MIN(ymin,polys[i].ymin);
    ymax = MAX(ymax,polys[i].ymax);
  }
  
  //If the limits are definied by the user:
  if(para->minDefinied[0]) xmin = para->min[0];
  if(para->maxDefinied[0]) xmax = para->max[0];
  if(para->minDefinied[1]) ymin = para->min[1];
  if(para->maxDefinied[1]) ymax = para->max[1];

  x0 = xmin - 1.0; y0 = ymin - 1.0;//Reference point. It must be outside the mask;
  
  gsl_histogram2d *mask=gsl_histogram2d_alloc(para->nx,para->ny);//Mask.
  gsl_histogram2d_set_ranges_uniform(mask,xmin,xmax,ymin,ymax);//Ranges of the mask.
  
  total = para->nx*para->ny;

  printf("Progress =     ");
  for(i=0; i<mask->nx; i++){
    for(j=0; j<mask->ny; j++){
      count = i*para->ny+j;
      printCount(&count,&total,1000);
      //Center of the pixel.
      x = (mask->xrange[i]+mask->xrange[i+1])/2.0;
      y = (mask->yrange[j]+mask->yrange[j+1])/2.0;
      //1 = outside the mask,0 = inside the mask
      if(!insidePolygon(polys,Npolys,x0,y0,x,y,&poly_id)) mask->bin[i*mask->ny+j] = 1.0;
    }
  }
  fflush(stdout);
  printf("\b\b\b\b100%%\n");
  
  //Writes the file
  printf("Writing %s...\n",para->fileOutName);
  for(j=0; j<mask->ny; j++){
    for(i=0; i<mask->nx; i++){
      fprintf(fileOut,"%f ",mask->bin[i*mask->ny+j]);
    }
    fprintf(fileOut,"\n");
  }
  fclose(fileRegIn);
  fclose(fileOut);

  return mask; 
}


void flagCat(const Config *para){
  /*Reads fileCatIn and add a flag at the end of the line. 1 is outside
 the mask and 0 is inside the mask. xcol and ycol are the column ids of resp. x coordinate and y coordinate.*/
  double x0,y0,x,y,xmin,xmax,ymin,ymax;
  size_t i,N,Ncol;
  int Npolys,poly_id,flag;
  char line[NFIELD*NCHAR], item[NFIELD*NCHAR],*str_end;

  FILE *fileRegIn = fopenAndCheck(para->fileRegInName,"r");
  FILE *fileCatIn = fopenAndCheck(para->fileCatInName,"r");
  FILE *fileOut = fopen(para->fileOutName,"w");
  
  //Reads the mask file (.reg DS9 file)
  //and construct the array of polygons.
  Polygon *polys=readPolygonFile(fileRegIn,&Npolys);
  xmin = polys[0].xmin;
  xmax = polys[0].xmax;
  ymin = polys[0].ymin;
  ymax = polys[0].ymax;
  for(i=1;i<Npolys;i++){
    xmin = MIN(xmin,polys[i].xmin);
    xmax = MAX(xmax,polys[i].xmax);
    ymin = MIN(ymin,polys[i].ymin);
    ymax = MAX(ymax,polys[i].ymax);
  }
  
  //If the limits are definied by the user:
  if(para->minDefinied[0]) xmin = para->min[0];
  if(para->maxDefinied[0]) xmax = para->max[0];
  if(para->minDefinied[1]) ymin = para->min[1];
  if(para->maxDefinied[1]) ymax = para->max[1];
  

  x0 = xmin - 1.0; y0 = ymin - 1.0;//Reference point. It must be outside the mask;
  
  
  N = 0;
  while(fgets(line,NFIELD*NCHAR,fileCatIn) != NULL)
    if(getStrings(line,item," ",&Ncol))  N++;
  rewind(fileCatIn);
  printf("Nobjects = %zu\n",N);
  
  i = 0;
  printf("Progress =     ");
  while(fgets(line,NFIELD*NCHAR,fileCatIn) != NULL){
    if(getStrings(line,item," ",&Ncol)){
      i++;
      printCount(&i,&N,1000);
      x = getDoubleValue(item,para->xcol);
      y = getDoubleValue(item,para->ycol);
      flag=0;
      if(!insidePolygon(polys,Npolys,x0,y0,x,y,&poly_id) && 
	 xmin < x && x <= xmax && ymin < y && y <= ymax) flag = 1;
      str_end = strstr(line,"\n");//adds cariage return to the end of the line.
      strcpy(str_end,"\0");//adds the end symbol to the line
      switch (para->format){
      case 1://Only objects outside the mask
	if(flag) fprintf(fileOut,"%s\n",line);  
	break;
      case 2: //Only objects inside the mask
	if(!flag) fprintf(fileOut,"%s\n",line);  
	break;
      case 3://All objects with the flag
	fprintf(fileOut,"%s %d\n",line,flag);    
      }
    }
  }
  fflush(stdout);
  printf("\b\b\b\b100%%\n");

  fclose(fileOut);
  fclose(fileCatIn);
  return;
}


void randomCat(const Config *para){
  /*Generates a random catalogue inside the mask (uniform PDF).
    If "all", it puts all objects and add a flag such as: 
    outside the mask:1, inside the mask:0. Otherwise it puts only objects
    outside the mask.*/
  int Npolys,poly_id,flag;
  size_t i;
  double x0,y0,x,y,xmin,xmax,ymin,ymax;
  gsl_rng *r=randomInitialize(1);
 
  FILE *fileOut = fopen(para->fileOutName,"w");

  //If no mask is provided, it reads xmin, xmax, ymin and ymax and
  //make a random catalogue with no mask with the format "outside"
  //(all objects written with no flag: x y)...
  if(!strcmp(para->fileRegInName,"\0")){
    printf("Generating catalogue with no mask...\n");
    xmin = para->min[0];
    xmax = para->max[0];
    ymin = para->min[1];
    ymax = para->max[1];
    printf("xmin = %f, xmax = %f, ymin = %f, ymax = %f\n",xmin,xmax,ymin,ymax);
    printf("Progress =     ");
    for(i=0;i<para->npart;i++){
      printCount(&i,&para->npart,1000);
      if(para->coordType == CART){
	x = gsl_ran_flat(r,xmin,xmax);
	y = gsl_ran_flat(r,ymin,ymax);
	fprintf(fileOut,"%.8f %.8f\n",x,y);
      }else{
	x = gsl_ran_flat(r,xmin,xmax);
	y = gsl_ran_flat(r,sin(ymin*PI/180.0),sin(ymax*PI/180.0));
	fprintf(fileOut,"%.8f %.8f\n",x, asin(y)*180.0/PI);
      }
    }
    fflush(stdout);
    printf("\b\b\b\b100%%\n");
    return;
  }
  
  //... otherwise it reads the mask file (.reg DS9 file)
  //and construct the array of polygons.
  FILE *fileRegIn = fopenAndCheck(para->fileRegInName,"r");
  Polygon *polys=readPolygonFile(fileRegIn,&Npolys);
  xmin = polys[0].xmin;
  xmax = polys[0].xmax;
  ymin = polys[0].ymin;
  ymax = polys[0].ymax;
  for(i=1;i<Npolys;i++){
    xmin = MIN(xmin,polys[i].xmin);
    xmax = MAX(xmax,polys[i].xmax);
    ymin = MIN(ymin,polys[i].ymin);
    ymax = MAX(ymax,polys[i].ymax);
  }  

  //If the limits are definied by the user:
  if(para->minDefinied[0]) xmin = para->min[0];
  if(para->maxDefinied[0]) xmax = para->max[0];
  if(para->minDefinied[1]) ymin = para->min[1];
  if(para->maxDefinied[1]) ymax = para->max[1];
    
  x0 = xmin - 1.0; y0 = ymin - 1.0;//Reference point. It must be outside the mask;
  printf("Creates a random catalogue with N = %zu objects. Format = %d\n",para->npart,para->format);
  printf("xmin = %f, xmax = %f, ymin = %f, ymax= %f\n",xmin,xmax,ymin,ymax);
  printf("Progress =     ");
  for(i=0;i<para->npart;i++){
    printCount(&i,&para->npart,1000);
    x = gsl_ran_flat(r,xmin,xmax);
    y = gsl_ran_flat(r,ymin,ymax);
    if(para->coordType == CART){
      x = gsl_ran_flat(r,xmin,xmax);
      y = gsl_ran_flat(r,ymin,ymax);
    }else{
      x = gsl_ran_flat(r,xmin,xmax);
      y = gsl_ran_flat(r,sin(ymin*PI/180.0),sin(ymax*PI/180.0));
      y = asin(y)*180.0/PI;
    }
    
    if(flag=0,!insidePolygon(polys,Npolys,x0,y0,x,y,&poly_id)) flag = 1;
    switch (para->format){
    case 1://Only objects outside the mask
      if(flag) fprintf(fileOut,"%f %f\n",x,y);
      break;
    case 2: //Only objects inside the mask
      if(!flag) fprintf(fileOut,"%f %f\n",x,y);
      break;
    case 3://All objects with the flag
      fprintf(fileOut,"%f %f %d\n",x,y,flag);
    }
  }
  fflush(stdout);
  printf("\b\b\b\b100%%\n");

  return;
}


int readParameters(int argc, char **argv, Config *para){
  int i,task,nomask;
    
  //default parameters
  nomask = 1;
  task = 1;
  para->nx = 512;
  para->ny = 512;
  para->xcol = 1;
  para->ycol = 2;
  para->npart = 1000000;
  para->format = 1;
  para->coordType = CART;
    
  for(i=0;i<2;i++){
    para->minDefinied[i] = 0;
    para->maxDefinied[i] = 0;
    para->min[i] = 0.0;
    para->max[i] = 0.0;
  }
  
  strcpy(para->fileOutName,"mask.out");
  strcpy(para->fileCatInName,"\0");
  strcpy(para->fileRegInName,"\0");
  


  for(i=0;i<argc;i++){
    //Help-------------------------------------------------------------------//
    if(!strcmp(argv[i],"-h") || !strcmp(argv[i],"--help") || argc == 1){
      printf("\n\n                   V E N I C E\n\n");
      printf("           mask utility program version 2.01 \n\n");
      printf("Usage: %s -m mask.reg [OPTIONS]\n",argv[0]);
      printf("    or %s -m mask.reg -cat file.cat [OPTIONS]\n",argv[0]);
      printf("    or %s -m mask.reg -r [OPTIONS]\n",argv[0]);
      printf("Notice: 0 means inside the mask, 1 outside\n");
      exit(-1);
    }
    //Polygon file in--------------------------------------------------------//
    if(!strcmp(argv[i],"-m")){
      if(argv[i+1] == NULL){
	printf("Missing argument after %s\nExiting...\n",argv[i]);
	exit(-1);
      }
      strcpy(para->fileRegInName,argv[i+1]);
      nomask = 0;
    }
    //Input catalogue (if -cat set)------------------------------------------//
    if(!strcmp(argv[i],"-cat")){
      if(argv[i+1] == NULL){
	printf("Missing argument after %s\nExiting...\n",argv[i]);
	exit(-1);
      }
      strcpy(para->fileCatInName,argv[i+1]);
      task = 2;
    }
    //Random catalogue-------------------------------------------------------//
    if(!strcmp(argv[i],"-r")){
      task = 3;
    }
    //Output file------------------------------------------------------------//
    if(!strcmp(argv[i],"-o")){
      if(argv[i+1] == NULL){
	printf("Missing argument after %s\nExiting...\n",argv[i]);
	exit(-1);
      }
      strcpy(para->fileOutName,argv[i+1]);
    }
    //Dimensions of the pixel mask--------------------------------------------//
    if(!strcmp(argv[i],"-nx")){
      if(argv[i+1] == NULL){
	printf("Missing argument after %s\nExiting...\n",argv[i]);
	exit(-1);
      }
      para->nx = atoi(argv[i+1]);
    }
    if(!strcmp(argv[i],"-ny")){
      if(argv[i+1] == NULL){
	printf("Missing argument after %s\nExiting...\n",argv[i]);
	exit(-1);
      }
      para->ny = atoi(argv[i+1]);
    }
    //Column ids for the coordinates in FILE_CAT_IN--------------------------//
    if(!strcmp(argv[i],"-xcol")){
      if(argv[i+1] == NULL){
	printf("Missing argument after %s\nExiting...\n",argv[i]);
	exit(-1);
      }
      para->xcol = atoi(argv[i+1]);
    }
    if(!strcmp(argv[i],"-ycol")){
      if(argv[i+1] == NULL){
	printf("Missing argument after %s\nExiting...\n",argv[i]);
	exit(-1);
      }
      para->ycol = atoi(argv[i+1]);
    }
    //NPART for the random catalogue------------------------------------------//
    if(!strcmp(argv[i],"-npart")){
      if(argv[i+1] == NULL){
	printf("Missing argument after %s\nExiting...\n",argv[i]);
	exit(-1);
      }
      para->npart = atoi(argv[i+1]);
    }
    //OPTION for the random catalogue-----------------------------------------//
    if(!strcmp(argv[i],"-f")){
      if(argv[i+1] == NULL){
	printf("Missing argument after %s\nExiting...\n",argv[i]);
	exit(-1);
      }
      if(!strcmp(argv[i+1],"outside")) para->format = 1;
      if(!strcmp(argv[i+1],"inside")) para->format = 2;
      if(!strcmp(argv[i+1],"all")) para->format = 3;
    }
    //Limits for the random catalogue-----------------------------------------//
    if(!strcmp(argv[i],"-xmin")){
      if(argv[i+1] == NULL){
	printf("Missing argument after %s\nExiting...\n",argv[i]);
	exit(-1);
      }
      para->minDefinied[0] = 1;
      para->min[0] = atof(argv[i+1]);
    }
    if(!strcmp(argv[i],"-xmax")){
      if(argv[i+1] == NULL){
	printf("Missing argument after %s\nExiting...\n",argv[i]);
	exit(-1);
      }
      para->maxDefinied[0] = 1;
      para->max[0] = atof(argv[i+1]);
    }
    if(!strcmp(argv[i],"-ymin")){
      if(argv[i+1] == NULL){
	printf("Missing argument after %s\nExiting...\n",argv[i]);
	exit(-1);
      }
      para->minDefinied[1] = 1;
      para->min[1] = atof(argv[i+1]);
    }
    if(!strcmp(argv[i],"-ymax")){
      if(argv[i+1] == NULL){
	printf("Missing argument after %s\nExiting...\n",argv[i]);
	exit(-1);
      }
      para->maxDefinied[1] = 1;
      para->max[1] = atof(argv[i+1]); 
    }
    if(!strcmp(argv[i],"-coord")) {
      if(!strcmp(argv[i+1],"spher")) para->coordType = RADEC;
      if(!strcmp(argv[i+1],"cart")) para->coordType = CART;
    }


  }
  //If no mask file is provided-----------------------------------------//
  
  if (nomask){
    //Checks if all the limits are definied in this case;
    if(task == 3 && (para->minDefinied[0] + para->minDefinied[1] + para->maxDefinied[0] + para->maxDefinied[1] < 4)) {
      printf("If you want to generate a random catalogue with no mask,\n");
      printf("please provide all the coordinate limits:\n");
      printf("%s -r -xmin value -xmax value -ymin value -ymax value [OPTIONS]\n",argv[0]);
      exit(-1);
    }
    //Checks if the limits are realistics (min < max)
    if(task == 3 && ( para->min[0] > para->max[0] || para->min[1] > para->max[1])){
      printf("Please put realistic limits (xmin < xmax and ymin < ymax).\n");
      exit(-1);
    }
    if(task != 3){
      printf("please provide a mask file.\n");
      printf("Usage: %s -m mask.reg [OPTIONS]\n",argv[0]);
      printf("    or %s -m mask.reg -cat file.cat [OPTIONS]\n",argv[0]);
      printf("    or %s -m mask.reg -r [OPTIONS]\n",argv[0]);
      exit(-1); 
    }
  }
  //----------------------------------------------------------------------//

  return task;
}

/*----------------------------------------------------------------*
 *Utils - geometric                                               *
 *----------------------------------------------------------------*/

int insidePolygon(Polygon *polys,int Npolys, double x0,double y0,double x,double y,int *poly_id){
  /*Returns 1 if the point (x,y) is inside one of the polygons in
    polys. The id of the first polygon in which the point is 
    found is also returned in poly_id. If the point is found to be outside 
    all polygons, the function returns 0 and poly_id is set to -1.
    The function first checks if the point is inside the square drawn 
    by the extrema of each polygon ("computationaly" more effecient). 
    If yes, it counts how many times the segment {(x0,y0),(x,y)} crosses 
    the sides of the polygon (Ncross). (x,y) inside the polygon 
    implies Ncross is an odd number if the point (x0,y0) is outside the polygon
    (then the point (x0,y0) must be chosen outside any polygon).*/
  int i,j,Ncross;
  double s,t,D;
  
  for(i=0;i<Npolys;i++){
    if(polys[i].xmin < x && x < polys[i].xmax && polys[i].ymin < y && y < polys[i].ymax){
      Ncross=0;
      for(j=0;j<polys[i].N;j++){
	if(j<polys[i].N-1){
	  D = (polys[i].x[j+1]-polys[i].x[j])*(y-y0)-(polys[i].y[j+1]-polys[i].y[j])*(x-x0);
	  s = ((x-x0)*(polys[i].y[j]-y)-(y-y0)*(polys[i].x[j]-x))/D;
	  t = ((polys[i].x[j]-x)*(polys[i].y[j+1]-polys[i].y[j])-(polys[i].y[j]-y)*(polys[i].x[j+1]-polys[i].x[j]))/D;
	}else{
	  D = (polys[i].x[0]-polys[i].x[j])*(y-y0)-(polys[i].y[0]-polys[i].y[j])*(x-x0);
	  s = ((x-x0)*(polys[i].y[j]-y)-(y-y0)*(polys[i].x[j]-x))/D;
	  t = ((polys[i].x[j]-x)*(polys[i].y[0]-polys[i].y[j])-(polys[i].y[j]-y)*(polys[i].x[0]-polys[i].x[j]))/D;
	} 
	if(0.0 < s && s < 1.0 && 0.0 < t && t < 1.0) Ncross++;	
      }
      if(GSL_IS_ODD(Ncross)){
	*poly_id = i;
	return 1;
      }
    }
  }
  
  *poly_id = -1;
  return 0;
}

Polygon *readPolygonFile(FILE *fileIn,int *Npolys){
  /*Reads the file file_in and returns the array of 
    polygons.*/
  char line[NFIELD*NCHAR], item[NFIELD*NCHAR],*str_begin,*str_end;
  int i,j;
  size_t N;
  
  *Npolys = 0;
  //Reads the entire file and count the total number of polygons, Npolys.
  while(fgets(line,NFIELD*NCHAR,fileIn) != NULL)
    if(strstr(line,"polygon") != NULL) *Npolys += 1;
  rewind(fileIn);
  Polygon *polys=malloc(*Npolys*sizeof(Polygon));
  
  i=0;
  //Reads the file and fills the array with polygons.
  while(fgets(line,NFIELD*NCHAR,fileIn) != NULL){
    if(strstr(line,"polygon") != NULL){
      str_begin = strstr(line,"(")+sizeof(char);
      str_end = strstr(line,")");
      strcpy(str_end,"\n\0");
      strcpy(line,str_begin);
      getStrings(line,item,",",&N);
      //------------------------------//
      //gets all coordinates separated by comas.
      polys[i].N = N/2;
      if(N/2 > 100){
	printf("%s: %zu = too many points for polygon %d (100 maxi). Exiting...\n",MYNAME,N/2,i);
	exit(-1);
      }      
      polys[i].xmin = atof(item);
      polys[i].xmax = atof(item);
      polys[i].ymin = atof(item+NCHAR);
      polys[i].ymax = atof(item+NCHAR);
      for(j=0;j<N/2;j++){
	polys[i].x[j] = atof(item+NCHAR*2*j);
	polys[i].y[j] = atof(item+NCHAR*(2*j+1));
	polys[i].xmin = MIN(polys[i].xmin, polys[i].x[j]);
	polys[i].xmax = MAX(polys[i].xmax, polys[i].x[j]);
	polys[i].ymin =  MIN(polys[i].ymin, polys[i].y[j]);
	polys[i].ymax =  MAX(polys[i].ymax, polys[i].y[j]);
      }
      
      i++;
      //------------------------------//
    }
  }
  if(i==0){
    printf("%s: 0 polygon found, check input file. Exiting...\n",MYNAME);
    exit(-1);
  }
  
  return polys;
}


/*----------------------------------------------------------------*
 *Utils - numeric                                                 *
 *----------------------------------------------------------------*/

gsl_rng *randomInitialize(int time_dependent){
  /*Define here which type of random generator you want. Set 
    time_dependent to 1 if you want the seed to change each time
    you run the program (change every 1 s). */
  
  /*
    WARNING: the time_dependent option actually doesn't work yet :-(, sorry !
    To initialize with a different seed, type :
    tcsh: setenv GSL_RNG_SEED `date +%s`
    bash: export GSL_RNG_SEED=`date +%s`
  */
  //time_dependent option
  //if(time_dependent) system("export GSL_RNG_SEED=20091982");
  //else system("export GSL_RNG_SEED=20091982");

  
  //Random number generator
  const gsl_rng_type *T;
  gsl_rng *r;
  struct tms Time;
  clock_t cl;

  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);

  // MK NEW:
  if (time_dependent) {
     cl = times(&Time);
     gsl_rng_set(r, cl);
  }

  return r;
}


double determineMachineEpsilon()
{
  double u, den;
  
  u = 1.0;
  do {
    u /= 2.0;
    den = 1.0 + u;
  } while(den>1.0);

  return(10.0 * u);
}


size_t determineSize_tError(){
  /*Equals size_t_max in fact.*/
  size_t count=0;
  count--;
  //count = (size_t)(-1);
  return count;
}

FILE *fopenAndCheck(const char *fileName,char *mode){
  /*Checks if fileName exists and opens it. Exits otherwise.;*/
  FILE *fileTmp = fopen(fileName,mode);
  if (fileTmp == NULL){
    printf("%s: %s not found. Exiting...\n",MYNAME,fileName);
    exit(-1);    
  }
  return fileTmp;
}

int getStrings(char *line, char *strings, char *delimit, size_t *N){
  /*Extract each word/number in line separated by delimit and returns 
    the array of items in strings.*/
  int i,j,begin,length;
  
  if(line == NULL || line[0] == '\n' || line[0] == '#' || line[0] == '\0') return 0;
  
  i = 0;
  j = 0;
  while(line[i] != '\0' && line[i] != '#' && line[i] != '\n'){
    begin = i;
    while(line[i] == *delimit || line[i] == '\t' && (line[i] != '\0' || line[i] != '#' || line[i] != '\n')) i++;
    begin = i;
    while(line[i] != *delimit && line[i] != '\t' && line[i] != '\0' && line[i] != '#' && line[i] != '\n') i++;
    length = i - begin;
    if(length > 0){
      strncpy(strings+NCHAR*j,&line[begin],length);
      strcpy(strings+NCHAR*j+length,"\0");
      j++;
    }
  }
 
  (*N) = j;

  if(*N > 0){
    return 1;
  }else{
    return 0;
  }
}


void printCount(const size_t *count, const size_t *total, const size_t step){
  if((*count)%step == 0){
    fflush(stdout);
    printf("\b\b\b\b%3.0f%%",100.0*(double)(*count)/(double)(*total));
  }
  return;
}

