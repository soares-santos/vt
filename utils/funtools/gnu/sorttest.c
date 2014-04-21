/* 
   export EVENTS="(dval:D,fval:E,ival:J,uival:V,sval:I,cval:B,cval2:B)"
   sorttest | _sort -B24 +... | fundisp "stdin[EVENTS()]"
*/

#include <stdio.h>
#include <stdlib.h>
#include <salloc.h>

extern char *optarg;
extern int optind;

#define NREC 10
#define SZ_LINE 1024

#define MAX(a,b) (((a)>(b))?(a):(b))

typedef struct brec{
  double dval;
  float fval;
  int ival;
  unsigned int uival;
  short sval;
  char cval;
  char cval2;
} Binary;

void
fillb(Binary *b, int i, int dodouble)
{
  int idiv;
  b->dval = (double)rand();
  b->fval = (float)rand();
  if( dodouble ){
    if( (i%2) == 0 )
      b->ival = rand();
  }
  else{
    b->ival = rand();
  }
  b->uival = (unsigned)rand();
  idiv = MAX(1,(RAND_MAX/32768));
  b->sval = rand()/idiv;
  idiv = MAX(1,(RAND_MAX/128));
  b->cval = rand()/idiv;
  b->cval2 = rand()/idiv;
}

int main(int argc, char **argv)
{
  int i;
  int c;
  int ival;
  int ichan, ochan, pid, status;
  int nrec=NREC;
  int dodouble=0;
  int doxeq=0;
  int narg=0;
  char *prog="_sort";
  char *slist=NULL;
  char *s;
  char tbuf[SZ_LINE];
  char cmd[SZ_LINE];
  char *args[SZ_LINE];
  FILE *fp=stdout;
  Binary b;

  /* process switch arguments */
  while ((c = getopt(argc, argv, "dn:s:x")) != -1){
    switch(c){
    case 'd':
      dodouble = 1;
      break;
    case 'n':
      nrec = atoi(optarg);
      break;
    case 's':
      doxeq = 1;
      slist=(char *)NewString(optarg);
      break;
    case 'x':
      doxeq = 1;
      break;
    }
  }
  
  /* see rand */
  srand(time(NULL));

  /* generate sort parameters and start sort program */
  if( doxeq ){
    args[narg++] = (char *)NewString(prog);
    sprintf(tbuf, "-B%d", sizeof(Binary));
    args[narg++] = (char *)NewString(tbuf);
    if( !slist ){
      slist = (char *)NewString("d");
    }
    for(s=slist; *s; s++){
      switch(*s){
      case 'd':
	sprintf(tbuf, "+d0");
	args[narg++] = (char *)NewString(tbuf);
	break;
      case 'f':
	sprintf(tbuf, "+f8");
	args[narg++] = (char *)NewString(tbuf);
	break;
      case 'i':
	sprintf(tbuf, "+i12");
	args[narg++] = (char *)NewString(tbuf);
	break;
      case 'I':
	sprintf(tbuf, "+I16");
	args[narg++] = (char *)NewString(tbuf);
	break;
      case 's':
	sprintf(tbuf, "+s20");
	args[narg++] = (char *)NewString(tbuf);
	break;
      case 'c':
	sprintf(tbuf, "+c22");
	args[narg++] = (char *)NewString(tbuf);
	break;
      default:
	fprintf(stderr, "ERROR: unknown sort type: %s\n", s);
	exit(1);
      }
      if( narg > (SZ_LINE-2) ){
	fprintf(stderr, "ERROR: too many arguments: %d\n", narg);
	exit(1);
      }
    }
    args[narg] = NULL;
    if( !ProcessOpen(prog, args, &ichan, &ochan, &pid) ){
      fprintf(stderr, "ERROR: can't start %s\n", prog);
      exit(1);
    }
    /* write unsorted records */
    for(i=0; i<nrec; i++){
      fillb(&b, i, dodouble);
      write(ochan, &b, sizeof(Binary));
    }
    /* signal EOF to sort program */
    close(ochan);
    while( read(ichan, &b, sizeof(Binary)) == sizeof(Binary) )
      write(1, &b, sizeof(Binary));
    if( slist ) free(slist);
    ProcessClose(pid, &status);
  }
  else{
    /* write unsorted records */
    for(i=0; i<nrec; i++){
      fillb(&b, i, dodouble);
      fwrite(&b, sizeof(Binary), 1, fp);
    }
  }
}
