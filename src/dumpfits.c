#include <string.h>
#include <stdio.h>
#include "fitsio.h"

int main(int argc, char *argv[])
{
    fitsfile *fptr;      /* FITS file pointer, defined in fitsio.h */
    char *val, value[1000], nullstr[]="*";
    char keyword[FLEN_KEYWORD], colname[FLEN_VALUE];
    int status = 0, hdutype, ncols, ii, anynul, dispwidth[1000];
    int firstcol, lastcol = 0, linewidth;
    long jj, nrows;

    if (argc != 2) {
      printf("Usage: fits2ascii filename[ext][col filter][row filter] \n");
      printf("\n");
      printf("Transform a FITS table into an ascii file\n");
      printf("\n");
      printf("Example: \n");
      printf("  fits2ascii tab.fits[GTI]           - list the GTI extension\n");
      printf("  fits2ascii tab.fits[1][#row < 101] - list first 100 rows\n");
      printf("  fits2ascii tab.fits[1][col -PI]    - don't list the PI column\n");
      printf("  fits2ascii tab.fits[1][col -PI][#row < 101]  - combined case\n");
      printf("Display formats can be modified with the TDISPn keywords.\n");
      return(0);
    }

    if (!fits_open_file(&fptr, argv[1], READONLY, &status))
    {
      val = value;
      fits_get_hdu_type(fptr, &hdutype, &status); /* Get the HDU type */

      if (hdutype == IMAGE_HDU)   /* primary array or image extension */
        printf("Error: a table extension must be specified.\n");

      else  /* a table extension */
      {
        fits_get_num_rows(fptr, &nrows, &status);
        fits_get_num_cols(fptr, &ncols, &status);

        while(lastcol < ncols) {
          linewidth = 0;
          firstcol = lastcol+1;
          for (lastcol = firstcol; lastcol <= ncols; lastcol++) {
             fits_get_col_display_width
	       (fptr, lastcol, &dispwidth[lastcol], &status);
             linewidth += dispwidth[lastcol] + 1;
	     /*             if (linewidth > 80)break;*/
          }
          if (lastcol > firstcol)lastcol--;  /* the last col didn't fit */

          /* print column names as column headers */
	  /* printf("#"); */
/*           for (ii = firstcol; ii <= lastcol; ii++) { */
/*             fits_make_keyn("TTYPE", ii, keyword, &status); */
/*             fits_read_key(fptr, TSTRING, keyword, colname, NULL, &status); */
/*             colname[dispwidth[ii]] = '\0';  /*
/* truncate long names */ 
/*             printf("%*s ",dispwidth[ii], colname);  */
/*           } */
/*           printf("\n");*/  
	     /* terminate header line */ 

          /* print each column, row by row (there are faster ways to do this) */
          for (jj = 1; jj <= nrows && !status; jj++) {
            for (ii = firstcol; ii <= lastcol; ii++)
            {
              /* read value as a string, regardless of datatype */
              fits_read_col_str
                   (fptr,ii,jj, 1, 1, nullstr, &val, &anynul, &status);
              if (!status)printf("%-*s ",dispwidth[ii], value);
            }
            printf("\n");
          }
        }
      }
      fits_close_file(fptr, &status);
    } 

    /* if error occured, print out error message */
    if (status) fits_report_error(stderr, status);
    return(status);
}
