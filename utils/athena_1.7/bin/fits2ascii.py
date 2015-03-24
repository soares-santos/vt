#!/usr/bin/env python
#
# Requires python v3.x, but is backward-compatible with 2.x for x>6

# Martin Kilbinger 2014
# fits2ascii.py


# Compability with python2.x for x>6
from __future__ import print_function


import sys
import pyfits
from optparse import OptionParser



def read_fits_file(input, verbose):

    keys = ['OBJECT', 'TYPE', 'NRESAMPLE', 'UNITS']
    

    hdulist = pyfits.open(input)
    nhdu    = len(hdulist)
    for hdu in range(nhdu):
        hdu_type = type(hdulist[hdu])
        if hdu_type == pyfits.hdu.table.BinTableHDU or hdu_type == pyfits.hdu.table.TableHDU:
            if verbose: print('Table found in fits file {0} (hdu #{1})'.format(input, hdu))
            header_all_keys = hdulist[hdu].header.keys()
            header_col_keys = [s for s in header_all_keys if 'TTYPE' in s]

            # Print some header keys 
            for key in keys:
               if key in hdulist[hdu].header:
                print('# {0} = {1} ({2})'.format(key, hdulist[hdu].header.get(key), hdulist[hdu].header.comments[key]))

            # Print header
            print('#', end=' ')
            for key in header_col_keys:
                header_col_val = hdulist[hdu].header.get(key)
                print('{0:>18}'.format(header_col_val), end='')
            print()

            # Print data
            nrows = hdulist[hdu].header.get('NAXIS2')
            for n in range(nrows):
                print('  ', end='')
                for key in header_col_keys:
                    header_col_val = hdulist[hdu].header.get(key)
                    print('{0:18g}'.format(hdulist[hdu].data[header_col_val][n]), end='')
                print()

    hdulist.close()


####################
### Main program ###
####################

def main(argv=None):

    # Command line options
    usage  = "%prog [OPTIONS]"

    parser = OptionParser(usage=usage)
    parser.add_option('-i', '--input', dest='input', type='string', help='input xi FITS file name')
    #parser.add_option('-o', '--output', dest='output', type='string', default=output_base,
                      #help='output base name (default = {0})'.format(output_base))
    parser.add_option('-V', '--verbose', dest='verbose', action='store_true',
                      help='Verbose')


    options, args = parser.parse_args()

    see_help = 'See option \'-h\' for help.'

    if options.input is None:
        print('Input xi file not given (use option \'-i\'). ' + see_help, file=sys.stderr)
        return


    read_fits_file(options.input, options.verbose)



if __name__ == "__main__":
    sys.exit(main(sys.argv))


