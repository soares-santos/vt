#!/usr/bin/env python
# Requires python v3.x, but is backwards compatible until python2.6

# Martin Kilbinger 2013
# test_suite_athena.py
# athena v1.7

# Runs tests to verify athena output to previous versions

# Compability with python2.x for x>5
from __future__ import print_function

import sys
import os
from optparse import OptionParser
import mkstuff


athena_version = 1.7


class test:
    """athena test runs
    """

    def __init__(self, dir, func, args, output, out_type, bindir):
        self.dir      = dir
        self.func     = func
        self.args     = args
        self.output   = output
        self.out_type = out_type
        self.bindir   = bindir


    def print(self):
        print(self.dir, ":", self.func, self.args)


    def chdir(self):
        pwd = os.getcwd()
        try:
            os.chdir(self.dir)
        except:
            print('chdir: Could not change into directory \'{0}\''.format(self.dir))
            return None

        return pwd


    def run(self):
        """Run test via shell command
        """

        pwd = self.chdir()
        if pwd is None: return -1
        res = mkstuff.run_cmd(self.bindir + '/' + self.func + ' ' + self.args)
        os.chdir(pwd)
        return res


    def write_file(self, i, path, fout):
        """Write test output file to results file
        """

        test_file = path + '/' + self.output[i]
        # Write file name
        print(test_file, file=fout, end='\n\n')

        extension = os.path.splitext(test_file)[1]
        if extension == '.fits' or extension == 'FITS':
            import subprocess
            prog = self.bindir + '/fits2ascii.py -i ' + test_file
            output = subprocess.check_output(prog.split(), shell=False)
            data = output.decode()
        else:
            fin  = open(test_file, 'r')
            data = fin.read()
            fin.close()
        #fout.write(data)
        print(data, file=fout)
        print(file=fout, end='\n')


    def append_output(self, b, fout):
        """Append test result to output file
        """

        print(b, self.dir, file=fout, end='\n\n')

        if self.output is None:
            print('No output file defined', file=fout)
            return

        pwd = self.chdir()
        for i in range(len(self.output)):
            self.write_file(i, '.', fout)
            self.write_file(i, 'results', fout)
            print(file=fout)
        os.chdir(pwd)



####################
### Main program ###
####################

def main(argv=None):

    # Command line options
    usage  = "%prog [OPTIONS]"

    ALL_TESTS = 1 | 2 | 4 | 8 | 16 | 32 | 64

    outbase_default = 'test_suite_results_files'

    # Path to this script
    bindir = os.path.abspath(os.path.dirname(__file__))

    # Path to test sub-directories
    testdir = './test'

    parser = OptionParser(usage=usage)
    parser.add_option('-l', '--list', dest='list', action='store_true', default=False, help='List all tests, do not run')
    parser.add_option('-t', '--test', dest='test', type='int', default=ALL_TESTS,
                      help='Run test TEST (bit-coded), see option \'-l\'). Default = {0}'.format(ALL_TESTS))
    parser.add_option('-P', '--path_athena', dest='path_to_athena', type='string', default=bindir,
                      help='Absolute path to programs \'athena\' and \'pallas.py\' (default: directory of this script)')
    parser.add_option('-T', '--path_test', dest='path_to_test', type='string', default=testdir,
                      help='Absolute path to test directory in which to find \'test_*\' (default: \'{0}\')'.format(testdir))
    parser.add_option('-o', '--output', dest='outbase', default=outbase_default,
                      help='Output file base (default:{0})'.format(outbase_default))
    options, args = parser.parse_args()


    # Definition of tests
    tests     = []
    bits      = []
    b         = 1


    ### Shear-shear (cosmic shear)

    # Shear two-point correlation function
    tests.append(test(options.path_to_test + '/test_xi', 'athena', '-c config_tree', ['xi'], '.ascii', options.path_to_athena))
    bits.append(b)
    b *= 2

    # Using a fits input catalogue
    tests.append(test(options.path_to_test + '/test_xi', 'athena', '-c config_fits', ['xi.fits'], 'fits', options.path_to_athena))
    bits.append(b)
    b *= 2

    # Aperture-mass dispersion and band-power spectrum
    tests.append(test(options.path_to_test + '/test_xi', 'pallas.py', '-i xi -w PbMp -v -o xi',
                      ['xi_map2_poly.txt', 'xi_pkappa_band.txt'], 'ascii', options.path_to_athena))
    bits.append(b)
    b *= 2

    # Using fits as input xi file
    tests.append(test(options.path_to_test + '/test_xi', 'pallas.py', '-i xi.fits -w PbMp -v -o xi.fits',
                      ['xi.fits_map2_poly.fits', 'xi.fits_pkappa_band.fits'], 'ascii', options.path_to_athena))
    bits.append(b)
    b *= 2


    ### Shear-position (galaxy-galaxy lensing)
    tests.append(test(options.path_to_test + '/test_g',  'athena', '-c config_tree', ['wgl'], 'ascii', options.path_to_athena))
    bits.append(b)
    b *= 2

    # Using fits input catalogues
    tests.append(test(options.path_to_test + '/test_g',  'athena', '-c config_fits', ['wgl.fits'], 'fits', options.path_to_athena))
    bits.append(b)
    b *= 2



    ### Position-position (galaxy angular clustering)
    output_w = []
    Nz       = 2
    estim    = ['LS', 'Ham']
    for i in range(Nz):
        for j in range(i, Nz):
            for e in estim:
                output_w.append('w_theta_{0}_{1}_{2}.dat'.format(i, j, e))

    tests.append(test(options.path_to_test + '/test_w',  'woftheta_xcorr.pl', '-a', output_w, '.ascii', options.path_to_athena))
    bits.append(b)
    b *= 2


    # List tests and return
    if options.list is True:
        for b, t in zip(bits, tests):
            #print(b, t.dir, t.func)
            print(b, " ", end="")
            t.print()
        return 0

    
    # Output text file
    outname = options.outbase + '.txt'
    fout = open(outname, 'w')
    print('test_suite_athena.py output file', file=fout)
    print('athena v.{0}'.format(athena_version), file=fout, end='\n\n\n')


    # Run tests
    for b, t in zip(bits, tests):
        if options.test & b:
            print(b, t.dir, t.func)
            res = t.run()
            if res == 0:
                t.append_output(b, fout)

    fout.close()


    # Finish
    print('Test results have been written to {0}'.format(outname))

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))


