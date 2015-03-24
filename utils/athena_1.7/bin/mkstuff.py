# mkstuff.py module
#
# Martin Kilbinger 2011
#


# Compability with python2.x for x>6
from __future__ import print_function

import re
import numpy as np
import math
import os
import errno
import subprocess
import shlex
import sys


###########
### I/O ###
###########


def read_table(fname, count_comment=False):
    """Reads the file 'fname', skips comment lines
    (starting with '#') and returns a numpy matrix.
    If count_comment=True, also returns number of comment
    lines and the header.
    """

    data = []
    header = []
    ncomment = 0
    f = open(fname, "r")
    all_lines = f.readlines()
    f.close()
    for line in all_lines:
        m = re.search("#", line)
        if m is None:
            data.append([])
            all_el = line.split()
            for el in all_el:
                data[-1].append(float(el))
        else:
            ncomment += 1
            header.append(line)
    mdata = np.matrix(data)

    if count_comment == True:
        return mdata, ncomment, header
    else:
        return mdata


# Returns the column number for field 'name', from header
def get_column(name, header, verbose=True):
    header_fields = header[0].replace('#', '').split()
    for i in range(len(header_fields)):
        if name == header_fields[i]: return i
    if verbose is True:
       warning('Field {0} not found in header'.format(name)) 
    return -1


# Runs shell command
def run_cmd(cmd, run=True, verbose=True, stop=False):
    cmds = shlex.split(cmd)

    if run is True:
        if verbose is True:
            print('Running command \'{0}\''.format(cmd))
        try:
            ex = subprocess.call(cmds)
        except OSError as e:
            print('Error: {0}'.format(e.strerror))
            ex = e.errno
    else:
        if verbose is True:
            print('Not running command \'{0}\''.format(cmd))
        return 0

    if ex != 0:
        if verbose is True:
            print('Last command returned {0}'.format(ex), end='')
            if stop is True:
                print(', stopping')
            else:
                print(', continuing')

        if stop is True:
            sys.exit(ex)

    return ex


# Like 'mkdir -p'
def mkdir_p(path):
    try:
        os.makedirs(path)
    #except OSError, exc: # Python <2.5
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST:
            pass
        else: raise


def write_matrix(a, name):
    """Writes the matrix a to the file name.
    """

    Nx, Ny = a.shape
    f = open(name, 'w')
    for x in range(Nx):
        for y in range(Ny):
            print(a[x,y] + ' ', file=f, end='')
        print('', file=f)
    f.close()


def log_command(argv, name=None):
    if name is None:
        name = 'log_' + os.path.basename(argv[0])
    f = open(name, 'w')
    for a in argv:
        print(a, end='', file=f)
        print(' ', end='', file=f)
    print('', file=f)
    f.close()


def error(str, val=1):
    print(str, file=sys.stderr)
    if val is not None: sys.exit(val)


def warning(str):
    error('Warning: ' + str, None)


############
### Maths ##
############


# Returns density value levels of L corresponding
# to confidence levels cl
def get_density_levels(L, cl=None):
    if cl == None: cl   = [0.6827, 0.9545, 0.9973]

    Ls   = np.sort(L, axis=None)[::-1]
    Lsum = Ls.sum()
    N    = Ls.shape[0]

    cum  = 0
    levels = np.zeros(shape=len(cl)) - 1
    for i in range(N):
        cum = cum + Ls[i]
        for l in range(len(cl)):
            if cum >= cl[l] * Lsum and levels[l] < 0:
                levels[l] = Ls[i]
    return levels


# x can be list or np.array
def mean(x):
    return float(sum(x))/len(x)


def corr_coeff(a):
    Nx, Ny = a.shape
    ra     = np.zeros(shape=a.shape) 
    for i in range(Nx):
        for j in range(Ny):
            ra[i,j] = a[i,j] / np.sqrt(a[i,i] * a[j,j])
    return ra


def frexp10(x):
    """Return the mantissa and exponent of x, as pair (m, e).
    m is a float and e is an int, such that x = m * 10.0**e.
    See math.frexp()"""
    if x == 0: return (0, 0)
    try:
        l = math.log10(abs(x))
    except:
        print('Error with math.log10(|' + (str(x)) + '|)')
        return None, None
    if l < 1: l = l - 1 + 1e-10
    exp = int(l)
    return x / 10**exp, exp 

