#!/usr/bin/env python

# Requires python v3.x, but is backward-compatible with 2.x for x>6

# Martin Kilbinger 2013
# pallas.py
# v1.0
#
# Power-spectrum for kAppa of eLL bAnd eStimator
# P_kAppa(eLL)^bAnd.eStimator
#
# Estimates the band-power spectrum C_ell of kappa from
# the shear 2PCF as in SvWKM02.
# Calculates E-, B- and EB-power spectra.
# Also calculates the aperture-mass dispersion from the 2PCF.
#
# pallas.py is part of athena, http://cosmostat.org/athena.html .
# This version of pallas.py is compatible with all athena versions


# Compability with python2.x for x>6
from __future__ import print_function


import sys
import numpy
from optparse import OptionParser
from mkstuff import *
from math import *
from scipy.special import j0, j1, jn

try:
    import pyfits as fits
except ImportError:
    pass
    try:
        from astropy.io import fits
    except ImportError:
        error("Could not import pyfits astropy.io/fits library")

import warnings



#################
### Constants ###
#################

units       = {'rad'    : 1,
               'deg'    : np.pi / 180.0,
               'arcmin' : np.pi / 180.0 / 60.0,
               'arcsec' : np.pi / 180.0 / 60.0 / 60.0,
               'none'   : 1}
unit_default = 'arcmin'


###############
### Classes ###
###############

class in_cols(object):
    """ Default names of input columns. Indices are for ascii files.
    """

    xi_names       = ['theta', 'xi_p', 'xi_m']
    xi_indices     = [0, 1, 2]

    xi_names_opt   = ['xi_x']
    xi_indices_opt = [3]

    gl_names       = ['theta', 'g_t']
    gl_indices     = [0, 1]

    gl_names_opt   = ['g_x']
    gl_indices_opt = [2]

    w_names        = ['theta', 'w']
    w_indices      = [0, 1]


class xi_data:
    """Correlation function data
    """

    def __init__(self, theta, xip, xim, xix):
        self.theta   = theta
        self.xip     = xip
        self.xim     = xim
        self.xix     = xix
        self.length  = len(theta)
        self.binning = None
        self.delta   = None


    def set_binning(self, force_reg=False, verbose=False):
        """Sets irregular (recommended), linear, or logarithmic binning
        """

        if (self.length < 2): error('Data vector has length {0}, has to be larger than 2'.format(self.length))

        if force_reg == False:
            self.binning, self.delta = 'ireg', 0.0
            return

        eps     = 1.0e-2

        dtheta1 = self.theta[1] - self.theta[0]
        dtheta2 = self.theta[self.length-1] - self.theta[self.length-2]
        drel = abs(dtheta2 - dtheta1)/self.theta[1]
        if drel < eps:
            self.binning, self.delta = 'lin', dtheta2
            return

        dlogtheta1 = log(self.theta[1]) - log(self.theta[0])
        dlogtheta2 = log(self.theta[self.length-1]) - log(self.theta[self.length-2])
        drel = abs(dlogtheta2 - dlogtheta1)
        if  drel < eps:
            self.binning, self.delta = 'log', dlogtheta2
            return

        if verbose == True:
            print('Bins seem not to be regular, falling back to \'ireg\'')
        self.binning, self.delta = 'ireg', 0.0

    def get_binning(self):
        """ Returns the binning type and delta (0 if 'ireg')
        """

        if self.binning is None or self.delta is None:
            error('Binning not set')
        return self.binning, self.delta

    def dtheta_theta(self, i):
        """Returns dtheta_i * theta_i
        """

        if self.binning is 'lin':
            # dtheta * theta
            d = self.delta * self.theta[i]
        elif self.binning is 'log':
            # d log theta * theta^2 = d theta / theta * theta^2 = dtheta * theta
            d = self.delta * self.theta[i] * self.theta[i]
        elif self.binning is 'ireg':
            if i != self.length-1:
                d = (self.theta[i+1] - self.theta[i]) * self.theta[i]
            else:           
                d = (self.theta[i] - self.theta[i-1]) * self.theta[i]

        return d



class pkappa_data:
    """Data for P_kappa as fct. of ell, the band-power spectrum P_i,
    and smooth real-space functions such as <M_ap^2>(theta).
    """

    def __init__(self, ell_min, ell_max, band=False, Nell=None, unit_out=None):
        """If Nell is not None, a band-power spectrum (or smooth real-space function)
        is initialised with Nell logarithmic bins.
        """

        if Nell is None:
            self.ell  = np.arange(ell_min, ell_max)
        else:
            self.ell  = np.logspace(log10(ell_min), log10(ell_max), Nell)

        self.pE   = np.zeros(np.shape(self.ell))
        self.pB   = np.zeros(np.shape(self.ell))
        self.pEB  = np.zeros(np.shape(self.ell))
        self.Nell = len(self.ell)
        self.band = band
        if unit_out is None:
            self.unit_out = 'none'
        else:
            self.unit_out = unit_out

    def write(self, fname, file_format, header=None, verbose=False):
        """Writes the power spectrum/smooth function to the ascii or fits
           file 'fname.[txt|fits]'
        """

        if file_format is 'ascii':
            ffname = fname + '.txt'
            self.write_ascii(ffname, header)
        elif file_format is 'fits':
            ffname = fname + '.fits'
            self.write_fits(ffname, header)
        else:
            error('Wrong file format \'{0}\', has to be \'ascii\' or \'fits\' (option \'-f\')'.format(file_format))
        
        if verbose == True: print('Writing output file \'{0}\''.format(ffname))


    def write_ascii(self, fname, header=None):
        """Writes the power spectrum/smooth function to the ascii file 'fname'
        """

        f = open(fname, 'w')
        if header is not None:
            if self.unit_out is not 'none':
                my_header = header.split()
                my_header[0] = '{0}[{1}]'.format(my_header[0], self.unit_out)
                header = ' '.join(my_header)
            f.write(header + '\n')

        for j in range(self.Nell):
            scale = rad_to_unit(self.ell[j], self.unit_out)
            f.write('{0:10.3f} {1: 12.5e} {2: 12.5e} {3: 12.5e}'.format(scale, self.pE[j], self.pB[j], self.pEB[j]))
            if self.band == True:
                ell_l, ell_u = self.ell_l_u(j)
                ell_l = rad_to_unit(ell_l, self.unit_out)
                ell_u = rad_to_unit(ell_u, self.unit_out)
                f.write('{0:10.3f} {1:10.3f}'.format(ell_l, ell_u))
            f.write('\n')

        f.close()


    def write_fits(self, fname, header=None):
        """Writes the power spectrum/smooth function to the ascii file 'fname'
        """

        n_col = 4
        if self.band is True: n_col += 2
        if header is None:
            my_header = ['' for x in range(n_col)]
        else:
            my_header = header.replace('#', '').split()

        scales = rad_to_unit(self.ell, self.unit_out)
        col0 = fits.Column(name=my_header[0], format='E', array=scales)
        col1 = fits.Column(name=my_header[1], format='E', array=self.pE)
        col2 = fits.Column(name=my_header[2], format='E', array=self.pB)
        col3 = fits.Column(name=my_header[3], format='E', array=self.pEB)
        if self.band is True:
            ell_l = []
            ell_u = []
            for j in range(self.Nell):
                my_ell_l, my_ell_u = self.ell_l_u(j)
                my_ell_l = rad_to_unit(my_ell_l, self.unit_out)
                my_ell_u = rad_to_unit(my_ell_u, self.unit_out)
                ell_l.append(my_ell_l)
                ell_u.append(my_ell_u)
            col4 = fits.Column(name=my_header[4], format='E', array=ell_l)
            col5 = fits.Column(name=my_header[5], format='E', array=ell_u)
            cols = fits.ColDefs([col0, col1, col2, col3, col4, col5])
        else:
            cols = fits.ColDefs([col0, col1, col2, col3])


        hdu  = fits.new_table(cols)
        if self.unit_out is not 'none':
            hdu.header.append(card=('UNITS', self.unit_out, 'Coordinate units'))
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            hdu.writeto(fname, clobber=True)
    

    def ell_l_u(self, j):
        """Returns the lower and upper band ell
        """

        dlogell = log(self.ell[1]) - log(self.ell[0])
        return self.ell[j] / sqrt(exp(dlogell)), self.ell[j] * sqrt(exp(dlogell))


#################
### Functions ###
#################


def unit_to_rad(theta, unit):
    return theta * units[unit]


def rad_to_unit(theta, unit):
    return theta / units[unit]


def check_fields(fields, strings, indices, exit, verbose=False):
    """Prints error (warning) message if header does not contain necessary (optional) fields
    """

    if fields is None:
        if exit is True:
            warning('No header found, using default columns (0:theta, 1:xi+, 2:xi-)')
        return False

    for i in range(len(strings)):
        if len(fields) <= indices[i] or not strings[i] in fields[indices[i]]:
            msg = 'Field \'{0}\' not found in header (expected in field #{1})'.format(strings[i], indices[i]-1)
            if exit is True:
                error(msg)
            else:
                if verbose is True:
                    warning(msg)
            return False

    return True


def get_unit(field, verbose=False):
    """Returns unit for angular scales, from header line with format e.g. 'theta[UNIT]'
    """

    for u in units.keys():
        if u in field:
            if verbose == True:
                print('Unit of input angular scales = {0}'.format(u))
            return u

    if verbose == True:
        warning('No unit for angular scales found, assuming \'{0}\''.format(unit_default))

    return unit_default


def get_unit_fits(header, verbose=False):
    if 'UNITS' in header.keys():
        unit = header['UNITS']
    else:
        if verbose == True:
            warning('No \'UNITS\' keyword in fits header found, assuming default unit for angular scales \'{0}\''.format(unit_default))
            return unit_default

    # Go through valid units and compare to assigned one (from header)
    for u in units.keys():
        if u == unit:
            if verbose == True:
                print('Unit of input angular scales = {0}'.format(u))
            return u

    if verbose == True:
        warning('Unit from fits header (\'UNITS\' = \'{0}\') not a valid unit, reverting to default unit \'{1}\''.format(unit, unit_default))
        return unit_default


def get_file_format(file_format, fname):
    """Returns file_format if not None, or 'fits' if extension is '.fits' or '.FITS', and 'ascii' otherwise.
    """

    if file_format is None:
        # Determine format using file extension
        extension = os.path.splitext(fname)[1]
        if extension == '.fits' or extension == '.FITS':
            file_format = 'fits'
        else:
            # Default
            file_format = 'ascii'

    elif file_format != 'fits' and file_format != 'ascii':
       error('Wrong file format \'{0}\', has to be \'ascii\' or \'fits\' (option \'-f\')'.format(file_format))

    return file_format


def read_xi(fname, file_format, force_reg, verbose=False):
    """Reads a shear-shear correlation function file (ascii or fits).
    """

    my_file_format = get_file_format(file_format, fname)

    if verbose:
        print("File format is {0}".format(my_file_format))

    if my_file_format is 'fits':
        xi = read_xi_fits(fname, force_reg, verbose)
    elif my_file_format is 'ascii':
        xi = read_xi_ascii(fname, force_reg, verbose)
    else:
       error('Wrong file format \'{0}\', has to be \'ascii\' or \'fits\' (option \'-f\')'.format(my_file_format))

    return xi, my_file_format


def read_xi_fits(fname, force_reg, verbose=False):
    """Read a shear-shear correlation function file in fits format.
    """

    hdulist = fits.open(fname)
    nhdu    = len(hdulist)
    for hdu in range(nhdu):
        hdu_type = type(hdulist[hdu])

        if hdu_type == fits.hdu.table.BinTableHDU or hdu_type == fits.hdu.table.TableHDU:
            if verbose: print('Table found in fits file {0} (hdu #{1})'.format(input, hdu))
            header_all_keys = hdulist[hdu].header.keys()
            header_col_keys = [s for s in header_all_keys if 'TTYPE' in s]
            header_col_vals = [hdulist[hdu].header.get(s) for s in header_col_keys]

            # Required fields
            check_fields(header_col_vals, in_cols.xi_names, in_cols.xi_indices, True, verbose)
            theta = hdulist[hdu].data[in_cols.xi_names[0]]
            xip   = hdulist[hdu].data[in_cols.xi_names[1]]
            xim   = hdulist[hdu].data[in_cols.xi_names[2]]

            # Optional fields
            has_xix = check_fields(header_col_vals, in_cols.xi_names_opt, in_cols.xi_indices_opt, True, verbose)
            if has_xix == True:
                xix = hdulist[hdu].data[in_cols.xi_names_opt[0]]
            else:
                xix = None

            # Units of angular scales
            unit  = get_unit_fits(hdulist[hdu].header, verbose)
            theta = unit_to_rad(theta, unit)

            break


    hdulist.close()

    xi = xi_data(theta, xip, xim, xix)
    xi.set_binning(force_reg, verbose)
    if verbose is True: print('Binning type of 2PCF is \'{0}\''.format(xi.binning))

    return xi


def read_xi_ascii(fname, force_reg, verbose=False):
    """Reads a shear-shear correlation function ascii file.
    """

    if verbose is True: print('Reading input xi file \'{0}\''.format(fname))

    xi, ncomment, header  = read_table(fname, True)

    if xi.shape[1] <= 2:
        error('Input file \'{0}\' has only {1} columns. Minimum columns are (theta, xi+, xi-)'.format(fname, xi.shape[1]))

    if ncomment > 0:
        fields = header[0].replace('#', '').split()
    else:
        fields = None

    # Required fields
    check_fields(fields, in_cols.xi_names, in_cols.xi_indices, True, verbose)

    theta = xi[:,in_cols.xi_indices[0]]
    xip   = xi[:,in_cols.xi_indices[1]]
    xim   = xi[:,in_cols.xi_indices[2]]

    # Optional fields
    has_xix = check_fields(fields, in_cols.xi_names_opt, in_cols.xi_indices_opt, False, verbose)
    if has_xix == True: 
        xix  = xi[:,in_cols.xi_indices_opt[0]]
    else:
        xix  = None

    # Units of angular scales
    if fields is not None:
        unit = get_unit(fields[in_cols.xi_indices[0]], verbose)
    else:
        unit = unit_default
    theta  = unit_to_rad(theta, unit)

    xi = xi_data(theta, xip, xim, xix)
    xi.set_binning(force_reg, verbose)
    if verbose is True: print('Binning type of 2PCF is \'{0}\''.format(xi.binning))

    return xi


def read_gl(fname, file_format, force_reg, verbose=False):
    """Reads a shear-position correlation function file (fits or ascii format).
    """

    my_file_format = get_file_format(file_format, fname)
    
    if verbose:
        print("File format is {0}".format(my_file_format))
    
    if my_file_format is 'fits':
        wgl = read_gl_fits(fname, force_reg, verbose)
    elif my_file_format is 'ascii':
        wgl = read_gl_ascii(fname, force_reg, verbose)
    else:
       error('Wrong file format \'{0}\', has to be \'ascii\' or \'fits\' (option \'-f\')'.format(my_file_format))

    return wgl, my_file_format


def read_gl_fits(fname, force_reg, verbose=False):
    """Reads a shear-position correlation function file in fits format.
    """

    hdulist = fits.open(fname)
    nhdu    = len(hdulist)
    for hdu in range(nhdu):
        hdu_type = type(hdulist[hdu])

        if hdu_type == fits.hdu.table.BinTableHDU or hdu_type == fits.hdu.table.TableHDU:
            if verbose: print('Table found in fits file {0} (hdu #{1})'.format(input, hdu))
            header_all_keys = hdulist[hdu].header.keys()
            header_col_keys = [s for s in header_all_keys if 'TTYPE' in s]
            header_col_vals = [hdulist[hdu].header.get(s) for s in header_col_keys]

            # Required fields
            check_fields(header_col_vals, in_cols.gl_names, in_cols.gl_indices, True, verbose)
            theta = hdulist[hdu].data[in_cols.gl_names[0]]
            gt    = hdulist[hdu].data[in_cols.gl_names[1]]

            # Optional fields
            has_gx = check_fields(header_col_vals, in_cols.gl_names_opt, in_cols.gl_indices_opt, True, verbose)
            if has_gx == True:
                gx = hdulist[hdu].data[in_cols.gl_names_opt[0]]
            else:
                gx = None

            # Units of angular scales
            unit  = get_unit_fits(hdulist[hdu].header, verbose)
            theta = unit_to_rad(theta, unit)

            break


    hdulist.close()

    wgl = xi_data(theta, gt, gx, None)
    wgl.set_binning(force_reg, verbose)
    if verbose is True: print('Binning type of <gt> is \'{0}\''.format(wgl.binning))

    return wgl


def read_gl_ascii(fname, force_reg, verbose=False):
    """Reads a shear-position correlation function file in ascii format.
    """

    if verbose is True: print('Reading input wgl file \'{0}\''.format(fname))

    wgl, ncomment, header  = read_table(fname, True)

    if wgl.shape[1] <= 1:
        error('Input file \'{0}\' has only {1} columns. Minimum columns are (theta, gt)'.format(fname, wgl.shape[1]))

    if ncomment > 0:
        fields = header[0].replace('#', '').split()
    else:
        fields = None

    # Required fields
    check_fields(fields, in_cols.gl_names, in_cols.gl_indices, True, verbose)

    theta = wgl[:,in_cols.gl_indices[0]]
    gt    = wgl[:,in_cols.gl_indices[1]]

    # Optional fields
    has_gx  = check_fields(fields, in_cols.gl_names, in_cols.gl_indices, False, verbose)
    if has_gx  == True:
        gx   = wgl[:,in_cols.gl_indices_opt[0]]
    else:
        gx   = None

    # Units of angular scales
    if fields is not None:
        unit = get_unit(fields[in_cols.gl_indices[0]], verbose)
    else:
        unit = unit_default
    theta  = unit_to_rad(theta, unit)

    wgl = xi_data(theta, gt, gx, None)
    wgl.set_binning(force_reg, verbose)
    if verbose is True: print('Binning type of <gt> is \'{0}\''.format(wgl.binning))

    return wgl


def read_w(fname, file_format, force_reg, verbose=False):
    """Reads a spatial correlation function file, format fits or ascii.
    """

    my_file_format = get_file_format(file_format, fname)
            
    if verbose:
        print("File format is {0}".format(my_file_format))
    
    if my_file_format is 'fits':
        w = read_w_fits(fname, force_reg, verbose)
    elif my_file_format is 'ascii':
        w = read_w_ascii(fname, force_reg, verbose)
    else:
       error('Wrong file format \'{0}\', has to be \'ascii\' or \'fits\' (option \'-f\')'.format(my_file_format))
            
    return w, my_file_format


def read_w_fits(fname, force_reg, verbose=False):
    """Read a spatial correlation function file in fits format.
    """

    hdulist = fits.open(fname)
    nhdu    = len(hdulist)
    for hdu in range(nhdu):
        hdu_type = type(hdulist[hdu])

        if hdu_type == fits.hdu.table.BinTableHDU or hdu_type == fits.hdu.table.TableHDU:
            if verbose: print('Table found in fits file {0} (hdu #{1})'.format(input, hdu))
            header_all_keys = hdulist[hdu].header.keys()
            header_col_keys = [s for s in header_all_keys if 'TTYPE' in s]
            header_col_vals = [hdulist[hdu].header.get(s) for s in header_col_keys]

            # Required fields
            check_fields(header_col_vals, in_cols.w_names, in_cols.w.indices, True, verbose)
            theta  = hdulist[hdu].data[in_cols.w_names[0]]
            w      = hdulist[hdu].data[in_cols.w_names[1]]

            # Units of angular scales
            unit  = get_unit_fits(hdulist[hdu].header, verbose)
            theta = unit_to_rad(theta, unit)

            break

    hdulist.close()

    w = w_data(theta, w, None, None)
    w.set_binning(force_reg, verbose)
    if verbose is True: print('Binning type of 2PCF is \'{0}\''.format(w.binning))

    return w


def read_w_ascii(fname, force_reg, verbose=False):
    """Reads a spatial correlation function file in ascii format fits or ascii.
    """

    if verbose is True: print('Reading input w file \'{0}\''.format(fname))

    w, ncomment, header  = read_table(fname, True)

    if w.shape[1] <= 1:
        error('Input file \'{0}\' has only {1} columns. Minimum columns are (theta, w)'.format(fname, w.shape[1]))

    if ncomment > 0:
        fields = header[0].replace('#', '').split()
    else:
        fields = None

    # Required fields
    check_fields(fields, in_cols.w_names, in_cols.w_indices, True, verbose)

    theta = w[:,in_cols.w_indices[0]]
    w     = w[:,in_cols.w_indices[1]]

    # Units of angular scales
    if fields is not None:
        unit = get_unit(fields[in_cols.w_indices[0]], verbose)
    else:
        unit = unit_default
    theta  = unit_to_rad(theta, unit)

    w = xi_data(theta, w, None, None)
    w.set_binning(force_reg, verbose)
    if verbose is True: print('Binning type of w is \'{0}\''.format(w.binning))

    return w

 
def pkappa_ell(xi, ell_min, ell_max):
    """Convergence power spectrum for every ell, SvWKM02 (45)
    """

    pkappa = pkappa_data(ell_min, ell_max)

    for j, ell in enumerate(pkappa.ell):

        for i, theta in enumerate(xi.theta):

            d = xi.dtheta_theta(i)
            pref = np.pi * d

            J0 = j0(theta * ell)
            J4 = jn(4, theta * ell)
            fp = xi.xip[i] * J0
            fm = xi.xim[i] * J4

            pkappa.pE[j]  += pref * (fp + fm)
            pkappa.pB[j]  += pref * (fp - fm)
            if xi.xix is not None:
                fx = xi.xix[i] * J4
                pkappa.pEB[j] += 2.0 * pref * fx

    return pkappa


def g(x, pm):
    """For pm=+1, -1: SvWKM02 (50)
    """

    if pm == +1:
        return x * j1(x)
    elif pm == -1:
        return (x - 8.0/x) * j1(x) - 8 * jn(2, x)
    elif pm == 0:
        return x * j1(x) - 2 * j0(x)
    else:
        error('g: pm has to be +1 or -1')


def pkappa_band_power(xi, ell_min, ell_max, Nell):
    """Band-power convergence spectrum <l^2/2pi P(l)>, SvWKM02 (49)
    """

    pkappa  = pkappa_data(ell_min, ell_max, True, Nell)
    dlogell = log(pkappa.ell[1]) - log(pkappa.ell[0])
    p0      = 1 / (2 * dlogell)     # 2pi/Delta / 2 / 2pi

    for j in range(pkappa.Nell):

        ell_l, ell_u = pkappa.ell_l_u(j)

        for i in range(xi.length):

            d    = xi.dtheta_theta(i)
            pref = p0 * d / xi.theta[i]**2

            Gp = g(xi.theta[i] * ell_u, +1) - g(xi.theta[i] * ell_l, +1)
            Gm = g(xi.theta[i] * ell_u, -1) - g(xi.theta[i] * ell_l, -1)
            fp = xi.xip[i] * Gp
            fm = xi.xim[i] * Gm

            pkappa.pE[j]  += pref * (fp + fm)
            pkappa.pB[j]  += pref * (fp - fm)
            if xi.xix is not None:
                fx = xi.xix[i] * Gm    # CHECK!
                pkappa.pEB[j] += 2.0 * pref * fx

    return pkappa


def pnkappa_band_power(xi, ell_min, ell_max, Nell):
    """Band-power overdensity-convergence cross-spectrum
    """

    pnkappa = pkappa_data(ell_min, ell_max, True, Nell)
    dlogell = log(pnkappa.ell[1]) - log(pnkappa.ell[0])
    p0      = 1 / dlogell 

    for j in range(pnkappa.Nell):

        ell_l, ell_u = pnkappa.ell_l_u(j)

        for i in range(xi.length):

            d    = xi.dtheta_theta(i)
            pref = p0 * d / xi.theta[i]**3

            G0 = g(xi.theta[i] * ell_u, 0) - g(xi.theta[i] * ell_l, 0)
            fp = xi.xip[i] * G0
            fm = xi.xim[i] * G0

            pnkappa.pE[j]  += pref * fp
            pnkappa.pB[j]  += pref * fm

    return pnkappa


def pn_band_power(xi, ell_min, ell_max, Nell):
    """Band-power galaxy overdensity spectrum, analogous to SvWKM02 (49)
    """

    pn      = pkappa_data(ell_min, ell_max, True, Nell)
    dlogell = log(pn.ell[1]) - log(pn.ell[0])
    p0      = 1 / dlogell 

    for j in range(pn.Nell):

        ell_l, ell_u = pn.ell_l_u(j)

        for i in range(xi.length):

            d    = xi.dtheta_theta(i)
            pref = p0 * d / xi.theta[i]**2

            Gp = g(xi.theta[i] * ell_u, +1) - g(xi.theta[i] * ell_l, +1)
            fp = xi.xip[i] * Gp

            pn.pE[j]  += pref * fp

    return pn



def T(x, pm, mfilter):
    """Filter function for aperture-mass dispersion
    """

    if pm != +1 and pm != -1:
        error('T: pm has to be +1 or -1')

    xsqr = x**2

    if mfilter is 'poly':
        if x >= 2: return 0

        if pm == +1:
            p = 6.0/5.0 * (2.0 - 15.0*xsqr) * (1.0 - 2.0/np.pi * math.asin(x/2.0)) \
                + x * sqrt(4.0 - xsqr) / (100.0 * np.pi) \
                      * (120.0 + xsqr * (2320.0 + xsqr * (-754.0 + xsqr * (132.0 - 9.0 * xsqr))))
        elif pm == -1:
            p = 192.0 / (35.0*np.pi) * xsqr * x * pow(1.0 - xsqr/4.0, 3.5)
        return p

    elif mfilter is 'Gauss':
        f = exp(-xsqr/4.0) / 128.0
        if pm == +1:
            p = xsqr * (xsqr - 16.0) + 32.0
        elif pm == -1:
            p = xsqr**2
        return p * f

    else:
        error('Filter {0} not supported'.format(mfilter))



def map2(xi, theta_min, theta_max, Ntheta, mfilter):
    """Aperture-mass dispersion
    """

    m2        = pkappa_data(theta_min, theta_max, False, Ntheta, unit_default)
    dlogtheta = log(m2.ell[1]) - log(m2.ell[0])

    for j, theta in enumerate(m2.ell):

        p = 1.0 / (2.0 * theta**2)

        for i, vartheta in enumerate(xi.theta):

            d    = xi.dtheta_theta(i)
            pref = p * d

            tp = T(vartheta / theta, +1, mfilter)
            tm = T(vartheta / theta, -1, mfilter)
            fp = xi.xip[i] * tp
            fm = xi.xim[i] * tm

            m2.pE[j]  += pref * (fp + fm)
            m2.pB[j]  += pref * (fp - fm)
            if xi.xix is not None:
                fx = xi.xix[i] * tm
                m2.pEB[j] += 2.0 * pref * fx

    return m2



def parse_options(f_theta_min, f_theta_max, Ntheta):
    """Parses command line options.
    """

    usage  = "%prog [OPTIONS]"

    # Default parameters
    output_base  = 'output'
    ell_info     = '100 1000 10'
    which        = 'PbMgMp'
    wcorr        = 1

    parser = OptionParser(usage=usage)
    parser.add_option('-i', '--input', dest='input', type='string', help='input xi file name (containing at least the columns theta, xi+, xi-)')
    parser.add_option('-o', '--output', dest='output_base', type='string', default=output_base,
                      help='output base name (default = {0})'.format(output_base))
    parser.add_option('-f', '--format', dest='format', type='string', default=None,
                      help='input (and output) file formats [fits|ascii]. (Default: fits if extension is \'.fits\', ascii otherwise.)')
    parser.add_option('-W', '--WCORR', dest='wcorr', type='int', default=wcorr,
                      help='type of correlation (bit-coded): 1=shear-shear, 2=shear-position, 4=position=position. Default={0}'.format(wcorr))
    parser.add_option('-w', '--which', dest='which', type='string', default=which,
                      help='which second-order function? \'Pl\' (P_kappa(ell)), \'Pb\' (band-power P_kappa), '
                            'Mp\' (<Map^2> polynomial filter), ' +
                           '\'Mg\' (<Map^2 Gaussian filter). Default = \'{0}\''.format(which))
    parser.add_option('-l', '--ell_info', dest='ell_info', type='string', default=ell_info,
                      help='ELL_INFO = \'ell_min ell_max N_ell_band\' (default = \'{0}\')'.format(ell_info))
    parser.add_option('-t', '--theta_info', dest='theta_info', type='string',
                      help='THETA_INFO = \'theta_min theta_max N_theta\' (default: xi.theta_min*{0} '
                           'xi.theta_max/{1} {2}'.format(f_theta_min, f_theta_max, Ntheta))
    parser.add_option('-H', '--no_header', dest='no_header', action='store_true', help='do not write header to output files')
    parser.add_option('-v', '--verbose', dest='verbose', action='store_true', help='verbose')
    parser.add_option('-F', '--force_reg', dest='force_reg', action='store_true', default=False,
                      help='force regular binning if possible')

    options, args = parser.parse_args()

    return options


def check_options(options):
    """Checks command line options. Returns False in case of an
    invalid options, and True otherwise.
    """

    see_help = 'See option \'-h\' for help.'

    if options.input is None:
        print('Input correlation file not given (use option \'-i\'). ' + see_help, file=sys.stderr) 
        return False

    return True



####################
### Main program ###
####################

def main(argv=None):
    """Main program of pallas.py.
    """

    
    # Default parameters needed in main
    f_theta_min  = 5
    f_theta_max  = 2
    Ntheta       = 20


    # Command line options
    options = parse_options(f_theta_min, f_theta_max, Ntheta)

    if check_options(options) is False:
        return 1

    # Save calling command
    log_command(argv)


    # For power spectra
    ell_min, ell_max, Nell_band = map(int, options.ell_info.split())

    # Headers
    if options.no_header:
        (header_M, header_Pl, header_Pb) = (None, None, None)
    else:
        header_M  = '# theta <M_ap^2> <M_x^2> <M_apM_x>'
        header_Pl = '# ell P_E P_B P_EB'
        header_Pb = header_Pl + ' ell_lower ell_upper'


    if options.wcorr & 1:

        # Input xi file
        xi, in_format = read_xi(options.input, options.format, options.force_reg, options.verbose)


        # Aperture scales
        if options.theta_info is None:
            theta_min = xi.theta[0] * f_theta_min
            theta_max = xi.theta[xi.length - 1] / f_theta_max
        else:
            theta_min, theta_max, Ntheta = map(int, options.theta_info.split())

        # <M_ap^2>
        if 'Mp' in options.which:
            if options.verbose is True:
                print('Calculating <M_ap^2> (polynomial filter)')
            m2 = map2(xi, theta_min, theta_max, Ntheta, 'poly')
            m2.write(options.output_base + '_map2_poly', in_format, header_M, options.verbose)

        if 'Mg' in options.which:
            if options.verbose is True:
                print('Calculating <M_ap^2> (Gaussian filter)')
            m2 = map2(xi, theta_min, theta_max, Ntheta, 'Gauss')
            m2.write(options.output_base + '_map2_gauss', in_format, header_M, options.verbose)

        # P_kappa
        if 'Pl' in options.which:
            if options.verbose is True:
                print('Calculating P_kappa(ell)')
            pkappa = pkappa_ell(xi, ell_min, ell_max)
            pkappa.write(options.output_base + '_pkappa_ell', in_format, header_Pl, options.verbose)

        if 'Pb' in options.which:
            if options.verbose is True:
                print('Calculating band-power P_kappa')
            pkappa = pkappa_band_power(xi, ell_min, ell_max, Nell_band)
            pkappa.write(options.output_base + '_pkappa_band', in_format, header_Pb, options.verbose)


    if options.wcorr & 2:

        # Input wgl file
        wgl, in_format = read_gl(options.input, options.format, options.force_reg, options.verbose)

        if 'Pb' in options.which:
            if options.verbose is True:
                    print('Calculating band-cross-power P_nkappa')
            pnkappa = pnkappa_band_power(wgl, ell_min, ell_max, Nell_band)
            pnkappa.write(options.output_base + '_pnkappa_band', in_format, header_Pb, options.verbose)


    if options.wcorr & 4:

        # Input w file
        w, in_format = read_w(options.input, options.format, options.force_reg, options.verbose)

        if 'Pb' in options.which:
            if options.verbose is True:
                print('Calculating band-power P_n')
            pn = pn_band_power(w, ell_min, ell_max, Nell_band)

            pn.write(options.output_base + '_pn_band', in_format, header_Pb, options.verbose)

        

if __name__ == "__main__":
    sys.exit(main(sys.argv))

