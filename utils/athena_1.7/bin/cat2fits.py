cat  = 'data_0.4_0.6-proj.dat'

import pyfits
from mkstuff import *


data = read_table(cat)

c1 = pyfits.Column(name='x', format='D', array=data[:,0])
c2 = pyfits.Column(name='y', format='D', array=data[:,1])

c3 = pyfits.Column(name='w', format='D', array=data[:,2])

#c3 = pyfits.Column(name='z', format='D', array=data[:,2])
#c4 = pyfits.Column(name='kappa', format='D', array=data[:,3])
#c5 = pyfits.Column(name='e1', format='D', array=data[:,4])
#c6 = pyfits.Column(name='e2', format='D', array=data[:,5])

#hdu = pyfits.new_table([c1, c2, c3, c4, c5, c6])
hdu = pyfits.new_table([c1, c2, c3])
hdu.name = cat

phdu = pyfits.PrimaryHDU()
hdulist = pyfits.HDUList([phdu, hdu])

hdulist.writeto(cat + '.fits')

