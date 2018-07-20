# coding: utf-8
import numpy as np
import matplotlib.image as mpimg
import pylab as pl
import sys


pngfname = 'herapixels.png'
npzfname = 'hera_txt.npz'
jdate = 2458098.27471265

im = mpimg.imread("herapixels.png")
arr = im[:,:,0]
arr = arr[::-1,:]

nx, ny = arr.shape
y, x = np.where(arr == 0)
dx = x - int(np.floor(nx/2.))
dy = y - int(np.floor(ny/2.))

pixsize_deg = 1.0

za = np.sqrt(dx**2 + dy**2)*pixsize_deg
aztmp= np.arctan2(dy,dx)*180./np.pi
az = np.arctan2(dy,dx) * 180./np.pi - 90.


alts = 90. - zas
Nsrcs = zas.size
fluxes = np.ones_like(azs)

catalog = []

source_coord = SkyCoord(alt=Angle(alts, unit=units.deg), az=Angle(azs, unit=units.deg),
                        obstime=time, frame='altaz', location=array_location)

icrs_coord = source_coord.transform_to('icrs')

np.savez('hera_txt.npz', az=az, za=za, units='degrees')

az = np.radians(az)
za = np.radians(za)
aztmp = np.radians(aztmp)
pl.scatter(za*np.cos(az), za*np.sin(az), label='Rotated')
pl.scatter(za*np.cos(aztmp), za*np.sin(aztmp), label='Original')
pl.legend()
pl.show()

#im = np.zeros((int(im.shape[1] * pixsize_deg+1), int(im.shape[0]* pixsize_deg)+1))
#
#im = np.zeros(100,100)
#extent = (im.shape[1]*pixsize_deg+1, im.shape[0].pixsize_deg+1)
#degpix = 
#
#xind, yind = (za*np.cos(az)*pixsize_deg + 6).astype(int), (za*np.sin(az)*pixsize_deg + 1).astype(int)
#
#im[xind, yind] += 1
