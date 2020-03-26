import numpy as np
from astropy.io import fits
import sys
from glob import glob
import os

"""
This code use filtermask to mask out outside part of corrected_xxx.fits and wt.fits, set them to zero

"""


file0 = sys.argv[1]
wtfile = file0.replace('.fits','.wt.fits')

f = fits.open(file0)
w = fits.open(wtfile)
filt = f[0].header['filter'].lower()
month = f[0].header['date'].split('-')[1]


optic_dir = os.environ['_90PRIME_OPTIC_DIR']
if month == '02':
    month_dir = '/Feb/'
elif month == '03':
    month_dir = '/Mar/'
maskfile = optic_dir+month_dir+'optic_flat_'+filt+'_filtermask.fits'
mask = fits.open(maskfile)


print('Masking '+'\t'+file0.split('/')[-1]+'\twith'+'\t'+'optic_flat_'+filt+'_filtermask.fits')
hlistI = []
hlistW = []

for j in range(1,17):

    dati = f[j].data
    datw = w[j].data
    maskj = mask[j].data
    dati[maskj == 0] = 0
    datw[maskj == 0] = 0

    hduI = fits.ImageHDU()
    hduW = fits.ImageHDU()
    hduI.data = dati
    hduW.data = datw
    hduI.header = f[j].header
    hduW.header = w[j].header
    hduI.header['BZERO'] = 0.0
    hduW.header['BZERO'] = 0
    hlistI.append(hduI)
    hlistW.append(hduW)

hduI0 = fits.PrimaryHDU()
hduI0.header = f[0].header
hlistI.insert(0,hduI0)

hduW0 = fits.PrimaryHDU()
hduW0.header = w[0].header
hlistW.insert(0,hduW0)


hduIA = fits.HDUList(hlistI)
hduIA.writeto(file0.replace('corrected_','science_'))

hduWA = fits.HDUList(hlistW)
hduWA.writeto(wtfile.replace('corrected_','science_'))

print('Created: \t'+file0.replace('corrected_','science_').split('/')[-1] + '\t' + wtfile.replace('corrected_','science_').split('/')[-1])



