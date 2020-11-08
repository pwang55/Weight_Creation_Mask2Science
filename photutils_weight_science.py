"""

Usage:

    In Script folder:
    $ python photutils_weight_science.py path_to_file/filename

    In Script folder, if running command line for loop, path to the files need to be provided so that the script can find the files:
    $ for files in `cat correctedlist.txt`; do
      > python photutils_weight_science.py $files path_to_files
      > done

    In files folder:
    $ python photutils_weight_science.py filename

    To change backsize, use option: backsize=64

This script takes corrected fits image as input, use photutils to calculate background RMS map, convert it to weight map, 
and mask corrected image and weight image with mask, output to science frames. Lightstreaks will be detected and masked in the weight image.


"""
import numpy as np
from astropy.io import fits
from astropy.stats import SigmaClip
from photutils import Background2D, MedianBackground
from photutils import detect_sources
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from photutils import source_properties
import sys
import os



backsize = 64
filtersize = 3
streak_sigma_threshold = 1.5
elongation_threshold = 10.0
npixels = 50
sigma_clip = SigmaClip(sigma=3.)
bkg_estimator=MedianBackground()


file0=sys.argv[1]
f=fits.open(file0)
filename=file0.split('/')[-1]
path=file0[: - len(filename)]

if len(sys.argv) > 2:
    for i in range(1, len(sys.argv)):
        arg = sys.argv[i]

        if len(arg.split('/')) > 1:
            path = arg
            if path[-1] != '/':
                path = path + '/'

        elif arg.split('=')[0] == 'backsize':
            backsize = int(arg.split('=')[1])


science_filename=filename.replace('corrected_', 'science_')
weight_filename=science_filename.replace('.fits', '.wt.fits')

optic_dir=os.environ['_90PRIME_OPTIC_DIR']
filt=f[0].header['filter'].lower()
month=f[0].header['date'].split('-')[1]
if month == '02':
    month_dir = '/Feb/'
elif month == '03':
    month_dir = '/Mar/'
maskfile = optic_dir + month_dir + 'optic_flat_' + filt + '_filtermask.fits'
mask = fits.open(maskfile)

print('Creat weight map for: ' + filename)
print('Filtermask: ' + month_dir + 'optic_flat_' + filt + '_filtermask.fits')

hlistI = []
hlistW = []
hlistS = []
streak_count = 0

for j in range(1, 17):
    dat = f[j].data
    hdr = f[j].header
    maskj = mask[j].data
    bkg = Background2D(dat, box_size=backsize, filter_size=filtersize, sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    wt = 1 / bkg.background_rms ** 2

    dat[maskj == 0] = 0
    wt[maskj == 0] = 0

    threshold = bkg.background + (streak_sigma_threshold * bkg.background_rms)
    sigma = 2.0 * gaussian_fwhm_to_sigma
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    kernel.normalize()
    segm = detect_sources(dat, threshold, npixels=npixels, filter_kernel=kernel)
    cat = source_properties(dat, segm)
    tab = cat.to_table()
    elongation = tab['elongation'].value
    streak_id = tab['id'][elongation > elongation_threshold]
    seg2 = segm.data.copy()
    if len(streak_id) > 0:
        streak_count = 1
        for i in range(len(streak_id)):
            wt[seg2 == streak_id[i]] = 0.0
            seg2[seg2 == streak_id[i]] = 10000.0

    hduI = fits.ImageHDU()
    hduW = fits.ImageHDU()
    hduI.data = dat
    hduW.data = wt
    hduI.header = hdr
    hduW.header = hdr
    hduI.header['BZERO'] = 0.0
    hduW.header['BZERO'] = 0.0
    hlistI.append(hduI)
    hlistW.append(hduW)

    seg2[seg2 < 10000.0] = 0.0
    hduS = fits.ImageHDU()
    hduS.data = seg2
    hduS.header = hdr
    hlistS.append(hduS)


hduI0 = fits.PrimaryHDU()
hduI0.header = f[0].header
hlistI.insert(0, hduI0)

hduW0 = fits.PrimaryHDU()
hduW0.header = f[0].header
hlistW.insert(0, hduW0)

hduIA = fits.HDUList(hlistI)
hduWA = fits.HDUList(hlistW)

hduIA.writeto(path + science_filename, overwrite=True)
hduWA.writeto(path + weight_filename, overwrite=True)


if streak_count > 0:
    hduS0 = fits.PrimaryHDU()
    hduS0.header = f[0].header
    hlistS.insert(0, hduS0)
    hduSA = fits.HDUList(hlistS)
    hduSA.writeto(path + science_filename.replace('.fits', '.lightstreak.fits'), overwrite=True)
    print('Lightstreak masked')

print('Science frame and weight image created.\n')

