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

    To change backsize, use option: backsize=64 or -backsize=64

This script takes corrected fits image as input, use photutils to calculate background RMS map, convert it to weight map, 
and mask corrected image and weight image with mask, output to science frames. Lightstreaks will be detected and masked in the weight image.


"""
import numpy as np
from astropy.io import fits
from astropy.stats import SigmaClip
from photutils import Background2D, MedianBackground
from photutils.segmentation import detect_sources
from photutils import deblend_sources
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from photutils.segmentation import SourceCatalog
import sys
import os
import time
import multiprocessing as mp

if len(sys.argv) == 1:
    print(__doc__)
    sys.exit()

t1 = time.time()

backsize = 32
backsize_lightstreak = 32
filtersize = 7
filtersize_lightstreak = 7
streak_sigma_threshold = 2.5
elongation_threshold = 10.0
npixels = 30
sigma_clip = SigmaClip(sigma=3.)
bkg_estimator=MedianBackground()
filtermask_radius = 1500

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

        elif arg.split('=')[0] == 'backsize' or arg.split('=')[0] == '-backsize':
            backsize = int(arg.split('=')[1])

backsize_lightstreak = backsize

science_filename=filename.replace('corrected_', 'science_')
weight_filename=science_filename.replace('.fits', '.wt.fits')

optic_dir=os.environ['_90PRIME_OPTIC_DIR']
filt=f[0].header['filter'].lower()
month=f[0].header['date'].split('-')[1]
if month == '02':
    month_dir = '/Feb/'
elif month == '03':
    month_dir = '/Mar/'
# elif month == '01':

if month == '02' or month == '03':
    maskfile = optic_dir + month_dir + 'optic_flat_' + filt + '_filtermask_{}.fits'.format(filtermask_radius)
    mask = fits.open(maskfile)
    print('Create weight map for: ' + filename)
    print('Filtermask: ' + month_dir + 'optic_flat_' + filt + '_filtermask_{}.fits'.format(filtermask_radius))

elif month == '01':
    print('Create weight map for: ' + filename)
    print('January data, no mask needed.')

hlistI = []
hlistW = []
# hlistS = []
# streak_count = 0


def jan_function():

    f = fits.open(file0)

    dat = f[0].data
    hdr = f[0].header

    # Use a large mesh size to detect lightstreak and extremely saturated stars
    bkg_lightstreak = Background2D(dat, box_size=backsize_lightstreak, filter_size=filtersize_lightstreak, sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    threshold = bkg_lightstreak.background + (streak_sigma_threshold * bkg_lightstreak.background_rms)
    sigma = 4.0 * gaussian_fwhm_to_sigma
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    kernel.normalize()
    segm = detect_sources(dat, threshold, npixels=npixels, kernel=kernel)
    # segm_deblend = deblend_sources(dat, segm, npixels=npixels, kernel=kernel, nlevels=32, contrast=0.001)
    cat = SourceCatalog(dat, segm)
    # cat = source_properties(dat, segm_deblend)
    # tab = cat.to_table()
    elongation = cat.elongation.value
    streak_id = cat.label[elongation > elongation_threshold]
    segd = segm.data
    # if len(streak_id) > 0:
    #     for i in range(len(streak_id)):
    #         # Fill lightstreak and extremely saturated stars with background value
    #         dat[segd == streak_id[i]] = bkg_lightstreak.background[segd == streak_id[i]]

    # Use a smaller mesh size to estimate reasonable background
    # bkg = Background2D(dat, box_size=backsize, filter_size=filtersize, sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    bkg = bkg_lightstreak
    wt = 1 / bkg.background_rms ** 2

    if len(streak_id) > 0:
        for i in range(len(streak_id)):
            wt[segd == streak_id[i]] = 0.0

    wt[dat == 0] = 0

    hduI = fits.PrimaryHDU()
    hduW = fits.PrimaryHDU()
    hduI.data = dat
    hduW.data = wt
    hduI.header = hdr
    hduW.header = hdr
    hduI.header['BZERO'] = 0.0
    hduW.header['BZERO'] = 0.0

    return hduI, hduW





def mpfunction(j):

    # hlistIj = []
    # hlistWj = []
    f = fits.open(file0)
    mask = fits.open(maskfile)

    dat = f[j].data
    hdr = f[j].header
    maskj = mask[j].data
    dat[maskj == 0] = 0

    # Use a large mesh size to detect lightstreak and extremely saturated stars
    bkg_lightstreak = Background2D(dat, box_size=backsize_lightstreak, filter_size=filtersize_lightstreak, sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    threshold = bkg_lightstreak.background + (streak_sigma_threshold * bkg_lightstreak.background_rms)
    sigma = 4.0 * gaussian_fwhm_to_sigma
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    kernel.normalize()
    segm = detect_sources(dat, threshold, npixels=npixels, kernel=kernel)
    # segm_deblend = deblend_sources(dat, segm, npixels=npixels, kernel=kernel, nlevels=32, contrast=0.001)
    cat = SourceCatalog(dat, segm)
    # cat = source_properties(dat, segm_deblend)
    # tab = cat.to_table()
    elongation = cat.elongation.value
    streak_id = cat.label[elongation > elongation_threshold]
    segd = segm.data
    # if len(streak_id) > 0:
    #     for i in range(len(streak_id)):
    #         # Fill lightstreak and extremely saturated stars with background value
    #         dat[segd == streak_id[i]] = bkg_lightstreak.background[segd == streak_id[i]]

    # Use a smaller mesh size to estimate reasonable background
    # bkg = Background2D(dat, box_size=backsize, filter_size=filtersize, sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    bkg = bkg_lightstreak
    wt = 1 / bkg.background_rms ** 2
    wt[maskj == 0] = 0

    if len(streak_id) > 0:
        for i in range(len(streak_id)):
            wt[segd == streak_id[i]] = 0.0

    hduI = fits.ImageHDU()
    hduW = fits.ImageHDU()
    hduI.data = dat
    hduW.data = wt
    hduI.header = hdr
    hduW.header = hdr
    hduI.header['BZERO'] = 0.0
    hduW.header['BZERO'] = 0.0
    # hlistIj.append(hduI)
    # hlistWj.append(hduW)

    return hduI, hduW

    # seg2[seg2 < 10000.0] = 0.0
    # hduS = fits.ImageHDU()
    # hduS.data = seg2
    # hduS.header = hdr
    # hlistS.append(hduS)


if month == '01':
    hduI0, hduW0 = jan_function()
    hlistI.append(hduI0)
    hlistW.append(hduW0)
    hduIA = fits.HDUList(hlistI)
    hduWA = fits.HDUList(hlistW)

    hduIA.writeto(path + science_filename, overwrite=True)
    hduWA.writeto(path + weight_filename, overwrite=True)



if month == '02' or month == '03':
    p1 = mp.get_context("fork").Pool(8)
    p2 = mp.get_context("fork").Pool(8)
    # jlist = [x for x in range(1, 17)]
    jlist1 = [x for x in range(1, 9)]
    jlist2 = [x for x in range(9, 17)]

    with p1:
        out1 = p1.map(mpfunction, jlist1)

    with p2:
        out2 = p2.map(mpfunction, jlist2)

    # hlistI = [out[x][0] for x in range(16)]
    # hlistW = [out[x][1] for x in range(16)]

    hlistI1 = [out1[x][0] for x in range(8)]
    hlistW1 = [out1[x][1] for x in range(8)]
    hlistI2 = [out2[x][0] for x in range(8)]
    hlistW2 = [out2[x][1] for x in range(8)]

    hlistI = hlistI1 + hlistI2
    hlistW = hlistW1 + hlistW2


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

    # if streak_count > 0:
        # hduS0 = fits.PrimaryHDU()
        # hduS0.header = f[0].header
        # hlistS.insert(0, hduS0)
        # hduSA = fits.HDUList(hlistS)
        # hduSA.writeto(path + science_filename.replace('.fits', '.lightstreak.fits'), overwrite=True)
        # print('Lightstreak masked')

t2 = time.time()
print('Science frame and weight image created, t={}\n'.format(t2-t1))

