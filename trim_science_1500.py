'''

Usage:

    In script directory:
    python trim_science_1500.py path_to_file/scienceframes.txt

    In files directory:
    python script_to_file/trim_science_1500.py scienceframes.txt


This code was created to further trim the science frames with a slightly smaller radius at the time before comprehensive as a quick fix. 
Both scienceframe and weight frame have to be present.
Under normal circumstances this code wouldn't be necessary as you can simply choose a smaller optic_flat_asu1_filtermask.fits.


'''
import numpy as np
from astropy.io import fits
import sys
import os

filtermask_radius = 1500

path_list = sys.argv[1]
listname = path_list.split('/')[-1]
path = path_list[: - len(listname)]

files = []
with open(path_list, 'r') as f:
    for l in f:
        files.append(l.strip())


optic_dir = os.environ['_90PRIME_OPTIC_DIR']



for i in range(len(files)):
    f = fits.open(path + files[i])
    w = fits.open(path + files[i].replace('.fits', '.wt.fits'))

    print(files[i])

    filt = f[0].header['filter'].lower()
    month = f[0].header['date'].split('-')[1]
    if month == '01':
        print('Jan data, skip')
        pass
    elif month == '02':
        month_dir = '/Feb/'
    elif month == '03':
        month_dir = '/Mar/'

    if month != '01':
        maskfile = optic_dir + month_dir + 'optic_flat_' + filt + '_filtermask_{}.fits'.format(filtermask_radius)
        mask = fits.open(maskfile)

        hlistI = []
        hlistW = []
        for j in range(1, 17):
            fj = f[j].data
            wj = w[j].data
            hdrI = f[j].header
            hdrW = w[j].header
            mj = mask[j].data
            h = mj == 0
            fj[h] = 0
            wj[h] = 0
            hduI = fits.ImageHDU()
            hduW = fits.ImageHDU()
            hduI.data = fj
            hduW.data = wj
            hduI.header = hdrI
            hduW.header = hdrW

            hlistI.append(hduI)
            hlistW.append(hduW)

        hduI0 = fits.PrimaryHDU()
        hduW0 = fits.PrimaryHDU()
        hduI0.header = f[0].header
        hduW0.header = w[0].header
        hlistI.insert(0, hduI0)
        hlistW.insert(0, hduW0)

        hduIA = fits.HDUList(hlistI)
        hduWA = fits.HDUList(hlistW)
        hduIA.writeto(path + files[i], overwrite=True)
        hduWA.writeto(path + files[i].replace('.fits', '.wt.fits'), overwrite=True)


