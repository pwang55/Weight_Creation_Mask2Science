from numpy import *
from astropy.io import fits
from glob import glob
import sys
import os


coef = 0.97             # correlation coefficient for finding streak
nofpixel = 50           # number of pixels criteria for finding streak

infile = sys.argv[1]	# bgrms file
segfile = sys.argv[2]	# segmentation file
f = fits.open(infile)
seg = fits.open(segfile)



# Get optic mask directory environment variable
optic_dir = os.environ['_90PRIME_OPTIC_DIR']
original_file = infile.replace('.bgrms.fits','.fits')
file_hdr = fits.getheader(original_file)    # get header from original fits file associated with the bgrms & seg fits
filemonth = file_hdr['date'].split('-')[1]  # get the month of file, either 02 or 03
filt = file_hdr['FILTER'].lower()           # get the filter of file
if filemonth == '02':
    month_dir = '/Feb/'
elif filemonth == '03':
    month_dir = '/Mar/'

ffiltermask = optic_dir+month_dir+'optic_flat_'+filt+'_filtermask.fits'
filtermask = fits.open(ffiltermask)


# Get the absolute directory path of where this python script is so that defect_seg_mask.fits can be found
script_path = sys.path[0]
defectmask = fits.open(script_path+'/../large_files/defect_segs_mask.fits') # For now put it outside the repository cause GitHub won't allow large files

print('\n')
print('Converting: '+'\t'+infile.split('/')[-1]+'\t'+'Using Weight Mask:'+'\t'+'optic_flat_'+filt+'_filtermask.fits')


nx = f[1].data.shape[0]
ny = f[1].data.shape[1]
dat = zeros((16, nx, ny), dtype=float)
streak = zeros((16, nx, ny), dtype=float)

hlist = []
slist = []

for j in range(1,17):
	streakcheck = []
	im = f[j].data
	s = seg[j].data

	dat[j-1] = 1.0/im**2
	streak[j-1] = s

	mask = filtermask[j].data
	h = mask > 0.0
	streak[j-1][~h] = 0.0       # Mask seg outside filter so later we only look for streaks inside filtermask

	hs = streak[j-1] > 0        # extracted sources selection, 2-D array of True/False
	xs = array(nonzero(hs))[0]  # x coordinates of all nonzero pixels in the seg fits, 1-D array
	ys = array(nonzero(hs))[1]  # y coordinates of all nonzero pixels in the seg fits, 1-D array
	value = streak[j-1][hs]     # segmentation value of extracted pixels, 1-D array, same order as xs, ys, such that xs[i], ys[i] has value value[i]

        # loop over all possible segmentation value
	for k in range(int(min(value)), 1+int(max(value))):
		idx = where(value == k)[0]                              # get the index numbers for the pixels that have value=k, 1-D, this index is for xs, ys and value
		if len(idx) > 0:                                        # if k actually exist, then check if it form a line
			a, b = corrcoef(xs[idx], ys[idx])               # use xs, ys that have value=k to see if they form a line by calculating correlation coefficient
			if (fabs(a[1]) > coef) and (len(idx)>nofpixel): # C_ij of covariant matrix is the correlation
				dat[j-1][xs[idx], ys[idx]] = 0.0        # set light streak = 0
				streakcheck.append(1.0)                 # a checker for later, the code will know if a streak is found
			else:
				streak[j-1][xs[idx], ys[idx]] = 0.0


	hduS = fits.ImageHDU()
	hduS.data = streak[j-1]
	hduS.header = seg[j].header
	slist.append(hduS)


	dh = defectmask[j].data > 0.0   # defect mask
	dat[j-1][dh] = 0.0		# defect = 0
	hduI = fits.ImageHDU()
	hduI.data = dat[j-1]
	hduI.header = f[j].header
	hlist.append(hduI)

hdu0 = fits.PrimaryHDU()
hdu0.header = f[0].header
hlist.insert(0,hdu0)
hduA = fits.HDUList(hlist)
hduA.writeto(infile.replace('.bgrms.fits','.wt.fits'))

if len(streakcheck) > 0:
	hduS0 = fits.PrimaryHDU()
	hduS0.header = seg[0].header
	slist.insert(0, hduS0)
	hduSA = fits.HDUList(slist)
	streakmaskname = infile.replace('.bgrms.fits','.lightstreakmask.fits')
	hduSA.writeto(streakmaskname)
	print('Light streak has been found, mask created:')
	print(streakmaskname)

print('\n')
print('Weight map created:'+'\t'+infile.replace('.bgrms.fits','.wt.fits').split('/')[-1])
print('\n')

f.close()
filtermask.close()
defectmask.close()
seg.close()
