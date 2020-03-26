from numpy import *
from astropy.io import fits
import sys
"""

outdated code just for reference only

"""



file0 = sys.argv[1]
coef = 0.97

seg = fits.open(file0)
#print 'fine'

hlist = []
for j in range(1,17):

	im = seg[j].data
	h = im > 0	# extracted sources
	x = array(nonzero(h))[0]
	y = array(nonzero(h))[1]
	value = im[h]
	
	for k in range(int(min(value)), int(max(value))):
		idx = where(value == k)[0]
		if len(idx) == 0:
			1*1
		else:
			a, b = corrcoef(x[idx], y[idx])
			if fabs(a[1]) < coef:
				im[x[idx], y[idx]] = 0

	hduI = fits.ImageHDU()
	hduI.data = im
	hduI.header = seg[j].header
	hlist.append(hduI)

hdu0 = fits.PrimaryHDU()
hdu0.header = seg[0].header
hlist.insert(0,hdu0)
hduA = fits.HDUList(hlist)
hduA.writeto('test.fits')


seg.close()

