import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
from astropy.io import fits
from astropy.wcs import WCS
import astropy
from astropy import units as u
import matplotlib.patches as patches

########## colours - Hogwart is my Home palette #######
#Hogwart is my Home palette
red='#b70000'
yellow='#ffc600'
blue='#000f9e'
green='#0f8400'
red='#b70000'
yellow='#ffc600'
blue='#000f9e'
green='#0f8400'

#########################################################################################

############## upload the cube & read parameters ####################0

cubename='uid___A001_X1ed_X4b0.cube_01.J0348-2749_B6.fits'
hdulist = fits.open(cubename)
data=hdulist[0].data
wcs = WCS(hdulist[0].header)

print wcs.wcs.cunit
if np.size(hdulist) == 2 :

	beams=hdulist[1].data
	### read the beam
	cell= astropy.wcs.utils.proj_plane_pixel_scales(wcs)*3600	#cell size | pix/deg -> arcsec
	beam_a = np.median(beams['BMAJ'])/cell[0] #*3600			 		#in pixel 2 x major axis 
	beam_b = np.median(beams['BMIN'])/cell[1] #*3600 					#in pixel  2 x minor axis
	beam_pa = np.median(beams['BPA'])							#for plotting

else:
	head = hdulist[0].header
	cell= astropy.wcs.utils.proj_plane_pixel_scales(wcs)#*3600	#cell size
	beam_a = head['BMAJ']/cell[0] #*3600			 		#in pixel 2 x major axis
	beam_b = head['BMIN']/cell[1] # *3600					#in pixel  2 x minor axis
	beam_pa = head['BPA']	
	print "header", beam_a/cell[0]/3600., beam_a*3600
####### cube parameters #######
x = np.shape(data)[3]	
y = np.shape(data)[2]
data_new = data*0
Z =np.shape(data)[1]

print x/2., y/2.

size=20
fig,ax = plt.subplots(1)
data_new = data
#rect = patches.Rectangle((x/2.-size/2.,y/2.-size/2.),size,size,linewidth=1,edgecolor='r',facecolor='none')
#ax.add_patch(rect)
for z in range(Z):
	for i in range(-size,size,1):
		for j in range(-size,size,1):
			data_new[0][z][int(x/2.) + i][int(y/2.)+j] = np.nan
ax.imshow(data[0][10][:][:])		
plt.plot()
plt.show()


### save fits file #####
hdu=fits.PrimaryHDU(data_new, header=hdulist[0].header)
if np.size(hdulist) == 2:
    hdu2=fits.TableHDU(beams, header=hdulist[1].header)
hdul = fits.HDUList([hdu])
hdu.writeto(cubename.replace('.fits','_mask.fits'))
hdulist.close()

