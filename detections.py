import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import fits
#import aplpy
from matplotlib.patches import Circle, Ellipse
from astropy.wcs import WCS
from astropy import units as u
from mpl_toolkits.axes_grid.anchored_artists import AnchoredAuxTransformBox
import astropy
import sys
from specutils import Spectrum1D,SpectralRegion
from specutils.analysis import line_flux, gaussian_fwhm, centroid, gaussian_sigma_width, fwhm
from specutils.fitting import fit_lines
from spectral_cube import SpectralCube
import matplotlib as mpl

mpl.rc ('xtick',labelsize=14)
mpl.rc ('ytick',labelsize=14)
mpl.rc('font',family='serif', size=20)
########## colours - Hogwart is my Home palette #######
#Hogwart is my Home palette
red='#b70000'
yellow='#ffc600'
blue='#000f9e'
green='#0f8400'

#Fruit salad	
beam='#cccccc'
contour=yellow#'#baffb0'
circ=red#'#ff6363'\
d=20.
cubes, clas= np.loadtxt('candidates_list.txt', usecols=(1,12),dtype='str',comments='#', unpack=True)
X,Y,Z,Width, SN, w50, f0, rms = np.loadtxt('candidates_list.txt', usecols=(3,4,5,7,6,11,10,14),dtype='float', unpack=True)
#print X,Y,Z
#print clas
N = np.size(cubes)
chan_width =  0.0155639648438 #GHz
vel_lim=1000.

fig=plt.figure(1,figsize=(12,8))

for i in range(N):
	plt.subplot(4,4,i+1)
	cube = cubes[i]+".fits"
	spec_cube = SpectralCube.read(cube)

	print(cube)
	hdul = fits.open(cube) #open original cube
	data = hdul[0].data
	wcs = WCS(hdul[0].header)
	wcs=wcs.dropaxis(2)
	wcs=wcs.dropaxis(2)

	if np.size(hdul) > 1.:
		beams=hdul[1].data

		### define beam to plot
		cell= astropy.wcs.utils.proj_plane_pixel_scales(wcs)*3600
		beam_a = np.mean(beams['BMAJ'])/cell[0]
		beam_b = np.mean(beams['BMIN'])/cell[1]
		beam_pa = np.mean(beams['BPA'])
		print("beam = ",beam_a, beam_b, beam_pa)
		print(np.mean(beams['BMAJ']), np.mean(beams['BMIN']))
	else: 
		beam_a, beam_b,beam_pa=0.,0.,0.

	#create zer-moment map - collapse the cube over channels with detection
	Zsize = np.shape(data)[1]
	collapsed=data[0,0]			

	for z in range(int(Z[i] - Width[i]/2.),int(Z[i] +Width[i]/2.)):
		if z < Zsize and z > 0.:
			collapsed += data[0,z]		

	#read the sectrum form the peak XY position from original cube
	#read the spectrum from the cube
	subcube_pix = spec_cube[:,int(Y[i]),int(X[i])]
	wav_pix,flux_pix = subcube_pix.spectral_axis, subcube_pix #flux [Jy/beam]
	spec_pix= Spectrum1D(spectral_axis=wav_pix.to(u.GHz), flux=flux_pix) #frequency spec 
	cube_vel_pix = subcube_pix.with_spectral_unit(u.km/u.s,velocity_convention='radio',rest_value=f0[i]*u.GHz ) 
	vel_pix = cube_vel_pix.spectral_axis
	spec_vel_pix = Spectrum1D(spectral_axis=vel_pix, flux=flux_pix)
	## spectrum stats
	med=np.median(flux_pix.value)
	RMS=rms[i]
	std=np.std(flux_pix.value)
	mean = np.mean(flux_pix.value)	

	plt.step(vel_pix,flux_pix.value*1e3, c='k', where='mid')
	plt.xlim([-vel_lim,vel_lim])
	# ax1.step(vel_pix,flux_pix.value*1e3, c='k', where='mid')
	# ax1.fill_between(vel_pix,flux_pix.value*1e3, step='mid', color='yellow')

	# secax = ax1.twiny()
	# secax.step(wav_pix,flux_pix.value*1e3, c='k')
	# f1,f2 = freq_from_vel(f0[i], -vel_lim),freq_from_vel(f0[i], vel_lim)
	# secax.set_xlim([f1,f2])
	# #ax1.axvline(Z[i], linestyle='--', c=blue)
	# #ax1.axhline(med, c=yellow, linestyle='--')
	# #ax1.axhline(med+3*std, c=green, linestyle='--')
	# #ax1.axhline(med-3*std, c=green, linestyle='--')
	# ax1.set_xlabel('velocity [km/s]')
	# secax.set_xlabel('frequency [GHz]', labelpad=10)
	# ax1.set_ylabel('flux [mJy/beam]')
	# ax1.set_xlim([-vel_lim,vel_lim])
	# #ax1.set_xlim(Z[i] - 3*Width[i],Z[i] +3*Width[i])
	# #ax1.set_title(cube+ " "+str(SN[i])+" "+clas[i])
	# #ax1.set_ylim(med - 3*std, med+6*std)
	# #ax1.axvline((Z[i]+w50[i]),c=red, linestyle='--')
	# #ax1.axvline((Z[i]-w50[i]),c=red, linestyle='--')
	# #### fits image ######

	# # plot the zoomed fits + beam

	# box = AnchoredAuxTransformBox(ax2.transData, loc=3,frameon=False,)
	# el = Ellipse((0,0), width=beam_b, height=beam_a, angle=beam_pa,fc="none",ec=beam,lw=.7,alpha=0.8) # in data coordinates! hatch='///'
	# box.drawing_area.add_artist(el)
	# ra = ax2.coords[0]
	# ra.set_major_formatter('hh:mm:ss.ss')
	# #plt.tick_labels.set_yformat('ddmmss.ss')

	# ax2.imshow(collapsed, cmap='gray') 
	# sigma_list = np.array([med+4*std,med+5*std,med+6*std,med+7*std, med+8*std,med+9*std, med+10*std,med+12*std])
	# ax2.contour(collapsed, colors=contour, levels=sigma_list, alpha=0.5) #levels=np.logspace(-4.7, -3., 10),
	# ax2.add_artist(box)
	# ax2.set_xlabel('RA')
	# ax2.set_ylabel('DEC')
	# ax2.set_xlim(X[i]-d,X[i]+d)
	# ax2.set_ylim(Y[i]-d,Y[i]+d)

	# #ax2.set_title('integration = '+str(round(time,3))+ ' s')
	# lon = ax2.coords[0]
	# lat = ax2.coords[1]
	# lon.set_ticks(exclude_overlapping=True)
	# lat.set_ticks(exclude_overlapping=True)
	# lat.set_ticklabel_position('r')
	# lat.set_axislabel_position('r')
	# #rad = max([abs(x_peak[n]-xmin[n]),abs(x_peak[n]-xmax[n]),abs(y_peak[n]-ymin[n]),abs(y_peak[n]-ymax[n])])
	# #ax2.plot(X[i],Y[i],'+',color=red,ms=7)
	# #ax2.add_patch(circ)
	# #print X[i],Y[i]

	# #ax2.set_title("RA = "+ra[n]+" DEC = "+dec[n]+"\n S/N = "+str(sn[n]) )

	# ## ax3 - full spectrum
	# ax3.axis('off')
	# #ax3.step(wav_pix.to(u.GHz),flux_pix.value*1e3,'-',c='k',where='mid') #real spectrum extracted form the cube

	# #ax3.axvline(Z[i], linestyle='--', c=blue)
	# #ax3.axhline(med, c=yellow, linestyle='--')
	# #ax3.axhline(med+3*std, c=green, linestyle='--')
	# #ax3.axvline((Z[i]+w50[i]),c=red, linestyle='--')
	# #ax3.axvline((Z[i]-w50[i]),c=red, linestyle='--')
	# #ax3.set_xlabel('channel')
	# #ax3.set_ylabel('flux')
	# #ax3.set_ylim(med - 3*std, med+6*std)
	# #plt.tight_layout()

#draw the beam

plt.show()

