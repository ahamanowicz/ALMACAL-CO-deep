import matplotlib.pyplot as plt
import numpy as np
from scipy import constants as const
from scipy.stats import skewnorm
from scipy.integrate import quad
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Circle, Ellipse
from scipy.signal import gaussian, fftconvolve
from astropy.io import fits
from astropy.wcs import WCS
import astropy
from scipy.ndimage.interpolation import rotate
import scipy.integrate as integrate
import sys
import os
from scipy.special import erf
from scipy.optimize import curve_fit, least_squares

import matplotlib as mpl

def Poiss_upper_limit(n, S):
    lambda_up = n + S * np.sqrt(n+1) + (S**2+2)/3.
    return lambda_up

def Poiss_lower_limit(n, S):
    lambda_down= n-S*np.sqrt(n)+ (S**2-1)/3.
    #Beta, Gamma =[0,0.062,0.222], [0,-2.19,-1.88]
    # if n == 0: gamma=0
    # else: gamma=Gamma[S-1]
    # lambda_down = n*(1 - 1/9*n - S /3*np.sqrt(n) + Beta[S-1]*n**gamma)**3
    return lambda_down

def error_func(x,C,s,a):
	#from walter+16
	f = 0.5 * erf((x - C) /s) +a
	return f

mpl.rc ('xtick',labelsize=15)
mpl.rc ('ytick',labelsize=15)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

red='#d11141'
yellow='#ffc425'
blue='#00aedb'
green='#00b159'
orange='#f37735'
#fname = sys.argv[1]
#cubes= np.loadtxt(fname, usecols=(0,),dtype='str')
matches_array=[]
fig = plt.figure(1,figsize=(8,8))

sn_p,w_p = [],[]
sn_n, w_n=[],[]
#Zmax = 700 km
chanw = 0.0155029296875 #GHz
positive_cats = np.loadtxt("positive/pos.in", dtype='str')
negative_cats = np.loadtxt("negative/neg.in", dtype='str')

for pref in positive_cats:

	stats = "../"+pref.replace("_nozero_cat.ascii", ".stats")
	with open(stats) as f:

		content = f.read().splitlines()
	a,b = content[3].split('\t')[1].split(' ')
	a,b =float(a),float(b)
	imsize = int(content[5].split()[3])
	f1,f2 = float(content[7].split()[2]), float(content[7].split()[3])
	Zmax = float(content[6].split()[1])
	f0 = (f1+f2)/2.
	dZmax = 700*1.e3/const.c *f0 #in GHz
	dZmax = dZmax * Zmax / abs(f2-f1)
	data=np.loadtxt("positive/"+pref,comments='#',dtype='float', usecols=(15,25,29,36, 5,6,7))

	M =  np.shape(data)[0]

	if np.size(data[0]) >1:

		for m in range(M):

			#print(data[m][1], data[m][2],data[m][3])

			sn = float(data[m][1])/float(data[m][2]) # f_peak / rms
			freq = data[m][3]/1.e9 #GHz
			width = int(data[m][0])*chanw * 3.e5/freq
			#print(width)
			X,Y,Z = float(data[m][4]),float(data[m][5]),float(data[m][6])
			pb = 1-np.sqrt((X - imsize/2.)**2	+ (Y - imsize/2.)**2	)/imsize/np.sqrt(2)
			
			if int(data[m][0]) < Zmax/4. and pb  >= 0.7:

				w_p = np.append(w_p, width)
				sn_p = np.append(sn_p, sn)
				if sn >= 7:
					print("high SN",pref, X,Y,Z, sn)
			

	else:
		

		
		sn = float(data[1])/float(data[2]) # f_peak / rms
		freq = data[3]/1.e9 #GHz
		width = int(data[0])*chanw * 3e5/freq

		X,Y,Z = float(data[4]),float(data[5]),float(data[6])
		pb = 1-np.sqrt((X - imsize/2.)**2	+ (Y - imsize/2.)**2	)/imsize/np.sqrt(2)

		if int(data[0]) < Zmax/4. and pb >= 0.7:

			w_p = np.append(w_p, width)
			sn_p = np.append(sn_p, sn)
			if sn >= 7:
					print("high SN",pref, X,Y,Z, sn)

for pref in negative_cats:
	stats = "../"+pref.replace("_nozero_cat.ascii", ".stats")
	with open(stats) as f:

		content = f.read().splitlines()
	a,b = content[3].split('\t')[1].split(' ')
	a,b =float(a),float(b)
	imsize = int(content[5].split()[3])
	f1,f2 = float(content[7].split()[2]), float(content[7].split()[3])
	Zmax = float(content[6].split()[1])
	f0 = (f1+f2)/2.
	dZmax = 700*1.e3/const.c *f0 #in GHz
	dZmax = dZmax * Zmax / abs(f2-f1)



	data=np.loadtxt("negative/"+pref,comments='#',dtype='float',  usecols=(15,25,29,36,5,6,7))

	M =  np.shape(data)[0]

	if np.size(data[0]) > 1:

		for m in range(M):
			if data[m][2] > 0:

				
				sn = float(data[m][1])/float(data[m][2]) # f_peak / rms
				freq = data[m][3]/1.e9 #GHz
				width = int(data[m][0])*chanw * 3e5/freq
		
				X,Y,Z = float(data[m][4]),float(data[m][5]),float(data[m][6])
				pb = 1-np.sqrt((X - imsize/2.)**2	+ (Y - imsize/2.)**2	)/imsize/np.sqrt(2)
				if int(data[m][0]) < Zmax/4. and pb  >= 0.7:
				
					w_n = np.append(w_n, width)
					sn_n = np.append(sn_n, sn)
					if sn >= 7:
						print("neg-high SN",pref, X,Y,Z, sn)
		
	else:

		if data[2] > 0:
			sn = float(data[1])/float(data[2]) # f_peak / rms
			freq = data[3]/1.e9 #GHz
			width = int(data[0])*chanw * 3e5/freq

			X,Y,Z = float(data[4]),float(data[5]),float(data[6])
			pb = 1-np.sqrt((X - imsize/2.)**2	+ (Y - imsize/2.)**2	)/imsize/np.sqrt(2)

			if int(data[0]) < Zmax/4. and pb >= 0.7:

			
				w_n = np.append(w_n, width)
				sn_n = np.append(sn_n, sn)
				if sn >= 7:
					print("neg-high SN",pref, X,Y,Z, sn)

fig=plt.figure(1)
plt.subplot(211)
po=plt.hist(sn_p, bins=np.arange(0,8,.5), histtype='step',log=True, label='positive')
ne=plt.hist(sn_n, bins=np.arange(0,8,.5),histtype='step', log=True,label='negative')
plt.xlabel("S/N", fontsize=15)
plt.ylabel("log N", fontsize=15)
plt.legend(fontsize=15)
plt.xlim([3,8.2])
#plt.xlim([3,8.2	])
print("reliability",1-ne[0]/po[0])
plt.subplot(212)
print(po[1])
X = po[1][2:]
Y = 1- ne[0][1:]/po[0][1:]
X = [  1. , 1.5, 2.,  2.5, 3. , 3.5, 4.,  4.5, 5.,  5.5 ,6. , 6.5 , 7.5] 
Y = np.array([  0.14285714,  0.140625  ,  0.11660777,  0.06006006,
        0.22648752,  0.22613566,  0.23576251,  0.29811996,  0.30481283,
        0.25714286,  0.65957447,  0.61538462,  1.        ])

popt,pcov=curve_fit(error_func, X,Y, p0=[6.0,2.0,0.7] )
perr=np.sqrt(np.diag(pcov))

print("rel", error_func(X, popt[0],popt[1],popt[2]))
print("err", error_func(X, popt[0],popt[1],popt[2])- Y)

a



plt.plot(X,Y, 'o', color='navy', label='data')
plt.plot(np.linspace(1,8,200), error_func(np.linspace(1,8,200), popt[0],popt[1],popt[2]), color='red', label='model')
# plt.hist(w_p, bins=25, histtype='step', log=True)
# plt.hist(w_n, bins=25,histtype='step', log=True)
# plt.xlabel("width [km/s]",  fontsize=15)
plt.ylabel(r"rel = 1 - $\frac{N_{neg}}{N_{pos}}$", fontsize=15)
# plt.axvline(3,lw=3,c=yellow)
# plt.xlim([0.8, 8.2])
# plt.ylim([0,500])
# plt.grid("on")
# plt.xlabel("S/N", fontsize=20)
# plt.ylabel("width [km/s]", fontsize=20)
plt.xlim([3,8.2])

plt.legend(fontsize=15)
# plt.subplot(313)
# hist_sn_p, bins_sn_p = np.histogram(sn_p, bins=20, range=(0,8))
# hist_sn_n, bins_sn_n = np.histogram(sn_n, bins=bins_sn_p,range=(0,8))

# hfrac = 1-np.array(hist_sn_n, dtype='float')/np.array(hist_sn_p, dtype='float')

# plt.bar(bins_sn_p[:-1], hfrac, color='none', edgecolor='navy', lw=3, align='edge', width=0.4)
# plt.plot(np.linspace(1,8,200), error_func(np.linspace(1,8,200), popt[0], popt[1], popt[2]), color='red')
plt.xlabel("S/N", fontsize=15)

# plt.xlim([3,8])
fig.savefig("reliability.png")


fig2 = plt.figure(2,figsize=(10,10))

Hn,xedges,yedges = np.histogram2d(sn_n,w_n,bins=[10,10], range=[[0,8],[10,500]])
H_p,xedges_p,yedges_p = np.histogram2d(sn_p,w_p,bins=[xedges,yedges],range=[[0,8],[10,500]])
H_frac = 1 - Hn/H_p
grid = H_frac.T

X,Y=np.meshgrid(xedges,yedges)
im=plt.pcolormesh(X,Y,grid)
cb=plt.colorbar(orientation='horizontal')
cb.set_label(label="reliability", fontsize=18)
plt.clim(0,1)
dets=[5.9, 5.9, 6.4, 5.1, 4.1, 4.4, 4.9, 6.7]
detw=[35,255,62,35,31,29,61,198]
plt.scatter(dets, detw, c='red')
plt.xlabel("S/N", fontsize=18)
plt.ylabel("width [km/s]", fontsize=18)
plt.xlim([2.5,8])
#N,M = range(np.size(w_a)), np.size(sn_a)
np.save("completnes",H_frac)
np.save("compeltnes_xedges", xedges)
np.save("compeltnes_yedges", yedges)

fig2.savefig("2dreliability.png")
print("C", popt[0], "s", popt[1], "a", popt[2])

candidates=[5.926,  5.14, 4.445, 4.891, 6.655]
for cand in candidates:
	print(cand, error_func(cand, popt[0],popt[1],popt[2]))

allcand = np.loadtxt("sofia_cand.txt", usecols=(1,))
names= np.loadtxt("sofia_cand.txt", usecols=(0,), dtype='str')

rel = error_func(allcand,  popt[0],popt[1],popt[2] )
plt.figure(3)
plt.hist(rel, bins=15)
for i  in range(np.size(names)):
	rel = error_func(allcand[i], popt[0],popt[1],popt[2] )
	#if rel > 0.3:
		#print(names[i],allcand[i],round(rel,2))
print(H_frac[5][0], H_frac[6][0], H_frac[6][1], H_frac[7][0], H_frac[7][1], H_frac[7][4],H_frac[8][3],H_frac[8][3], H_frac[7][3])





# cb_ax = fig2.add_axes([0.83, 0.38, 0.02, 0.52])
# cbar = fig2.colorbar(im, cax=cb_ax)
# cbar.set_label(label='Fraction recovered', size=20)
# cbar.set_clim(0,1.0)
# ax2 = plt.subplot2grid((3,3),(2,1), colspan=2)

# hist_sn_m, bins_sn_m = np.histogram(sn_m, bins=10)
# hist_sn_det, bins_sn_det = np.histogram(sn_a, bins=bins_sn_m)
# bin_width = bins_sn_m[1] - bins_sn_m[0]

# fraction=np.array(hist_sn_det,dtype='float')/np.array(hist_sn_m, dtype='float')

# error= fraction-pl
# error=np.vstack((error,pu-fraction))
# error_sn=error
# np.save("sn_fraction.npy",fraction)
# np.save("sn_error.npy",error_sn)
# print(pu-fraction, pl-fraction)
# plt.errorbar(bins_sn_m[1:]- bin_width/2.,fraction,yerr=error,fmt='o--', c='navy', capsize=5.)
# plt.ylim([0.0,1.0])
# plt.xlim([xedges[0], xedges[-1]])
# plt.xlabel("S/N", fontsize=20)
# plt.ylabel("Detection fraction", fontsize=20)
# ax3 = plt.subplot2grid((3,3),(0,0), rowspan=2)

# hist_w_m, bins_w_m = np.histogram(w_m, bins=10)
# hist_w_det, bins_w_det = np.histogram(w_a, bins=bins_w_m)
# bin_width = bins_w_m[1] - bins_w_m[0]
# fractionw=np.array(hist_w_det,dtype='float')/np.array(hist_w_m, dtype='float')
# plw=Poiss_lower_limit(hist_w_det, 3)/np.array(hist_w_m, dtype='float')
# puw = Poiss_upper_limit(hist_w_det,3)/np.array(hist_w_m, dtype='float')
# error= fractionw-plw
# error=np.vstack((error,puw-fractionw))
# error_w = error
# np.save("w_fraction.npy",fractionw)
# np.save("w_error.npy",error_w)
# plt.errorbar(fraction,bins_w_m[1:]- bin_width/2.,xerr=error,fmt='o--', c='navy', capsize=5.)
# plt.xlim([1.0,0.0])
# plt.ylabel("width [km/s]", fontsize=20)
# plt.xlabel("Detection fraction", fontsize=20)
# plt.ylim([yedges[0], yedges[-1]])

# #plt.tight_layout()
# fig.savefig("Mock_detected_all.pdf")
# fig2.savefig("Detection_fraction2.pdf")
plt.show()
