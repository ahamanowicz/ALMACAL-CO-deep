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
import scipy.integrate as integrate
from scipy.ndimage.interpolation import rotate
from astropy import units as u
import sys
import astropy
import os
import glob
### signal shape #####

## gaussian ##

def gauss(x, mu, fwhm,A):
	
	"""
	mu - peak frequency
	sig - width (FWHM) / [ 2 * sqrt(2*ln(2)) ]
	A - peak flux

	"""
	sig = fwhm/(2*np.sqrt(2*np.log(2)))
	return A*np.exp(-np.power((x - mu)/sig, 2.)/2)

def boxy(x,mu,wid,A):
	"""
	mu - central frequency
	wid -  width of full box ( borders mu - wid/2 ; mu + wid/2)
	A - peak flux
	"""
	box = np.array([])
	for f in x:
	
		if f < mu-wid/2 or f > mu+wid/2:
			box=np.append(box,0)
		elif f >= mu-wid/2 and f<= mu+wid/2:
			box=np.append(box,A)
	return box

def horn(x,mu,wid,A):
	"""
	mu - central frequency
	wid -  width of full box ( borders mu - wid/2 ; mu + wid/2)
	A - peak flux
	used x^4 polynomial
	"""
	po = 4
	a = A/3./(wid/2.)**po
	b= 1.*A/3.
	
	box = np.array([])
	for f in x:
	
		if f < mu-wid/2 or f > mu+wid/2:
			box=np.append(box,0)
		elif f >= mu-wid/2 and f<= mu+wid/2:
			box=np.append(box,abs((f-mu)**po)*a+b)
	return box
 
# 1) get the cube rms
def get_cube_rms(cube, directory, plot=True, save=True):

    # returns the rms of the cube (rms) and the array of outlier channels (outl), can print out and save the pdf plot
   
    c = directory+cube
    print(cube)
    #open the fits cube
    hdulist = fits.open(c)
    data=hdulist[0].data

    #whole cube rms
    rms = np.sqrt(np.mean(np.square(data[0]))) 
    print(rms)
    
    #rms by channel, flag outliers
    Z =np.shape(data)[1]
    rms_chan=np.array([],dtype=float)
    for z in range(Z):	

        rms_single = np.sqrt(np.mean(np.square(data[0][z][:][:])))
        rms_chan = np.append(rms_chan,rms_single)

    med = np.median(rms_chan)
    std = np.std(rms_chan)
    outl=np.array([], dtype='int')

    #choose outliers
    for z in range(Z):
        if abs(rms_chan[z]-med) > std:

            outl = np.append(outl, int(z))
            
    print("outliers:", outl)
    
    #plot rms to pdf file
    if plot == True:
        fig=plt.figure(1,figsize=(6,6))  
        plt.plot(rms_chan, 'bo')
        plt.axhline(med,c='r',linewidth=2)
        plt.axhline(med+std,linestyle='--',c='r',linewidth=2)
        plt.axhline(med-std,linestyle='--',c='r',linewidth=2)
        plt.xlim(-2,Z+2)
        plt.xlabel('channel')
        plt.ylabel('rms')
        plt.title(cube.replace(".fits", ""))
        plt.tight_layout()
        plt.show()
        if save == True:
            fig.savefig(cube.replace('.fits','_rms.pdf'))

        
    hdulist.close()
    
    return rms, outl

#mask the flux =0 

def mask_zero(cube, dire, save=False):
    #returns the corrected array
    
    c = dire+cube
    print(cube)
    
    #open the cube
    hdulist = fits.open(c)
    data = hdulist[0].data
    wcs = WCS(hdulist[0].header)
    
    data_new = data*0
    Z =np.shape(data)[1]
    
    for z in range(Z):

        a = data[0][z][:][:]
        a = np.ma.array(a, mask=np.isnan(a))
        flux = np.sum(a)  
        #m,ask the flux=0
        if flux == 0.0:

            data_new[0][z][:][:] = np.nan
        else:
            data_new[0][z][:][:] = data[0][z][:][:]

    if save == True:    
        hdu=fits.PrimaryHDU(data_new, header=hdulist[0].header)
        hdu.writeto(c.replace('.fits', '_nozero.fits'))
        hdulist.close()
    
    return data_new

### insert mock sources

def get_hdu(cube, dire):
    c = dire + cube
    
    #open the fits
    hdulist = fits.open(c)
    data=hdulist[0].data
    
    return hdulist

def read_cube_params(cube, dire, rms):
    c = dire + cube
    rms_cube=rms
    
    #open the fits
    hdulist = fits.open(c)
    data=hdulist[0].data
    wcs = WCS(hdulist[0].header)
    
    if np.size(hdulist) == 2 :

        beams=hdulist[1].data
        ### read the beam
        cell= astropy.wcs.utils.proj_plane_pixel_scales(wcs)	#cell size
        beam_a = np.median(beams['BMAJ'])/cell[0]			 		#in pixel 2 x major axis
        beam_b = np.median(beams['BMIN'])/cell[1] 					#in pixel  2 x minor axis
        beam_pa = np.median(beams['BPA'])							#for plotting

    else:
        head = hdulist[0].header
        cell= astropy.wcs.utils.proj_plane_pixel_scales(wcs)	#cell size
        beam_a = head['BMAJ']/cell[0]			 		#in pixel 2 x major axis
        beam_b = head['BMIN']/cell[1] 					#in pixel  2 x minor axis
        beam_pa = head['BPA']				
    
    print("beam", beam_a, beam_b, beam_pa, cell[0])
    
    ####### cube parameters #######

    ### XYZ size ###

    # edges of the cube - from cube formation X = Y, imsize is always defiend as the square
    x1 = y1 = 0
    x2 = np.shape(data)[3]
    y2 = np.shape(data)[2]
    pix =  np.shape(data)[3]						# nuber of pixels in one axis (same for X and Y in this cube formation pattern )
    X = Y = range(x1,x2)							# number of pixels in each direction
    chn = np.shape(data)[1]							# number of channels , Z

    Z = range(1,chn)								# frequency slices - channels
    
    #print("X,Y,Z",x2,y2,chn)
    
    ### frequency coverage ###

    freq_1 = hdulist[0].header['CRVAL3'] 			# Hz | edge of the cube in frequency
    df = hdulist[0].header['CDELT3']		
    freq_2 =(freq_1 + (chn-1)*df) 					# Hz | edge of the cube in frequency

    freq_1 = freq_1/1.e9							# convert to GHz
    freq_2 = freq_2/1.e9							# convert to GHz

    if freq_2 > freq_1:

        freq = np.linspace(freq_1,freq_2,chn) 		# GHz 

    else: 

        freq = np.linspace(freq_2,freq_1,chn) 		# GHz 

    chan_width = abs(freq_1 - freq_2)/chn 			# width of a channel in GHz

    #print("chan_width: "+str(chan_width) + " GHz")
    #print("freq ["+str(freq_1)+" ; "+str(freq_2)+ "] GHz")

    #print("noise sigma: "+str(round(sigma_n*1e3,3))+ " mJy/beam")

    g = open(cube.replace(".fits",".stats"), "w")
    g.write("cubename"+'\t'+ cube+"\n")
    g.write("rms"+'\t'+ str(rms_cube)+"\n")
    g.write("cell [arcsec, arcsec]"+'\t'+ str(cell[0:2])+"\n")
    g.write("beam_axes [pix, pix]"+'\t'+ str(beam_a)+" "+str(beam_b) +"\n")
    g.write("beam_PA [deg]"+'\t'+ str(beam_pa) +"\n")
    g.write("image_size [pix pix]"+'\t'+ str(x2)+" "+str(y2) +"\n")
    g.write("channels"+'\t'+ str(chn) +"\n")
    g.write("frequency_range [GHz]"+'\t'+ '\t'+ str(freq_1)+" "+str(freq_2) +"\n")
    g.write("channel_width [GHz]"+'\t'+ str(chan_width) +"\n")
    g.close()
    
    hdulist.close()
    
    return cube, rms_cube, cell[0:2], beam_a, beam_b, beam_pa, x2,y2,chn, freq_1,freq_2, chan_width

############## upload the cube & read parameters ####################
def mock_insert(cube, rms, cube_params, hdulist, n_mock=20, no=0):
    
    ################ signal parameters ##################

    #n_mock = 20 # how many mock signals plugged - 20 is plenty for a single cube

    #detection array

    ### create the list of mock signals 

    Xs,Ys,Zs,Freqs,FPs,SNs,Ws,Shs=[],[],[],[],[],[],[],[] # saving moc signal parameters

    #cube name 
    c = cube_params[0]
    f = open(c.replace('.fits','_mock_'+str(no)+'.txt'), 'w')
    f.write("OBJ_ID" + '\t' + "X" + '\t' + "Y" + '\t' + "Z" + '\t' + "Frequency" + '\t' + "F_peak" + '\t' + "S/N" + '\t' + " width" + '\t' +'\t' + "n_chann" + '\t'+"shape" + '\t'+ "F_integrated" + '\n')	
    f.write("[]" + '\t' + "[pix]" + '\t' + "[pix]" + '\t' + "[pix]" + '\t' + "[Hz]" + '\t'+'\t' + "[mJy/beam]"  + '\t' +'[]'+'\t'+ " [km/s]"  + '\t'+"[]" +'\t'+"[]" + '\t'+ "[mJy/beam]" + '\n')	



    #positions
    x1,y1 = 0,0
    x2,y2 = int(cube_params[6]), int(cube_params[7])
    Xposition = range(x1+3,x2-3)				# XY position space for random (not centralized on any pixel)
    Yposition = range(y1+3,y2-3)				# XY position space for random (not centralized on any pixel)
    # 3 pix frame included - ommitting the region
    Xs,Ys,Zs,Freqs,FPs,SNs,Ws,Shs=[],[],[],[],[],[],[],[] # saving moc signal parameters
    detection= np.zeros(np.shape(cube))				# detection array - same size as the cube

    freq_1,freq_2, chn = float(cube_params[9]),float(cube_params[10]), int(cube_params[8])

    sigma_n=float(cube_params[1])
    chan_width = float(cube_params[11])
    
    X = Y = range(x1,x2)							# number of pixels in each direction
    Z = range(1,chn)
    
    
    if freq_2 > freq_1:

        freq = np.linspace(freq_1,freq_2,chn) 		# GHz 

    else: 

        freq = np.linspace(freq_2,freq_1,chn) 		# GHz 
        
    beam_a, beam_b, beam_pa = float(cube_params[3]), float(cube_params[4]), float(cube_params[5])

    print(beam_a, beam_b, beam_pa)
          
    for N in range(n_mock):

        ### position XYZ ###

        ## mean freq ## - soaced with the channel width

        freq_space = np.linspace(freq_1,freq_2,chn) 	# freq space for random, spaced by the channel number, central frequency must be in some channel
        #choose the chanel for central freq
        Z0 = np.random.choice(range(0,chn))
        #find corresponding frequncy
        f0 = freq_space[Z0] 				# central frequency
        #print("mean freq: "+str(f0)+ " GHz" )

        ### F_peak ###
        sig = np.arange(1,8.1,0.1) # sigma space for random; step every 0.1
        A = np.random.choice(sig) # choose sigma level of signal
        F_peak = A*sigma_n 								# Jy/beam 

        #print("F peak = "+str(round(F_peak*1e3,3))+" S/N = " + str(round(A,3)))

        ### width and integrated flux (ideal) ###
        
        #width in number of channels, but from the km/s. The step 5 km/s

        #km/s # width [km/s] to width [GHz]   
        #df/f0 = dv/c => df = dv/c * f0

        width = np.arange(50,801,5) 				# [km/s] velocity space for random
        width0 = np.random.choice(width)				# choose the width [km/s]
        fwidth0 = width0*1.e3/const.c *f0 				# conver [km/s] to [GHz]

        n_chan = fwidth0/chan_width 						#widdth of detection in channels

        #print("width = "+str(round(width0,3))+" km/s","N channels: ", int(n_chan))

        ## shape of signal ##

        #0 - gaussian
        #1 - boxy
        #2 - 2-horn

        shape = [0,1,2]									# shape space for random
        #shape0 = np.random.choice(shape)				# choose the signal shape
        #fix gaussian
        shape0=0
        
        if shape0 == 0:
            signal = gauss(freq,f0,fwidth0,F_peak)		# model 
            F_int,err = integrate.quad(gauss, freq[0], freq[-1],args=(f0,fwidth0,F_peak))
            #print "shape: gaussian"
        elif shape0 == 1:
            signal = boxy(freq,f0,fwidth0,F_peak)		# model
            F_int = np.trapz(signal,freq)
            #print "shape: box"
        else:
            signal = horn(freq,f0,fwidth0,F_peak)		# model 
            F_int = np.trapz(signal,freq)
            #print "shape: double horn"


        ### XY position ###
        #choosing the pixel representing unresolved source

        X0=np.random.choice(Xposition)						# choose x0 					
        Y0=np.random.choice(Yposition)						# choose y0
        #print X0
            #check if the position isn't already taken +/- 10 pix around

        if N > 0:
            for n in range(N):
                X1=int(Xs[n])
                Y1 =int(Ys[n])
                if np.linalg.norm(np.array([X1,Y1])-np.array([X0,Y0])) < 10:
                    while np.linalg.norm(np.array([X1,Y1])-np.array([X0,Y0])) < 10:

                        X0=np.random.choice(Xposition)						# choose x0 					
                        Y0=np.random.choice(Yposition)	

        #plug-in detection - delta function

        #####. delta funcion representing unresolved detection #####
        for z in Z:
            detection[0][z][Y0][X0] = signal[z]


        Xs.append(str(X0))
        Ys.append(str(Y0))
        Zs.append(str(Z0))
        Freqs.append(str(f0*1e9))
        FPs.append(str(round(F_peak*1e3,5)))
        SNs.append(str(round(A,3)))
        Ws.append(str(round(width0,3)))
        Shs.append(str(shape0))

        f.write(str(N+1) + '\t' + str(X0) + '\t' + str(Y0) + '\t' + str(Z0) + '\t' + str(f0*1e9) + '\t' + str(round(F_peak*1e3,5)) + '\t' + str(round(A,3)) + '\t'+ str(round(width0,3)) + '\t'+ '\t' + str(int(n_chan))+'\t'  + str(shape0) +'\t'+str(round(F_int*1e3,3))+'\n')

    ### beam convolution ####
    kernel = np.outer(gaussian(x2, beam_a/2.35), gaussian(y2, beam_b/2.35))		#kernel FWHM to sigma
    kernel = rotate(kernel,beam_pa/2.) #rotation of the beam according to the PA
    data_new = cube*0
    
    convolved =  np.zeros(np.shape(cube))
    for z in Z:


            Slice = detection[0][z]
            blurred = fftconvolve( Slice,kernel, mode='same')
            convolved[0][z]=blurred

    data_new = cube + convolved

    f.close()
    
    ### save fits file #####
    hdu=fits.PrimaryHDU(data_new, header=hdulist[0].header)
#    if np.size(hdulist) == 2:
#            beams=hdulist[1].data
#            hdu2=fits.TableHDU(beams, header=hdulist[1].header)
    hdul = fits.HDUList([hdu])
    hdu.writeto(c.replace('.fits','_mock_'+str(no)+'.fits'), overwrite=True)
    hdulist.close()
    
    return data_new

def make_flagged_array(cube_params, outl):

    imsize=cube_params[6]
    #prepae the flagged array format
    print(outl)
    flagged_array=[]
    edges = [outl[0]]
    for i in range(1,np.size(outl)):
        if (int(outl[i]) - int(outl[i-1]) )!= 1:
            edges.append(outl[i-1])
            edges.append(outl[i])

    edges.append(outl[-1])

    ## create flagged array
    print(edges)
    for i in range(0,np.size(edges),2):
        ar= [0,imsize,0,imsize,int(edges[i]),int(edges[i+1])+1]
        flagged_array.append(ar)

    print("flagged aray",flagged_array)
    return flagged_array

def sofia_files(c, N, flagged_array):

    for no in range(N):

        #change parameter file to include the file name
        fin=open('/Users/ahamanow/Desktop/ALMA/CO_search/ALMA_fields/template.param','r') # template

        par = c.replace('.fits','_mock_'+str(no)+'.fits')
        fout=open(par.replace(".fits",".sofia"),'wt') # new file for sofia
        #print(par)
        for line in fin:
            a=line.split()
            if "import.inFile" in a:
                fout.write(line.replace(a[2], par))
            elif "flag.regions" in a:
                fout.write(line.replace(a[2], str(flagged_array)))
            else:
                fout.write(line)
        fin.close()
        fout.close()

        os.system("/Users/ahamanow/SoFiA/SoFiA-master/sofia_pipeline.py "+ par.replace(".fits",".sofia"))
    
    files_in_dir = glob.glob('*.fits')
    for file in files_in_dir:
        os.remove(file)
    
