import numpy as np
import matplotlib.pyplot as plt

import sys
import matplotlib as mpl
from astropy.cosmology import Planck15 as cosmo
from astropy.cosmology import z_at_value
import astropy.units as u
from astropy import constants as const
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import matplotlib.patches as patches
import numpy.ma as ma    
import matplotlib.ticker as ticker

mpl.rc ('xtick',labelsize=18)
mpl.rc ('ytick',labelsize=18)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

yellow='#ffdb00'
red='ff2500'
#################################################

def CO_Lum_prime(z,fobs, Int_flux):
	
	dL = cosmo.luminosity_distance(z).value
	Lline_prime = 3.257e7 * Int_flux *  dL**2/(1+z)**3 / fobs**2 #K km/s pc^2

	return Lline_prime


def Mmol(Lco10,alpha=4.):
	#SLED choice 0- MW (normal SF, 1-M82 (ULIRG), 2 - Lagos+12 z=2 (SIMS) 3 - Lagos+12 z=2
	#returns LCO(1-0)' (K km/s pc^2) and  the H2 mas in Msol

	Mmol = alpha * Lco10
	return Mmol

def CO10_Lum_prime(z,Int_flux):
	co = 115.27
	
	dL = cosmo.luminosity_distance(z).value
	Lline_prime = 3.257e7 * Int_flux *  dL**2/(1+z)/co**2 #K km/s pc^2

	return Lline_prime

### Poissonian limits from Gehrels et al. 2018 APJ 303
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


#################################################

######
# Take input flux + cenrtral frequency of the b=observed line and calculate the probability for the certain CO transition from Shark SAMS
######
print("################################################################################ \n Returns probability of the observed line to be a certain CO transition. \n Probabilities are based on the population of galaxies generated from Shark SAMS (Lagos et al 2018)\n for a galaxies with log(Molecular Mass) > 9.5.\n################################################################################")

## CO transitions
CO=[115.27,230.538,345.796,461.041,576.268,691.473,806.652, 921.7997,1036.912393,1151.985452] #GHz
CO_label = ["CO(1-0)","CO(2-1)","CO(3-2)","CO(4-3)","CO(5-4)","CO(6-5)","CO(7-6)", "CO(8-7)", "CO(9-8)", "CO(10-9)"]

Volume = 282823. #Mpc3 LAMACAL deep volume


#J "CO(1-0) - 0 (10) ","CO(2-1) - 1 (11) ","CO(3-2) - 2 (12)","CO(4-3) - 3 (13) ","CO(5-4) -4","CO(6-5) -5 ","CO(7-6) -6", "CO(8-7) -7", "CO(9-8) -8", "CO(10-9) -9"
# in Shark file 10, 11, 12 
cubes, names= np.loadtxt('CO_candidates_list.txt', usecols=(1,2),dtype='str',comments='#', unpack=True)
SN,Width,f0, int_flux, flux_err= np.loadtxt('CO_candidates_list.txt', usecols=(6,8,10,12,13),dtype='float', unpack=True)

####
#print "central frequency [GHz]", freq0, "Integrated flux [mJy km/s]", flux
#### Create the table of possible transitions 
transition_prob= np.zeros((np.size(CO_label), 4))
#col 0 - J, col1 - z col2 - flux col3 - probability for the transition

##load transition bins
xedges= np.load("xbins_linear.npy")
yedges=np.load("ybins_linear.npy")
nbins=np.size(xedges)

#COVol = np.loadtxt("ALMACALdeep_Volumes_CO.txt", skiprows=3)
#luminosities = np.linespace(5,13,50)

## load completness

compeltnesH =np.load("completnes.npy")
completnes_x = np.load("compeltnes_xedges.npy") #SN
completnes_y = np.load("compeltnes_yedges.npy") #width km/s

Nprob,COprob,Z=[],[],[]
M=0
S = np.size(cubes)
Nproberr= 100
transitions_array = np.zeros((S,Nproberr,10,5)) #[n_det][fluxerr][Jtrans][prob, L',z,L'(1-0),completnes]
alpha = 4.3
#relaibility coeficient
rel_coeff = 0.9 #comes from sofia | to be improved

print("Read the probability function for each candidate (including MC over the flux errors)")

for s in range(S): #loop over the candidates
	print('Candidates {} / {}'.format(s+1,S))
	flux,ferr = int_flux[s], flux_err[s] #mJy # read the flux
	freq0 = f0[s]
	
	Nprob,COprob,Z,LCO10=[],[],[],[]
	M=0

	#check completnes for detection parametes \
	sn, w  = SN[s], Width[s]
	if w < 50. : w= min(completnes_y )
	for k in range(np.size(completnes_x)-1):
		x1,x2 = completnes_x[k],completnes_x[k+1]
		if x1 <= sn <= x2: 
			x=k
			for l in range(np.size(completnes_y)-1):

				y1,y2 = completnes_y[l],completnes_y[l+1]
				if y1 <= w <= y2: y=l

	comp_coef = compeltnesH[x][y] #completnes coefficient for a given detections 1/coefficent is the actual nr of objects like that
	# include the flux error, draw fropm the error range

	# create the range of possible flux values
	flux_range = np.linspace(start=flux-ferr, stop=flux+ferr, num=1000) #1000 realisations
	#draw a random realisation - 100 times
	for zz in range(Nproberr):
		#print('Error {} / 100'.format(zz+1))
		index=np.random.randint(0,Nproberr) #draw the random index
		Flux = flux_range[index] #current flux
		Nprob,COprob,Z,LCO10=[],[],[],[]

		# create the probability function
		for j in range(np.size(CO_label)): 
			#print('Transition {} / 10'.format(j+1))
			jup=j+1
			z = CO[j]/freq0 -1.
			flu = np.log10(flux)
			# read in the histogram 
			H = np.load(CO_label[j]+"_linear_hist.npy")
			M += np.sum(H)
			if jup > 1 :
				LCO10_H = np.load("I_CO10_"+str(jup)+".npy")
			if z > 0 and z < 10:
				#check which bin you are in for k in range(nbins-1):
				for k in range(nbins-1):
					x1,x2 = xedges[k],xedges[k+1]
					if x1<=z <=x2: 
						x = k
					
						for l in range(nbins-1):
							y1,y2 = yedges[l],yedges[l+1]
							if y1<=flu <=y2: 
								
								y = l
								
								Nprob.append(H[y][x])
								COprob.append(j)
								Z.append(z)
								if jup > 1:
									Ico10 = LCO10_H[y][x] #it is a log
									lco10_prime = CO10_Lum_prime(z,Ico10) #from mJy to Jy
									LCO10.append(lco10_prime)
								else:
									
									lco10_prime = CO10_Lum_prime(z,flux/1.e3)
									LCO10.append(lco10_prime)
			elif z < 0:
				COprob.append(j+1)
				Z.append(str(0))
				Nprob.append(str(0))
				LCO10.append(str(0))
		Nprob = np.array(Nprob, dtype='float')
		N = np.sum(Nprob)	
		Mnew=np.sum(Nprob/M)
		Nprob_all = Nprob/M
		LCO10 = np.array(LCO10,dtype='float')
		for i in range(np.size(COprob)):
			j = COprob[i]
			z = CO[j]/freq0 -1.
			prob = Nprob[i]/M/Mnew
			jup = int(j) +1
			if prob > 0: 
				Lprime = CO_Lum_prime(z,freq0, flux)
				
				L10prime = LCO10[i]
			else: 
				prob , Lprime,L10prime= 0,0,0
			array=[prob,Lprime,z,L10prime,comp_coef]
			#results_array[l][n][:]
			for k in range(5):
				array[k]
				transitions_array[s][zz][i][k] = array[k]

### use transition array to draw the probable samples
ntry= 1000 # muber of tries for probability drawing
results_array=np.zeros((ntry*Nproberr,S,3)) #[ntries],[nprob][z,L'1-0,MH2, comp]
nr_detections= np.shape(transitions_array)[0]
nr_trans = np.shape(transitions_array)[2]

fid=0.9
####create redhsift bins for roH2 - 

## !!!! TEST !!!
# ASPECS - like bins egdes
z_bins = np.array([0,0.5,1.0,2.0,3.,4.])
ZB=np.size(z_bins)-1
errors_array=np.zeros((ntry,ZB,6)) # [ntries][bin] [ errors 1u,2u,3u,1l,2l,3l]
# corresponding volume: in these redhsift bins sum the volumes coming from each transition to get Volume [Mpc^3] per z bin
Volumes = np.array([1355.2,  6558.1,32788.45,42469.9,47429.59]) #bvolume per redshfi bin

# luminosity bins 0.5 dex
luminosity_bins = np.arange(8,13,1) #8-8.5, 8.5-9, 9-9.5, 9.5-10, 10-10.5, 10.5-11, 11-11.5, 11.5-12, 12-12.5

LFa1, LFa2, LFa3, LFa4, LFa5=np.zeros((1,np.size(luminosity_bins)-1)), np.zeros((1,np.size(luminosity_bins)-1)),np.zeros((1,np.size(luminosity_bins)-1)),np.zeros((1,np.size(luminosity_bins)-1)),np.zeros((1,np.size(luminosity_bins)-1))
LFa1err, LFa2err, LFa3err, LFa4err, LFa5err=np.zeros((6,np.size(luminosity_bins)-1)), np.zeros((6,np.size(luminosity_bins)-1)),np.zeros((6,np.size(luminosity_bins)-1)),np.zeros((6,np.size(luminosity_bins)-1)),np.zeros((6,np.size(luminosity_bins)-1))
MH1err, MH2err, MH3err, MH4err, MH5err=np.zeros((6,np.size(z_bins)-1)), np.zeros((6,np.size(z_bins)-1)),np.zeros((6,np.size(z_bins)-1)),np.zeros((6,np.size(z_bins)-1)),np.zeros((6,np.size(z_bins)-1))

Lum1,Lum2,Lum3,Lum4,Lum5 = [],[],[],[],[]
mh1,mh2,mh3,mh4,mh5 = [],[],[],[],[]
LFall = np.zeros((1,np.size(luminosity_bins)-1))
print("Calculate the  CO luminosity function and H2 masses including the probability function for each detection")


for l in range(ntry): #realsations of the sample
	#print('Realisations {} / {}'.format(l+1,ntry*Nproberr))
	lco10_sample, zsmample=[],[]
	for zz in range(Nproberr):
		for n in range(nr_detections):
						
			P,prob = [],[]
			#for each detection

			if names[n] == 'J0334-4008':
				for m in range(nr_trans):
					# 1) create a table of transitions probabilities and table of indexes
					if m==0:
						P = np.append(P,1.)
					else: 
						P = np.append(P,0.)

			else:
				for m in range(nr_trans):
					# 1) create a table of transitions probabilities and table of indexes
					
					P = np.append(P,transitions_array[n][zz][m][0])

			Indexes = np.array(range(10))
			# 2) draw a transition to ba assigned to the detection, according to the probabilities
			#draw a transition index I
			I = np.random.choice(Indexes,size=1,p=P)[0]

			for k in range(3):
				results_array[l*Nproberr+zz][n][k] = transitions_array[n][zz][I][k+2] #[z,L'(1-0), completnes]
Z1,Z2,Z3,Z4,Z5=[],[],[],[],[]	
for l in range(ntry*Nproberr): 
	#### L'CO(1-0) Redshift bins
	l1,l2,l3,l4,l5=np.array([[0,0]]),np.array([[0,0]]),np.array([[0,0]]),np.array([[0,0]]),np.array([[0,0]])

	z1,z2  = z_bins[0],z_bins[1]
	for n in range(nr_detections):
			if  z_bins[0] <= results_array[l][n][0] < z_bins[1]: #
				Z1 = np.vstack(Z1,results_array[l][n][0])
				l1 = np.vstack((l1, [np.log10(results_array[l][n][1]),results_array[l][n][2]]))
			elif z_bins[1]  <= results_array[l][n][0] < z_bins[2]:
				l2 = np.vstack((l2, [np.log10(results_array[l][n][1]),results_array[l][n][2]]))
				Z2 = np.vstack(Z2,results_array[l][n][0])
			elif z_bins[2]  <= results_array[l][n][0] < z_bins[3]:
				l3 = np.vstack((l3, [np.log10(results_array[l][n][1]),results_array[l][n][2]]))
				Z3 = np.vstack(Z3,results_array[l][n][0])			
			elif z_bins[3]  <= results_array[l][n][0] < z_bins[4]:
				l4 = np.vstack((l4, [np.log10(results_array[l][n][1]),results_array[l][n][2]]))		
				Z4 = np.vstack(Z4,results_array[l][n][0])	
			elif z_bins[4]  <= results_array[l][n][0] < z_bins[5]:
				l5 = np.vstack((l5, [np.log10(results_array[l][n][1]),results_array[l][n][2]]))	
				Z5 = np.vstack(Z5,results_array[l][n][0])

	#### Bin the luminocities in each z bin 
	## secondary bin the molecular masses
	## zbin 1 0.0 - 0.5

	if np.shape(l1)[0] > 1:
		lum,lum_raw=0,0
		ncounts= np.zeros((1,np.size(luminosity_bins)-1))
		for j in range(np.size(luminosity_bins)-1): # binning in luminosity
			for i in range(np.size(l1.T[0])):
				if l1.T[0][i] > 0:
					if luminosity_bins[j] <= l1.T[0][i] <  luminosity_bins[j+1]:
						ncounts.T[j] += fid/l1.T[1][i]
						lum += 10**l1.T[0][i]/l1.T[1][i]
				
		LF1 = ncounts*fid/Volumes[0]
		LFa1 = np.vstack((LFa1,LF1))

		ncounts_ma = np.ma.masked_where(ncounts == 0, ncounts)

		#poissonina errors
		errors=np.zeros((6,np.size(luminosity_bins)-1))
		for t in range(3):
			pl=Poiss_lower_limit(ncounts_ma,t+1)
			pu=Poiss_upper_limit(ncounts_ma,t+1)
			errors[t],errors[t+3]=pl*fid/Volumes[0],pu*fid/Volumes[0]
		LFa1err = np.dstack((LFa1err,errors))

		## molecular masses
		mh = lum*alpha*fid/Volumes[0]
		if mh > 0:
			mh1.append(mh)	
		mh1_n = np.size(mh1)

		errors_mh=np.zeros((6,np.size(z_bins)-1))
		for t in range(6):
		 	for l, ll in zip(range(np.size(luminosity_bins)-1), [8.5,9.5,10.5,11.5]):
		 		errors_mh[t]+= alpha*errors[t][l]*10**ll
		# 	pl=Poiss_lower_limit(mh1_n,t+1)
		# 	pu=Poiss_upper_limit(mh1_n,t+1)
			
		# 	errors_mh[t],errors_mh[t+3]=pl*lum*fid/Volumes[0]/mh1_n,pu*lum*alpha*fid/Volumes[0]/mh1_n
		MH1err = np.dstack((MH1err,errors_mh))

		## zbin 1 0.5 - 1.0
	if np.shape(l2)[0] > 1:
		lum, lum_raw=0,0
		ncounts= np.zeros((1,np.size(luminosity_bins)-1))
		for j in range(np.size(luminosity_bins)-1): # binning in luminosity
			for i in range(np.size(l2.T[0])):
				if l2.T[0][i] > 0:
					if luminosity_bins[j] <= l2.T[0][i] <  luminosity_bins[j+1]:
						ncounts.T[j] += fid/l2.T[1][i]
						lum += 10**l2.T[0][i]/l2.T[1][i]
						lum_raw+=10**l2.T[0][i]
		LF2 = ncounts*fid/Volumes[1]
		LFa2 = np.vstack((LFa2,LF2))
		ncounts_ma = np.ma.masked_where(ncounts == 0, ncounts)

		#poissonina errors
		errors=np.zeros((6,np.size(luminosity_bins)-1))
		for t in range(3):
			pl=Poiss_lower_limit(ncounts_ma,t+1)
			pu=Poiss_upper_limit(ncounts_ma,t+1)
			errors[t],errors[t+3]=pl*fid/Volumes[1],pu*fid/Volumes[1]
		LFa2err = np.dstack((LFa2err,errors))

		## molecular masses
		mh=lum*alpha*fid/Volumes[1]
		if mh > 0:
			mh2.append(mh)	
		mh2_n = np.size(mh2)
		errors_mh=np.zeros((6,np.size(z_bins)-1))
		for t in range(6):
		 	for l, ll in zip(range(np.size(luminosity_bins)-1), [8.5,9.5,10.5,11.5]):
		 		errors_mh[t]+= alpha*errors[t][l]*10**ll	## molecular mass errors
		# errors_mh=np.zeros((6,np.size(z_bins)-1))
		# for t in range(3):
		# 	pl=Poiss_lower_limit(mh2_n,t+1)
		# 	pu=Poiss_upper_limit(mh2_n,t+1)
		# 	errors_mh[t],errors_mh[t+3]=pl*lum*alpha*fid/Volumes[1]/mh2_n,pu*lum*alpha*fid/Volumes[1]/mh2_n
		MH2err = np.dstack((MH2err,errors_mh))

		## zbin 1 1.0 - 2.0
	if np.shape(l3)[0] > 1:
		lum, lum_raw=0,0
		ncounts= np.zeros((1,np.size(luminosity_bins)-1))
		for j in range(np.size(luminosity_bins)-1): # binning in luminosity
			for i in range(np.size(l3.T[0])):
				if l3.T[0][i] > 0:
					if luminosity_bins[j] <= l3.T[0][i] <  luminosity_bins[j+1]:
						ncounts.T[j] += fid/l3.T[1][i]
						lum += 10**l3.T[0][i]/l3.T[1][i]
						lum_raw+=10**l3.T[0][i]
		LF3 = ncounts*fid/Volumes[2]
		LFa3 = np.vstack((LFa3,LF3))
		ncounts_ma = np.ma.masked_where(ncounts == 0, ncounts)

		#poissonina errors
		errors=np.zeros((6,np.size(luminosity_bins)-1))
		for t in range(3):
			pl=Poiss_lower_limit(ncounts_ma,t+1)
			pu=Poiss_upper_limit(ncounts_ma,t+1)
			errors[t],errors[t+3]=pl*fid/Volumes[2],pu*fid/Volumes[2]

		LFa3err = np.dstack((LFa3err,errors))

		## molecular masses
		mh=lum*alpha*fid/Volumes[2]
		if mh > 0:
			mh3.append(mh)	
		mh3_n = np.size(mh3)
		## molecular mass errors
		errors_mh=np.zeros((6,np.size(z_bins)-1))
		for t in range(6):
		 	for l, ll in zip(range(np.size(luminosity_bins)-1), [8.5,9.5,10.5,11.5]):
		 		errors_mh[t]+= alpha*errors[t][l]*10**ll
		# for t in range(3):
		# 	pl=Poiss_lower_limit(mh3_n,t+1)
		# 	pu=Poiss_upper_limit(mh3_n,t+1)
		# 	errors_mh[t],errors_mh[t+3]=pl*lum*alpha*fid/Volumes[2]/mh3_n,pu*lum*alpha*fid/Volumes[2]/mh3_n
		MH3err = np.dstack((MH3err,errors_mh))

		## zbin  2.0 - 3.0
	if np.shape(l4)[0] > 1:
		lum, lum_raw=0,0
		ncounts= np.zeros((1,np.size(luminosity_bins)-1))
		for j in range(np.size(luminosity_bins)-1): # binning in luminosity
			for i in range(np.size(l4.T[0])):
				if l4.T[0][i] > 0:
					if luminosity_bins[j] <= l4.T[0][i] <  luminosity_bins[j+1]:
						ncounts.T[j] += fid/l4.T[1][i]
						lum += 10**l4.T[0][i]/l4.T[1][i]
						lum_raw+=10**l4.T[0][i]
		LF4 = ncounts*fid/Volumes[3]
		LFa4 = np.vstack((LFa4,LF4))
		#poissonina errors
		ncounts_ma = np.ma.masked_where(ncounts == 0, ncounts)
		errors=np.zeros((6,np.size(luminosity_bins)-1))
		for t in range(3):


			pl=Poiss_lower_limit(ncounts_ma,t+1)
			pu=Poiss_upper_limit(ncounts_ma,t+1)
			errors[t],errors[t+3]=pl*fid/Volumes[3],pu*fid/Volumes[3]

		LFa4err = np.dstack((LFa4err,errors))

		
		## molecular masses
		mh=lum*alpha*fid/Volumes[3]
		if mh > 0:
			mh4.append(mh)	
		mh4_n = np.size(mh4)
		## molecular mass errors
		errors_mh=np.zeros((6,np.size(z_bins)-1))
		for t in range(6):
		 	for l, ll in zip(range(np.size(luminosity_bins)-1), [8.5,9.5,10.5,11.5]):
		 		errors_mh[t]+= alpha*errors[t][l]*10**ll
		# 	pl=Poiss_lower_limit(mh4_n,t+1)
		# 	pu=Poiss_upper_limit(mh4_n,t+1)
		# 	errors_mh[t],errors_mh[t+3]=pl*lum*alpha*fid/Volumes[3]/mh4_n,pu*lum*alpha*fid/Volumes[3]/mh4_n
		MH4err = np.dstack((MH4err,errors_mh))

		## zbin  3.0 - 4.0
	if np.shape(l5)[0] > 1:
		lum,lum_raw=0,0
		ncounts= np.zeros((1,np.size(luminosity_bins)-1))
		for j in range(np.size(luminosity_bins)-1): # binning in luminosity
			for i in range(np.size(l5.T[0])):
				if l5.T[0][i] > 0:
					if luminosity_bins[j] <= l5.T[0][i] <  luminosity_bins[j+1]:
						ncounts.T[j] += fid/l5.T[1][i]
						lum += 10**l5.T[0][i]/l5.T[1][i]
						lum_raw+=10**l5.T[0][i]
		LF5 = ncounts*fid/Volumes[4]
		LFa5 = np.vstack((LFa5,LF5))
		ncounts_ma = np.ma.masked_where(ncounts == 0, ncounts)

		#poissonina errors
		errors=np.zeros((6,np.size(luminosity_bins)-1))
		for t in range(3):
			pl=Poiss_lower_limit(ncounts_ma,t+1)
			pu=Poiss_upper_limit(ncounts_ma,t+1)
			errors[t],errors[t+3]=pl*fid/Volumes[4],pu*fid/Volumes[4]
		LFa5err = np.dstack((LFa5err,errors))

		## molecular masses
		mh=lum*alpha*fid/Volumes[4]
		if mh > 0:
			mh5.append(mh)	
		mh5_n = np.size(mh5)
		## molecular mass errors
		errors_mh=np.zeros((6,np.size(z_bins)-1))
		for t in range(6):
			for l, ll in zip(range(np.size(luminosity_bins)-1), [8.5,9.5,10.5,11.5]):
		 		errors_mh[t]+= alpha*errors[t][l]*10**ll
		# 	pl=Poiss_lower_limit(mh5_n,t+1)
		# 	pu=Poiss_upper_limit(mh5_n,t+1)
		# 	errors_mh[t],errors_mh[t+3]=pl*lum*alpha*fid/Volumes[4]/mh5_n,pu*lum*alpha*fid/Volumes[4]/mh5_n
		MH5err = np.dstack((MH5err,errors_mh))

print(np.log10([np.mean(mh1),np.mean(mh2),np.mean(mh3),np.mean(mh4),np.mean(mh5)]))

##### save the parameters for the plots

for lf, lferr,mh, mherr, i in zip([LFa1,LFa2,LFa3,LFa4,LFa5], [LFa1err,LFa2err,LFa3err,LFa4err,LFa5err],[mh1,mh2,mh3,mh4,mh5],[MH1err, MH2err,MH3err, MH4err, MH5err],range(1,6,1),):
	np.save("LF_{}".format(i),lf)
	np.save("LF_err_{}".format(i),lferr)
	np.save("MH2_{}".format(i),mh)
	np.save("MH2_err_{}".format(i),mherr)


fig2 = plt.figure(2, figsize=(12,6))
ax2 = plt.axes(yscale='log')
X2,Y2 = np.array([]),np.array([])
for i, mherr in zip(range(5), [MH1err, MH2err,MH3err, MH4err, MH5err]):
		X2 = np.append(X2,np.mean([z_bins[i], z_bins[i+1]]))
		err1l,err1u = np.mean(mherr[0]),np.mean(mherr[3])
		err2l,err2u= np.mean(mherr[1]),np.mean(mherr[4])
		err3l,err3u = np.mean(mherr[2]),np.mean(mherr[5])
		print(err3l,err3u,err2l,err2u,np.log10(err1l),np.log10(err1u))
		rect2 = patches.Rectangle((z_bins[i],err2l),z_bins[i+1]-z_bins[i],err2u-err2l,linewidth=1,edgecolor='r',facecolor='none') #2 sigma
		ax2.add_patch(rect2)
		rect1 = patches.Rectangle((z_bins[i],err3l),z_bins[i+1]-z_bins[i],err3u-err3l,linewidth=1,edgecolor='b',facecolor='none') #1 sigma
		ax2.add_patch(rect1)
		#rect3 = patches.Rectangle((z_bins[i],err1l),z_bins[i+1]-z_bins[i],err1u-err1l,linewidth=1,edgecolor='b',facecolor='none') #1 sigma
		#ax2.add_patch(rect3)


Y2=[np.mean(mh1),np.mean(mh2),np.mean(mh3),np.mean(mh4),np.mean(mh5)]
ax2.plot(X2,Y2,'ro')
ax2.set_ylabel(r" $\rho$(M$_{\rm H2}$) [M$_{\odot}$ Mpc$^{\rm -3}$]", fontsize=15)
ax2.set_xlabel("Redshift", fontsize=15)
ax2.xaxis.set_ticks_position('both')
ax2.yaxis.set_ticks_position('both')
ax2.set_ylim([1.5*10**6,2*10**9])

## aspecs
zaspecs = [np.mean([0.003,0.369]), np.mean([1.006,1.738]), np.mean([2.008, 3.107]), np.mean([3.011,4.475])]
Zerr_as=[zaspecs[0]-0.003,zaspecs[1]-1.006,zaspecs[2]-2.008,zaspecs[3]-3.011]
MH_aspecs1=[np.mean([5.89,6.80]), np.mean([7.74,7.96]), np.mean([7.50,7.96]), np.mean([7.20,7.62])]
MH_err_as1=[[MH_aspecs1[0]-5.89,MH_aspecs1[1]-7.74,MH_aspecs1[2]-7.50,MH_aspecs1[3]-7.20 ],[6.80-MH_aspecs1[0],7.96-MH_aspecs1[1],7.96-MH_aspecs1[2],7.62-MH_aspecs1[3] ]]
MH_aspecs2=[np.mean([5.40,7.01]), np.mean([7.63,8.05]), np.mean([7.26,8.10]), np.mean([6.97,7.77])]
MH_err_as2=[[MH_aspecs2[0]-5.40,MH_aspecs2[1]-7.63,MH_aspecs2[2]-7.26,MH_aspecs1[3]-6.97 ],[7.01-MH_aspecs1[0],8.05-MH_aspecs1[1],8.10-MH_aspecs1[2],7.77-MH_aspecs1[3] ]]
for z, mh1,mh2 in zip(zaspecs,MH_aspecs1,MH_aspecs2):

		rect2 = patches.Rectangle((z_bins[i],err2l),z_bins[i+1]-z_bins[i],err2u-err2l,linewidth=1,edgecolor='r',facecolor='none') #2 sigma
		ax2.add_patch(rect2)
		rect1 = patches.Rectangle((z_bins[i],err3l),z_bins[i+1]-z_bins[i],err3u-err3l,linewidth=1,edgecolor='b',facecolor='none') #1 sigma
		ax2.add_patch(rect1)
		#rect3 = patches.Rectangle((z_bins[i],err1l),z_bins[i+1]-z_bins[i],err1u-err1l,linewidth=1,edgecolor='b',facecolor='none') #1 sigma
		#ax2.add_patch(rect3)

# plt.errorbar(zaspecs,MH_aspecs2, xerr=Zerr_as, yerr=MH_err_as2)
# plt.errorbar(zaspecs,MH_aspecs1, xerr=Zerr_as, yerr=MH_err_as1)
# fig3=plt.figure(3)
# plt.hist(l1.T[1],color='r', histtype='step',)
# plt.hist(l2.T[1], color='g',histtype='step')
# plt.hist(l3.T[1], color='b',histtype='step')
# plt.hist(l4.T[1], color='y', histtype='step')
# plt.hist(l5.T[1], color='k', histtype='step')
# fig3.savefig("bins.png")
plt.show()

