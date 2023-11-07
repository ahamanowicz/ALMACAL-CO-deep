import numpy as np
import matplotlib.pyplot as plt
import sys
import matplotlib as mpl

#################################################
#Code parameters
mpl.rc ('xtick',labelsize=18)
mpl.rc ('ytick',labelsize=18)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
######
# Take input flux + cenrtral frequency of the b=observed line and calculate the probability for the certain CO transition from Shark SAMS
######p
print("################################################################################ \n Returns probability of the observed line to be a certain CO transition. \n Probabilities are based on the population of galaxies generated from Shark SAMS (Lagos et al 2018)\n for a galaxies with log(Molecular Mass) > 9.5.\n################################################################################")

freq0 =float(raw_input("Observed frequency of the line (central in GHz): "))
flux = np.log10(float(raw_input("Observed integrated flux (in mJy km/s): ")))

#print M
## CO transitions
CO=[115.27,230.538,345.796,461.041,576.268,691.473,806.652, 921.7997,1036.912393,1151.985452] #GHz
CO_label = ["CO(1-0)","CO(2-1)","CO(3-2)","CO(4-3)","CO(5-4)","CO(6-5)","CO(7-6)", "CO(8-7)", "CO(9-8)", "CO(10-9)"]
#J "CO(1-0) - 0 (10) ","CO(2-1) - 1 (11) ","CO(3-2) - 2 (12)","CO(4-3) - 3 (13) ","CO(5-4) -4","CO(6-5) -5 ","CO(7-6) -6", "CO(8-7) -7", "CO(9-8) -8", "CO(10-9) -9"
# in Shark file 10, 11, 12 

print("central frequency [GHz]", freq0, "Integrated flux [mJy km/s]", flux)

#### Create the table of possible transitions 
transition_prob= np.zeros((np.size(CO_label), 4))
#col 0 - J, col1 - z col2 - flux col3 - probability for the transition

##load transition bins
xedges= np.load("xbins_new.npy")
yedges=np.load("ybins_new.npy")
nbins=np.size(xedges)

Nprob,COprob,Z, SLED=[],[],[],[]
M=0

for j in range(np.size(CO_label)):

	z = CO[j]/freq0 -1.
	flu = np.log10(flux)
	# read in the histogram 
	H = np.load(CO_label[j]+"_new_hist.npy")
	M += np.sum(H)

	#print z, flu
	if z > 0 and z < 10:
		#print z
		#check which bin you are in for k in range(nbins-1):
		for k in range(nbins-1):
			#print xedges[k],xedges[k+1], yedges[l],yedges[l+1],H.T[k][l]
			x1,x2 = xedges[k],xedges[k+1]
			if x1<=z <=x2: 
				x = k+1
			
				for l in range(nbins-1):
					y1,y2 = yedges[l],yedges[l+1]
					print(y2-y1, x2-x1)
					if y1<=flux <=y2: 
						#print y1,y2,flux,l+1
						y = l+1
						#plt.plot(x,y,'or',ms=10)
						#plt.plot(y,x,'ob',ms=10)
						Nprob.append(H[y][x])
						COprob.append(j)
						Z.append(z)
						#print H.T[y][x], j 

Nprob = np.array(Nprob, dtype='float')
N = np.sum(Nprob)	
print(Nprob, M)
Mnew=np.sum(Nprob/M)
Nprob_all = Nprob/M
for i in range(np.size(COprob)):
	j = COprob[i]
	z = CO[j]/freq0 -1.
	prob = Nprob[i]/M/Mnew
	print("transition", CO_label[j], "z =", np.round(z,3), "probability", np.round(prob*100,3))
	COprob[i] +=1
fig=plt.figure(1, figsize=(7,7))

#ASPECS

plt.plot(COprob,Nprob/M/Mnew, "o--", c='navy',ms=10)
plt.xlabel("CO transition Jup",fontsize=18)
plt.ylabel("Probability",fontsize=18)
plt.ylim([-0.05, 1.05])
plt.xlim([-1,11])
plt.show()


