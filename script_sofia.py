#!/bin/bash
import numpy as np
import os
import shutil

Fidelity=[0.3,0.5,0.7,0.8,0.95]
SN = [2.5,3,4,5,6,7]

for sn in SN:
	os.system("mkdir SN_"+str(sn))
	director = "SN_"+str(sn)
	print(sn)
	h = open("sofia.in",'r')
	for line in h:
		if '#' not in line:
			N= np.size(line.split())
			par = line.split()[0]
			flag= line.split()[4]
			imsize=int(line.split()[2])
			zsize=int(line.split()[3])
			flag =flag.split(',')
			
			print flag
			##find the flagging edges
			edges = [flag[0]]

			for i in range(1,np.size(flag)):

				if (int(flag[i]) - int(flag[i-1]) )!= 1:
					
				# 	print flag[i], pr

				# 	edges.append(pr)
				 	edges.append(flag[i-1])
				 	edges.append(flag[i])
				# 	pr =flag[i]
				# 	print pr
			edges.append(flag[-1])
			## create flagged array
			flagged_array=[]
			for i in range(0,np.size(edges),2):
			    ar= [0,imsize,0,imsize,int(edges[i]),int(edges[i+1])+1]
			    flagged_array.append(ar)

			print "flagged aray",flagged_array

			print par
			fin=open('template.param','r') #we gan give here an uniform directory
			fout=open(par.replace("_mock_nozero.fits",".sofia"),'wt')

			for line2 in fin:
				a=line2.split()
				if "import.inFile" in a:
					fout.write(line2.replace(a[2], par))
				elif "flag.regions" in a:
					fout.write(line2.replace(a[2], str(flagged_array)))
				elif "SCfind.threshold" in a:
					fout.write(line2.replace(a[2], str(sn)))
				elif  "writeCat.outputDir" in a:
					line = "=  "+director
					fout.write(line2.replace(a[1], line))
				# elif "import.subcube" in a:
				# #find channels with 0 flux 
				# 	if N > 5:
				# 		zeros = line.split()[5]
				# 		zeros = zeros.split(',')
				# 		# for creating zero array there cannot be any zeroes inside the cube - cut the cubes
				# 		top,bottom=[],[]
				# 		for z in zeros:
				# 			if int(z) < zsize/2.:
				# 				top=np.append(top,z)
				# 			else:
				# 				bottom=np.append(bottom,z)

				# 		top,bottom=np.array(top, dtype=int), np.array(bottom,dtype=int)
				# 		#print np.size(top),np.size(bottom)
				# 		if np.size(bottom) < 1:
				# 			z1,z2=max(top)+1, zsize
				# 		elif np.size(top) < 1:
				# 			z1,z2=0,min(bottom)-1

				# 		else:
				# 			z1,z2=max(top)+1, min(bottom)-1

				# 		zerojy_array = [0,imsize,0,imsize,z1,z2]
				# 		print "subcube",zerojy_array
				# 		fout.write(line2.replace(a[2], str(zerojy_array)))
					# else:
					# 	fout.write(line2)
				else:
					fout.write(line2)
			fin.close()
			fout.close()

			os.system("/Users/ahamanow/SoFiA/SoFiA-master/sofia_pipeline.py "+ par.replace("_mock_nozero.fits",".sofia"))		

