################################################################
# SQUAD-MATS - second-order PT solver for the sc quantum dot   #
# Copyright (C) 2018-2020  Vladislav Pokorny; pokornyv@fzu.cz  #
# homepage: github.com/pokornyv/SQUAD-MATS                     #
# spectrum_ddot_map.py - generate color maps of spec. func.    #
# method described in                                          #
#    Sci. Rep. 5, 8821 (2015).                                 #
#    Phys. Rev. B 93, 024523 (2016).                           #
################################################################

import numpy as np
import scipy as sp
from sys import argv
from time import ctime

xmin = float(argv[1])
xmax = float(argv[2])
dx   = float(argv[3])

fname = str(argv[4])		## e.g. 'gr'

try:
	gf = argv[5]
except IndexError:
	gf = 'normal'

try:
	dot = argv[6]
except IndexError:
	dot = 'L'	## left or right dot

print('# '+ctime())
print('# Reading data for dot '+str(dot))

emin = -2.0	## max and min energy for output to avoid huge files
emax =  2.0

if gf == 'normal':
	fname_up  = 'Gup_normal'+str(dot)+'.dat'
	fname_tot = 'Gtot_normal'+str(dot)+'.dat'
elif  gf == 'anomalous':
	fname_up  = 'Gup_anom'+str(dot)+'.dat'
	fname_tot = 'Gtot_anom'+str(dot)+'.dat'
elif  gf == 'interdot':
	fname_up  = 'Gup_inter.dat'
	fname_tot = 'Gtot_inter.dat'
else:
	print('Wrong GF block. Select normal, anomalous or interdot.')

print('# Reading data for '+str(gf)+' Green function.')

## blocks:
# 0  1  2  3   0:1  1:2    2:3     x
# 4  5  6  7   4:5  5:6      x   7:8
# 8  9 10 11   8:9    x  10:11 11:12
#12 13 14 15     x 13:14 14:15 15:16

with open(fname_up,'w') as fout1:
	with open(fname_tot,'w') as fout2:
	
		x = xmin
		while x <= xmax+dx:
			E = 0
			fn = fname+str(np.around(x,3))
			try:
				Data_A = np.loadtxt(fn)
				print('Reading '+fn)
			except IOError:
				#print('Error reading '+fn)
				E = 1
			if E:
				try:
					Data_A = np.loadtxt(fn+'0')
					print('Reading '+fn+'0')
				except IOError:
					print('File '+fn+'0'+' not found.')
					x += dx
					continue
			N = len(Data_A)
			if gf == 'normal':
				if dot == 'L': 
					Data1_A = Data_A[:int(N/16)]
					Data2_A = np.flipud(Data_A[int(5*N/16):int(6*N/16)])
				else: # dot=='R'
					Data1_A = Data_A[int(10*N/16):int(11*N/16)]
					Data2_A = np.flipud(Data_A[int(15*N/16):])
			elif  gf == 'anomalous':
				if dot == 'L': Data1_A = Data_A[int(N/16):int(2*N/16)]     ## left dot
				else:          Data1_A = Data_A[int(11*N/16):int(12*N/16)] ## right dot
			elif  gf == 'interdot':
				Data1_A = Data_A[int(2*N/16):int(3*N/16)]
			for i in range(int(N/16)):
				if Data1_A[i,0]>emin and Data1_A[i,0]<emax:
					fout1.write('{0: .6f}\t{1:.12f}\t{2:.12f}\n'\
					.format(x,Data1_A[i,0],Data1_A[i,2]))
					if gf == 'normal': fout2.write('{0: .6f}\t{1:.12f}\t{2:.12f}\n'\
					.format(x,Data1_A[i,0],Data1_A[i,2]+Data2_A[i,2]))
			fout1.write('\n')
			fout2.write('\n')
			x += dx

## spectrum_map_ddot.py END

