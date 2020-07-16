################################################################
# SQUAD-MATS - second-order PT solver for the sc quantum dot   #
# Copyright (C) 2018-2020  Vladislav Pokorny; pokornyv@fzu.cz  #
# homepage: github.com/pokornyv/SQUAD-MATS                     #
# spectrum_map.py - generate color maps of spectral functions  #
# method described in                                          #
#    Sci. Rep. 5, 8821 (2015).                                 #
#    Phys. Rev. B 93, 024523 (2016).                           #
################################################################

import numpy as np
import scipy as sp
from sys import argv

xmin = float(argv[1])
xmax = float(argv[2])
dx   = float(argv[3])

fname = str(argv[4])		## e.g. 'gr'

emin = -2.0	## max and min energy for output to avoid huge files
emax =  2.0

fout1 = open('spectrum_Gup.dat','w')
fout2 = open('spectrum_full.dat','w')

with open(fname_up,'w') as fout1:
	with open(fname_tot,'w') as fout2:

		x = xmin
		while x <= xmax+dx:
			try:
				Data_A = sp.loadtxt(fname+str(x))
			except IOError:
				x += dx
				continue
			print('Reading '+fname+str(x))
			N = len(Data_A)
			Data1_A = Data_A[:N/4]
			Data4_A = np.flipud(Data_A[3*N/4:N])
			for i in range(N/4):
				if Data1_A[i,0]>emin and Data1_A[i,0]<emax:
					fout1.write('{0: .6f}\t{1:.12f}\t{2:.12f}\n'\
					.format(float(x),float(Data1_A[i,0]),float(Data1_A[i,2])))
					fout2.write('{0: .6f}\t{1:.12f}\t{2:.12f}\n'\
					.format(float(x),float(Data1_A[i,0]),float(Data1_A[i,2]+Data4_A[i,2])))
			fout1.write('\n')
			fout2.write('\n')
			x += dx

## spectrum_ddot.py END

