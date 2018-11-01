## second-order PT solver for the superconducting quantum dot
## Matsubara frequency formulation
## uses TRIQS 1.4
## Vladislav Pokorny; 2018; pokornyv@fzu.cz

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
	Data4_A = sp.flipud(Data_A[3*N/4:N])
	for i in range(N/4):
		if Data1_A[i,0]>emin and Data1_A[i,0]<emax:
			fout1.write('{0: .6f}\t{1:.12f}\t{2:.12f}\n'\
			.format(float(x),float(Data1_A[i,0]),float(Data1_A[i,2])))
			fout2.write('{0: .6f}\t{1:.12f}\t{2:.12f}\n'\
			.format(float(x),float(Data1_A[i,0]),float(Data1_A[i,2]+Data4_A[i,2])))
	fout1.write('\n')
	fout2.write('\n')
	x += dx

fout1.close()
fout2.close()

