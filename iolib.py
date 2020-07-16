################################################################
# SQUAD-MATS - second-order PT solver for the sc quantum dot   #
# Copyright (C) 2018-2020  Vladislav Pokorny; pokornyv@fzu.cz  #
# homepage: github.com/pokornyv/SQUAD-MATS                     #
# iolib.py - input/output library                              #
# method described in                                          #
#    Sci. Rep. 5, 8821 (2015).                                 #
#    Phys. Rev. B 93, 024523 (2016).                           #
################################################################

from config import *

def PrintAndWrite(line,fname):
	'''	print the same line to stdout and to file fname '''
	print(line)
	with open(fname,'a') as f: f.write(line+'\n')


def WriteG_iw(GF,fname):
	''' writes a Matsubara function to a file '''
	MatsFreq_F = sp.zeros(len(GF.data))
	k = 0
	for iw in GF.mesh:
		MatsFreq_F[k] = sp.imag(iw)
		k += 1
	with open(fname,'w') as fout:
		fout.write('# file written on '+ctime()+'\n')
		for i,j in product(GF.indices[0], repeat = 2):
			for nw in range(len(GF.data)):
				fout.write('{0:.8f}\t{1:.12f}\t{2:.12f}\n'\
				.format(float(MatsFreq_F[nw]),float(sp.real(GF[i,j].data[nw]))\
				,float(sp.imag(GF[i,j].data[nw]))))
			fout.write('\n\n')
	PrintAndWrite('  File '+fname+' written.',logfname)


def WriteG_tau(GF,fname):
	''' writes a Matsubara function to a file '''
	Tau_F = sp.zeros(len(GF.data))
	k = 0
	for tau in GF.mesh:
		Tau_F[k] = sp.real(tau)
		k += 1
	with open(fname,'w') as fout:
		fout.write('# file written on '+ctime()+'\n')
		for i,j in product(GF.indices[0], repeat = 2):
			for tau in range(len(GF.data)):
				fout.write('{0:.12f}\t{1:.12f}\t{2:.12f}\n'\
				.format(float(Tau_F[tau]),float(sp.real(GF[i,j].data[tau]))\
				,float(sp.imag(GF[i,j].data[tau]))))
			fout.write('\n\n')
	PrintAndWrite('  File '+fname+' written.',logfname)


def WriteG_real(GF,fname):
	''' writes a function of a real frequency to a file '''
	Freq_F = np.zeros(len(GF.data))
	k = 0
	for w in GF.mesh:
		Freq_F[k] = np.real(w)
		k += 1
	with open(fname,'w') as fout:
		fout.write('# file written on '+ctime()+'\n')
		for i,j in product(GF.indices[0], repeat = 2):
			for nw in range(len(Freq_F)):
				fout.write('{0:.12f}\t{1:.12f}\t{2:.12f}\n'\
				.format(float(Freq_F[nw]),float(np.real(GF[i,j].data[nw]))\
				,float(np.imag(GF[i,j].data[nw]))))
			fout.write('\n\n')
	PrintAndWrite('  File '+fname+' written.',logfname)


def WriteMatrix(X_A):
	''' writes a matrix to output in a formatted way '''
	out_text = ''
	for i in range(X_A.shape[0]):
		for j in range(X_A.shape[1]):
			X = X_A[i][j]
			if sp.imag(X) < 1e-9: out_text += '{0: .6f}\t'.format(float(np.real(X)))
			else: out_text += '{0: .6f}{1:+0.6f}i\t'.format(float(np.real(X)),float(np.imag(X)))
		out_text += '\n'
	PrintAndWrite(out_text,logfname)

## iolib.py end ##

