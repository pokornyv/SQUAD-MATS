## second-order PT solver for the superconducting quantum dot
## Matsubara frequency formulation
## uses TRIQS 1.4
## Vladislav Pokorny; 2018; pokornyv@fzu.cz

import scipy as sp
from time import ctime
from itertools import product
from pytriqs.gf.local import *
from pytriqs.gf.local.descriptors import Function

###########################################################
## general functions ######################################

def OmegaN(n,beta):
	''' calculates the n-th fermionic Matsubara frequency '''
	return (2.0*n+1.0)*sp.pi/beta


def NuM(m,beta):
	''' calculates the m-th bosonic Matsubara frequency '''
	return 2.0*m*sp.pi/beta


def Brillouin(J,x):
	''' Brillouin function B_J(x)	'''
	a = (2.0*J+1.0)/(2.0*J)
	b = 1.0/(2.0*J)
	return 0.0 if x==0.0 else a/tanh(a*x)-b/tanh(x)


def Langevin(x):
	''' Langevin function: classical limit of the J=1/2 Brillouin function '''
	return 0.0 if x==0.0 else 1.0/sp.tanh(x)-1.0/x


def PrintAndWrite(line,fname):
	'''	print the same line to stdout and to file fname '''
	print(line)
	f = open(fname,'a')
	f.write(line+'\n')
	f.close()


def TotalDensity(G):
	''' calculates the density from a Green function '''
	bands_T = G.indices
	NBand = len(bands_T)
	N_F = sp.zeros([NBand,NBand])
	for i,j in product(range(NBand), repeat = 2):
		N_F[i][j] = sp.real(G[bands_T[i],bands_T[j]].total_density())
	return N_F


def Occupation(N_A):
	''' calculates n,m,mu from occupation matrix N '''
	n  = 0.5*(N_A[0][0]-N_A[1][1]+1.0)
	m  = 0.5*(N_A[0][0]+N_A[1][1]-1.0)
	mu = 0.5*(N_A[0][1]+N_A[1][0])
	return [n,m,mu]

###########################################################
## non-interacting Green function #########################

def IntFiniteBW(Delta,W,iw):
	''' integral enering hybridizations for finite bandwidth 2W '''
	sq = lambda x: sp.sqrt(Delta**2+sp.imag(x)**2)
	return 2.0/(sp.pi*sq(iw))*sp.arctan2(W,sq(iw))


def HybDiag(GammaL,GammaR,Delta,W,iw):
	''' diagonal part of the finite-bandwidth sc lead hybridization, real 
	    caution, does not contain the iw_n factor from the matrix!!! '''
	return (GammaL+GammaR)*IntFiniteBW(Delta,W,iw)


def HybOffDiag(GammaL,GammaR,Delta,P,W,iw):
	''' off-diagonal part of the finite-bandwidth sc lead hybridization 
	PhiS angle keeps the hybridization real '''
	Phi = P*sp.pi
	PhiS = sp.arctan((GammaL-GammaR)/(GammaL+GammaR+1e-12)*sp.tan(Phi/2.0))
	return Delta*sp.exp(1.0j*PhiS)*\
	(GammaL*sp.exp(-1.0j*Phi/2.0)+GammaR*sp.exp(1.0j*Phi/2.0))*IntFiniteBW(Delta,W,iw)


def GFzero(params_A,bands_T,N_A,BW,NMats,FitMin,FitMax,NTerms):
	''' constructs the non-interacting Green function as the input '''
	[beta,U,Delta,GL,GR,GN,eps,P,B] = params_A
	[n,m,mu] = Occupation(N_A)
	V2 = BW*GN/sp.pi  # hybridization with normal lead
	## define lambdas (hybridizations are real)
	GFinv11 = lambda x: x*(1.0 + HybDiag(GL,GR,Delta,BW,x)) - (eps + U*(n-0.5)) + (B + U*m) 
	GFinv12 = lambda x: HybOffDiag(GL,GR,Delta,P,BW,x) - U*mu
	GFinv21 = lambda x: sp.conj(HybOffDiag(GL,GR,Delta,P,BW,x)) - U*mu
	GFinv22 = lambda x: x*(1.0 + HybDiag(GL,GR,Delta,BW,x)) + (eps + U*(n-0.5)) + (B + U*m)
	## define GF objects
	GFinv = GfImFreq(indices = bands_T,beta = beta,n_points = NMats)
	GF    = GfImFreq(indices = bands_T,beta = beta,n_points = NMats)
	## fill GF objects with lambdas
	## BW is half-bandwidth in Flat() descriptor from TRIQS
	GFinv[bands_T[0],bands_T[0]] << Function(GFinv11) - V2*Flat(BW)
	GFinv[bands_T[0],bands_T[1]] << Function(GFinv12)
	GFinv[bands_T[1],bands_T[0]] << Function(GFinv21)
	GFinv[bands_T[1],bands_T[1]] << Function(GFinv22) - V2*Flat(BW)
	## fit the tail
	GFinv.tail.zero()
	fixed_tail = TailGf(2,2,1,-1)
	fixed_tail[-1] = sp.eye(2)
	GFinv.fit_tail(fixed_tail,NTerms,FitMin,FitMax)
	## calculate inverse
	GF << inverse(GFinv)
	## refit the tail, just in case
	GF.tail.zero()
	fixed_tail = TailGf(2,2,3,-1)
	fixed_tail[-1] = sp.zeros([2,2])
	fixed_tail[ 0] = sp.zeros([2,2])
	fixed_tail[ 1] = sp.eye(2)
	GF.fit_tail(fixed_tail,NTerms,FitMin,FitMax)
	return GF


###########################################################
## Pade analytic continuation #############################

def PadeContinuation(GFw,emax,NRealP,NMats,izero):
	''' Pade analytic continuation for the Matsubara Green function '''
	GFr = GfReFreq(indices = GFw.indices, window = (-emax,emax), n_points = NRealP)
	GFr.set_from_pade(GFw, n_points = NMats, freq_offset = izero)
	return GFr

###########################################################
"""
def RollGF(GF,m,NM):
	''' roll the data in Green function by a given bosonic frequency 
	return the shifted data array extended to [-NMats,NMats] interval '''
	NBand = len(GF.indices)
	Data_A = sp.zeros([2*NM,NBand,NBand],dtype=complex)
	offset = NM-len(GF.data)/2
	for i in range(len(Data_A)):
		if i+m < offset or i+m >= NM+len(GF.data)/2:
			Data_A[i] = GF.tail(1.0j*OmegaN(i+m-NM,GF.beta))
		else:
			Data_A[i] = GF.data[i-offset+m]
	return Data_A


def TwoParticleBubble(GF1,GF2):
	''' calculate the two-particle bubble as a Matsubara sum '''
	NMats = len(GF1.data)/2
	Chi = GfImFreq(indices = GF1.indices, beta = GF1.beta, n_points = NMats, statistic = 'Boson')
	Int = GfImFreq(indices = GF1.indices, beta = GF1.beta, n_points = 2*NMats)
	m = 0
	for inu in Chi.mesh:  ## loop over bosonic frequency
		Int.zero()
		G1data_A = RollGF(GF1,0,2*NMats)
		G2data_A = RollGF(GF2,m-NMats+1,2*NMats)
		Int.data[:] = G1data_A[:]*G2data_A[:]
		Chi.data[m] = Int.density()-sp.array([[0.0,0.0],[0.0,0.0]])
		m += 1
	return Chi
"""

def TwoParticleBubbleFFT(GF1,GF2t):
	''' calculate the two-particle bubble using FFT '''
	NMats = len(GF1.data)/2
	GF1tau = GfImTime(indices = GF1.indices, beta = GF1.beta, n_points = 4*NMats+1)
	GF2tau = GfImTime(indices = GF1.indices, beta = GF1.beta, n_points = 4*NMats+1)
	Chitau = GfImTime(indices = GF1.indices, beta = GF1.beta, n_points = 4*NMats+1, statistic = 'Boson')
	Chi    = GfImFreq(indices = GF1.indices, beta = GF1.beta, n_points = NMats, statistic = 'Boson')
	GF2 = GF2t.copy()
	GF2['up','up'].data[:] =	GF2t['dn','dn'].data[:]
	GF2['up','dn'].data[:] =	GF2t['dn','up'].data[:]
	GF2['dn','up'].data[:] =	GF2t['up','dn'].data[:]
	GF2['dn','dn'].data[:] =	GF2t['up','up'].data[:]
	GF1tau << InverseFourier(GF1)
	GF2tau << InverseFourier(GF2)
	flip_A = sp.array([[-1.0,1.0],[1.0,-1.0]])
	Chitau.data[:] = flip_A*GF1tau.data[:]*GF2tau.data[:]
	#WriteG_tau(Chitau,'chit','log.log')
	Chi << Fourier(Chitau)
	## fitting the tail
	Chi.tail.zero()
	fixed_tail = TailGf(2,2,2,-1)
	fixed_tail[-1] = sp.zeros([2,2])
	fixed_tail[ 0] = sp.zeros([2,2])
	Chi.fit_tail(fixed_tail,6,int(0.8*NMats),int(0.9*NMats))
	return Chi


def SelfEnergy(GF,Psi):	
	''' calculate the self-energy from SD equation using FFT '''
	NMats = len(GF.data)/2
	Sigma = GF.copy()
	Sigmatau = GfImTime(indices = GF.indices, beta = GF.beta, n_points = 4*NMats+1)
	GFtau = GfImTime(indices = GF.indices, beta = GF.beta, n_points = 4*NMats+1)
	Ktau  = GfImTime(indices = GF.indices, beta = GF.beta, n_points = 4*NMats+1, statistic = 'Boson')
	Kernel = GfImFreq(indices = GF.indices, beta = GF.beta, n_points = NMats, statistic = 'Boson')
	for i,j in product(GF.indices, repeat = 2):
		Kernel[i,j] << Psi
	GFtau << InverseFourier(GF)
	Ktau  << InverseFourier(Kernel)
	flip_A = sp.array([[-1.0,-1.0],[-1.0,-1.0]])
	Sigmatau.data[:] = flip_A*GFtau.data[:]*Ktau.data[:]
	#WriteG_tau(Sigmatau,'st','log.log')
	Sigma << Fourier(Sigmatau)
	## fitting the tail
	Sigma.tail.zero()
	fixed_tail = TailGf(2,2,2,-1)
	fixed_tail[-1] = sp.zeros([2,2])
	fixed_tail[ 0] = sp.zeros([2,2])
	Sigma.fit_tail(fixed_tail,6,int(0.8*NMats),int(0.9*NMats))
	return Sigma


###########################################################
## functions for writing data files #######################

def WriteG_iw(GF,fname,logfname):
	''' writes a Matsubara function to a file '''
	fout = open(fname,'w')
	fout.write('# file written on '+ctime()+'\n')
	NMats = len(GF.data)
	bands_T = GF.indices
	NBand = len(bands_T)
	MatsFreq_F = sp.zeros(NMats)
	k = 0
	for iw in GF.mesh:
		MatsFreq_F[k] = sp.imag(iw)
		k += 1
	for i,j in product(range(NBand), repeat = 2):
		for nw in range(NMats):
			fout.write('{0:.12f}\t{1:.12f}\t{2:.12f}\n'\
			.format(float(MatsFreq_F[nw]),float(sp.real(GF[bands_T[i],bands_T[j]].data[nw][0][0]))\
               ,float(sp.imag(GF[bands_T[i],bands_T[j]].data[nw][0][0]))))
		fout.write('\n\n')
	fout.close()
	PrintAndWrite('  File '+fname+' written.',logfname)


def WriteG_tau(GF,fname,logfname):
	''' writes a Matsubara function to a file '''
	fout = open(fname,'w')
	fout.write('# file written on '+ctime()+'\n')
	bands_T = GF.indices
	NBand = len(bands_T)
	NTau = len(GF.data)
	Tau_F = sp.zeros(NTau)
	k = 0
	for tau in GF.mesh:
		Tau_F[k] = sp.real(tau)
		k += 1
	for i,j in product(range(NBand), repeat = 2):
		for tau in range(NTau):
			fout.write('{0:.12f}\t{1:.12f}\t{2:.12f}\n'\
			.format(float(Tau_F[tau]),float(sp.real(GF[bands_T[i],bands_T[j]].data[tau][0][0]))\
               ,float(sp.imag(GF[bands_T[i],bands_T[j]].data[tau][0][0]))))
		fout.write('\n\n')
	fout.close()
	PrintAndWrite('  File '+fname+' written.',logfname)


def WriteG_real(GF,fname,logfname):
	''' writes a function of a real frequency to a file '''
	fout = open(fname,'w')
	fout.write('# file written on '+ctime()+'\n')
	bands_T = GF.indices
	NBand = len(bands_T)
	Freq_F = sp.zeros(len(GF.data))
	k = 0
	for w in GF.mesh:
		Freq_F[k] = sp.real(w)
		k += 1
	for i,j in product(range(NBand), repeat = 2):
		for nw in range(len(Freq_F)):
			fout.write('{0:.12f}\t{1:.12f}\t{2:.12f}\n'\
			.format(float(Freq_F[nw]),float(sp.real(GF[bands_T[i],bands_T[j]].data[nw][0][0]))\
               ,float(sp.imag(GF[bands_T[i],bands_T[j]].data[nw][0][0]))))
		fout.write('\n\n')
	fout.close()
	PrintAndWrite('  File '+fname+' written.',logfname)


def WriteMatrix(X_A,bands_T,Xtype,logfname):
	''' writes a dict or matrix with two indices to output in a matrix form 
	input: Xtype = "D" = dict, "M" = matrix '''
	out_text = ''
	NBand = len(bands_T)
	for i,j in product(range(NBand), repeat = 2):
		if sp.imag(X_A[i][j]) == 0: out_text += '{0: .6f}\t'.format(float(sp.real(X_A[i][j])))
		else: out_text += '{0: .6f}{1:+0.6f}i\t'.format(float(sp.real(X_A[i][j])),float(sp.imag(X_A[i][j])))
		if j==len(bands_T)-1: out_text += '\n'
	PrintAndWrite(out_text,logfname)


## matslib.py end ##

