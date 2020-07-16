################################################################
# SQUAD-MATS - second-order PT solver for the sc quantum dot   #
# Copyright (C) 2018-2020  Vladislav Pokorny; pokornyv@fzu.cz  #
# homepage: github.com/pokornyv/SQUAD-MATS                     #
# matslib.py - library of functions                            #
# method described in                                          #
#    Sci. Rep. 5, 8821 (2015).                                 #
#    Phys. Rev. B 93, 024523 (2016).                           #
################################################################

from config import *
from iolib import *

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


def Occupation(N_A):
	''' calculates n,m,mu from the occupation matrix N '''
	n  =  0.5*(N_A[0][0]-N_A[1][1]+1.0)
	m  =  0.5*(N_A[0][0]+N_A[1][1]-1.0)
	mu =  0.5*(N_A[0][1]+N_A[1][0])
	return [n,m,mu]


###########################################################
## non-interacting Green function #########################

def IntFiniteBW(Delta,iw):
	''' integral enering hybridizations for finite bandwidth 2W '''
	sq = lambda x: np.sqrt(Delta**2+np.imag(x)**2)
	return 2.0/(np.pi*sq(iw))*np.arctan2(BW,sq(iw))


def HybDiag(iw):
	''' diagonal part of the finite-bandwidth sc lead hybridization, real 
	does not contain the iw_n factor from the matrix '''
	return (GammaL+GammaR)*IntFiniteBW(Delta,iw)


def HybOffDiag(iw):
	''' off-diagonal part of the finite-bandwidth sc lead hybridization 
	PhiS angle keeps the hybridization real '''
	PhiS = np.arctan((GammaL-GammaR)/(GammaL+GammaR+1e-12)*np.tan(Phi/2.0))
	return Delta*np.exp(1.0j*PhiS)*\
	(GammaL*np.exp(-1.0j*Phi/2.0)+GammaR*np.exp(1.0j*Phi/2.0))*IntFiniteBW(Delta,iw)


def GFzero(N_A,FitMin,FitMax,NTerms):
	''' constructs the Hartree-Fock Green function as the input '''
	[n,m,mu] = Occupation(N_A)
	V2 = BW*GammaN/sp.pi  # hybridization with normal lead
	## define lambdas (hybridizations are real)
	[n,m,mu] = Occupation(N_A)
	GFinv11 = lambda x: x*(1.0 + HybDiag(x)) - (eps + U*(n-0.5)) + (h + U*m)
	GFinv12 = lambda x: HybOffDiag(x) - U*mu
	GFinv21 = lambda x: sp.conj(HybOffDiag(x)) - U*mu
	GFinv22 = lambda x: x*(1.0 + HybDiag(x)) + (eps + U*(n-0.5)) + (h + U*m)
	## define GF objects
	GFinv = GfImFreq(indices = bands_T,beta = beta,n_points = NMats)
	GF    = GfImFreq(indices = bands_T,beta = beta,n_points = NMats)
	## fill GF objects with lambdas
	## BW is half-bandwidth in Flat() descriptor from TRIQS
	GFinv[bands_T[0],bands_T[0]] << Function(GFinv11) - V2*Flat(BW)
	GFinv[bands_T[0],bands_T[1]] << Function(GFinv12)
	GFinv[bands_T[1],bands_T[0]] << Function(GFinv21)
	GFinv[bands_T[1],bands_T[1]] << Function(GFinv22) - V2*Flat(BW)
	## calculate inverse
	GF << inverse(GFinv)
	## fit the tail
	moments_A = np.zeros((0,2,2),dtype=complex)
	## tail object is a tuple!
	tail = fit_tail_on_window(GF, n_min=FitMin,n_max=FitMax,n_tail_max=10*FitMax,known_moments=moments_A,expansion_order=8)
	replace_by_tail(GF,sp.array(tail[0]),n_min=FitMin)
	if GF.is_gf_hermitian() and mpi.is_master_node():
		PrintAndWrite(' - Non-interacting Green function is Hermitean',outfile)
	return GF

###########################################################
## Pade analytic continuation #############################

def PadeContinuation(GFw):
	''' Pade analytic continuation for the Matsubara Green function '''
	GFr = GfReFreq(indices = GFw.indices, window = (-pade_emax,pade_emax), n_points = pade_NReal)
	GFr.set_from_pade(GFw, n_points = pade_NMats, freq_offset = pade_izero)
	return GFr

###########################################################
## convolutions ###########################################

def TwoParticleBubbleFFT(GF1,GF2t):
	''' calculate the two-particle bubble using FFT '''
	NMats = len(GF1.data)/2
	GF1tau = GfImTime(indices = GF1.indices, beta = GF1.mesh.beta, n_points = 6*NMats+1)
	GF2tau = GfImTime(indices = GF1.indices, beta = GF1.mesh.beta, n_points = 6*NMats+1)
	Chitau = GfImTime(indices = GF1.indices, beta = GF1.mesh.beta, n_points = 6*NMats+1, statistic = 'Boson')
	Chi    = GfImFreq(indices = GF1.indices, beta = GF1.mesh.beta, n_points = NMats, statistic = 'Boson')
	GF2 = GF2t.copy()
	GF2['up','up'].data[:] =	GF2t['dn','dn'].data[:]
	GF2['up','dn'].data[:] =	GF2t['dn','up'].data[:]
	GF2['dn','up'].data[:] =	GF2t['up','dn'].data[:]
	GF2['dn','dn'].data[:] =	GF2t['up','up'].data[:]
	GF1tau << Fourier(GF1)
	GF2tau << Fourier(GF2)
	flip_A = sp.ones([NBand,NBand])-2*sp.eye(NBand)
	Chitau.data[:] = flip_A*GF1tau.data[:]*GF2tau.data[:]
	WriteG_tau(Chitau,'chitau')
	Chi << Fourier(Chitau)
	## fitting the tail	
	moments_A = np.zeros((0,2,2),dtype=complex)
	## Chi contains artefacts from FFT at high frequencies, do not fit up to NMats
	fitmin = int(0.8*NMats)
	fitmax = int(0.9*NMats)
	tail = fit_hermitian_tail_on_window(Chi,n_min=fitmin,n_max=fitmax,n_tail_max=10*fitmax,known_moments=moments_A,expansion_order=6)
	#print(tail[0])
	#replace_by_tail(Chi,sp.array(tail[0]),n_min=6)
	return Chi


def SelfEnergy(GF,Psi):	
	''' calculate the self-energy from SD equation using FFT '''
	NMats = len(GF.data)/2
	Sigma = GF.copy()
	Sigmatau = GfImTime(indices = GF.indices, beta = GF.mesh.beta, n_points = 6*NMats+1)
	GFtau  = GfImTime(indices = GF.indices, beta = GF.mesh.beta, n_points = 6*NMats+1)
	Ktau   = GfImTime(indices = GF.indices, beta = GF.mesh.beta, n_points = 6*NMats+1, statistic = 'Boson')
	Kernel = GfImFreq(indices = GF.indices, beta = GF.mesh.beta, n_points = NMats, statistic = 'Boson')
	for i,j in product(GF.indices[0], repeat = 2): 
		Kernel[i,j].data[:] = Psi.data[:,0,0]
	GFtau << Fourier(GF)
	Ktau  << Fourier(Kernel)
	Sigmatau.data[:] = -GFtau.data[:]*Ktau.data[:]
	WriteG_tau(Sigmatau,'sigmatau')
	Sigma << Fourier(Sigmatau)
	## fitting the tail
	moments_A = np.zeros((0,2,2),dtype=complex)
	## Sigma contains artefacts from FFT at high frequencies, do not fit up to NMats
	fitmin = int(0.8*NMats)
	fitmax = int(0.9*NMats)
	tail = fit_hermitian_tail_on_window(Sigma,n_min=fitmin,n_max=fitmax,n_tail_max=10*fitmax,known_moments=moments_A,expansion_order=6)
	#print(tail[0])
	#replace_by_tail(Sigma,sp.array(tail[0]),n_min=6)
	return Sigma

## matslib.py end ##

