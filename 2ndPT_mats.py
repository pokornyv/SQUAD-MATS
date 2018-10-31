## second-order PT solver for the superconducting quantum dot
## Matsubara frequency formulation
## uses TRIQS 1.4
## Vladislav Pokorny; 2018; pokornyv@fzu.cz

import scipy as sp
from sys import argv,exit
from os import getcwd
from time import time,ctime
from matslib import *

from pytriqs.gf.local import *
from pytriqs.archive import *

stars  = '*'*60
hashes = '#'*60

## HDF5 output
hdf5file = 'squad_mats.h5'

## convergence criteria
eps_hf  = 1e-4
eps_2nd = 1e-4

U      = float(argv[1])
beta   = float(argv[2])
Delta  = 1.0
GammaL = 0.5
GammaR = 0.5
GammaN = 0.0
P      = 0.5
eps    = float(argv[3])
try:
	B = float(argv[4])
except IndexError:
	B = 0.0
iw_cut = 100.0	 # energy cutoff in Matsubaras
BW     = 1000.0 # half-bandwidth

## mixing parameters
alpha_hf  = 0.2
alpha_2nd = 0.2

## Pade continuation parameters
emax   = 10.0
NRealP = 20000
NMats_cont = 50
izero = 1e-3

bands_T  = ['up','dn']
params_A = [beta,U,Delta,GammaL,GammaR,GammaN,eps,P,B]

NMats = int(0.5*(beta*iw_cut/sp.pi)+1.0)
logfname = '2ndPT.log'

###########################################################
PrintAndWrite(hashes+'\nWorking directory: '+str(getcwd()),logfname)
PrintAndWrite('output generated by '+str(argv[0])+' on '+ctime(),logfname)
PrintAndWrite(' U ={0: .3f}, Delta ={1: .3f},\n GammaR ={2: .3f}, \
GammaL ={3: .3f}, GammaN ={4: .3f},\n eps ={5: .3f}, B ={6: .3f}, Phi/pi ={7: .3f}'\
.format(U,Delta,GammaR,GammaL,GammaN,eps,B,P),logfname)
PrintAndWrite('Inverse temperature beta ={0: .3f}, half-bandwidth W = {1: .3f}'\
.format(beta,BW),logfname)
PrintAndWrite('Using {0: 3d} Matsubara frequencies, cutoff: {1: .3f}'\
.format(int(NMats),float(iw_cut)),logfname)

###########################################################
## Hartree-Fock calculation ###############################

## initial condition
## it is important to change m while searching for polarized solutions
[n,m,mu] = [0.5,0.001,0.2]

PrintAndWrite('\nCalculating Hartree-Fock solution:',logfname)
[nold,mold,muold] = [1e5,1e5,1e5]
N0_A = sp.array([[(m+n)/2.0,mu],[mu,(m-n)/2.0+1.0]])
PrintAndWrite('\n    iter\tn\t\tm\t\tmu',logfname)
t = time()
niter = 0
while any([sp.fabs(n-nold) > eps_hf,sp.fabs(m-mold) > eps_hf,sp.fabs(mu-muold) > eps_hf]):
	Nold_A = N0_A
	[nold,mold,muold] = Occupation(Nold_A)
	G0_iw = GFzero(params_A,bands_T,N0_A,BW,NMats,int(0.8*NMats),NMats,8)
	N0_A = TotalDensity(G0_iw)
	[n,m,mu] = Occupation(N0_A)
	N0_A = alpha_hf*N0_A+(1.0-alpha_hf)*Nold_A
	niter += 1
	PrintAndWrite('  {0: 3d}\t{1: .8f}\t{2: .8f}\t{3: .8f}'.format(niter,n,m,mu),logfname)
PrintAndWrite('  ...done in {0: 3d} seconds.'.format(int(time()-t)),logfname)

WriteG_iw(G0_iw,'gwhf',logfname)

NZero_A = TotalDensity(G0_iw)
PrintAndWrite('\n  Occupation matrix from G_HF(iw):',logfname)
WriteMatrix(NZero_A,bands_T,'M',logfname)

G0_real = PadeContinuation(G0_iw,emax,NRealP,NMats_cont,izero)
WriteG_real(G0_real,'grhf',logfname)
WriteG_real(G0_real,'grhf'+str(B),logfname)

PrintAndWrite('{0: .4f}\t{1: .4f}\t{2: .4f}\t{3: .8f}\t{4: .8f}\t{5: .8f}'.format(U,eps,B,n,m,mu),logfname)

###########################################################
## second order correction ################################

PrintAndWrite('\nCalculating the two-particle bubble...',logfname)
Chi = TwoParticleBubbleFFT(G0_iw,G0_iw)
WriteG_iw(Chi,'chiw',logfname)
'''
if 0:
	PrintAndWrite('  Hartree-Fock data saved to '+hdf5file+' archive.',logfname)
	R = HDFArchive(hdf5file,'w')
	R['g0_iw'] = G0_iw
	R['chi']   = Chi
	R['n']     = n
	R['mu']    = mu
	del R

if 0:
	PrintAndWrite('  Hartree-Fock data read from '+hdf5file+' archive.',logfname)
	R = HDFArchive(hdf5file,'r')
	Chi   = R['chi']
	G0_iw = R['g0_iw']
	n     = R['n']
	mu    = R['mu']
	del R
'''
Chi_real = PadeContinuation(Chi,emax,NRealP,50,izero)
WriteG_real(Chi_real,'chir',logfname)

PrintAndWrite('\nCalculating the kernel of the Schwinger-Dyson equation...',logfname)
Psi = GfImFreq(indices = [0], beta = beta, n_points = NMats, statistic = 'Boson')
Psi << U**2*(Chi['up','up']+Chi['up','dn'])
WriteG_iw(Psi,'psiw',logfname)
Psi_real = PadeContinuation(Psi,emax,NRealP,50,izero)
WriteG_real(Psi_real,'psir',logfname)

PrintAndWrite('\nCalculating the 2ndPT dynamic self-energy...',logfname)
Sigma = SelfEnergy(G0_iw,Psi)
WriteG_iw(Sigma,'sw',logfname)
Sigma_real = PadeContinuation(Sigma,emax,NRealP,50,izero)
WriteG_real(Sigma_real,'sr',logfname)

Gint_iw = G0_iw.copy()
Gint_iw << inverse(inverse(G0_iw)-Sigma)	## Dyson equation

PrintAndWrite('\nCorrecting the Hartree-Fock self-energy...',logfname)
N_A = N0_A.copy()
[nold,mold,muold] = [1e5,1e5,1e5]
[n,m,mu] = Occupation(N_A)
#[n,mu] = [N0_A[0][0],N0_A[0][1]]
PrintAndWrite('\n    iter\tn\tm\t\tmu',logfname)
niter = 0
while any([sp.fabs(n-nold) > eps_2nd,sp.fabs(m-mold) > eps_2nd,sp.fabs(mu-muold) > eps_2nd]):
	Nold_A = N_A
	[nold,mold,muold] = Occupation(N_A)
	G0_iw = GFzero(params_A,bands_T,N_A,BW,NMats,int(0.9*NMats),NMats,6)
	Gint_iw << inverse(inverse(G0_iw)-Sigma)	## Dyson equation
	N_A = TotalDensity(Gint_iw)
	N_A = alpha_2nd*N_A+(1.0-alpha_2nd)*Nold_A
	[n,m,mu] = Occupation(N_A)
	niter += 1
	PrintAndWrite('  {0: 3d}\t{1: .8f}\t{2: .8f}\t{3: .8f}'.format(niter,n,m,mu),logfname)

PrintAndWrite('\nCalculating the 2ndPT interacting Green function...',logfname)

WriteG_iw(Gint_iw,'gw',logfname)

N_A = TotalDensity(Gint_iw)
PrintAndWrite('\n  Occupation matrix from G(iw):',logfname)
WriteMatrix(N_A,bands_T,'M',logfname)
[n,m,mu] = Occupation(N_A)

Gint_real = PadeContinuation(Gint_iw,emax,NRealP,50,izero)
WriteG_real(Gint_real,'gr'+str(B),logfname)
#WriteG_real(Gint_real,'gr',logfname)

PrintAndWrite('{0: .4f}\t{1: .4f}\t{2: .4f}\t{3: .8f}\t{4: .8f}\t{5: .8f}'.format(U,eps,B,n,m,mu),logfname)

PrintAndWrite('# '+str(argv[0])+' END, '+ctime(),logfname)

## 2ndPT_mats.py end ##


