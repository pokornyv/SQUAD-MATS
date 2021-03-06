################################################################
# SQUAD-MATS - second-order PT solver for the sc quantum dot   #
# Copyright (C) 2018-2020  Vladislav Pokorny; pokornyv@fzu.cz  #
# homepage: github.com/pokornyv/SQUAD-MATS                     #
# 2PT_mats.py - solver for a single quantum dot                #
# method described in                                          #
#    Sci. Rep. 5, 8821 (2015).                                 #
#    Phys. Rev. B 93, 024523 (2016).                           #
################################################################

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

from matslib import *
from iolib import *

###########################################################
PrintAndWrite(hashes+'\nWorking directory: '+str(getcwd()),logfname)
PrintAndWrite('Output generated by '+str(argv[0])+' on '+ctime(),logfname)
PrintAndWrite(' U ={0: .3f}, Delta ={1: .3f},\n GammaR ={2: .3f}, \
GammaL ={3: .3f}, GammaN ={4: .3f},\n eps ={5: .3f}, h ={6: .3f}, Phi/pi ={7: .3f}'\
.format(U,Delta,GammaR,GammaL,GammaN,eps,h,P),logfname)
PrintAndWrite('Inverse temperature beta ={0: .3f}, kBT = {1: .4f}, half-bandwidth W = {2: .3f}'\
.format(beta,1.0/beta,BW),logfname)
PrintAndWrite('Using {0: 3d} Matsubara frequencies, cutoff: {1: .3f}'\
.format(int(NMats),float(iw_cut)),logfname)
PrintAndWrite('Initial magnetization: m ={0: .3f}'.format(float(magzero)),logfname)
PrintAndWrite('Pade continuation done from first {0: 3d} Matsubara frequencies, iw_max = {1: .5f}'\
.format(int(pade_NMats),float(OmegaN(pade_NMats,beta))),logfname)
if run_FSC: PrintAndWrite('Running the fully self-consistent (FSC) calculation.',logfname)

###########################################################
## Hartree-Fock calculation ###############################

## initial condition
## it is important to change 'magzero' in squad.in while searching for polarized solutions
muzero = 0.2 if magzero == 0.0 else -0.05 ## initial value of mu, negative for polarized case
[n,m,mu] = [0.5,magzero,muzero]

PrintAndWrite('\nCalculating the Hartree-Fock solution:',logfname)
[nold,mold,muold] = [1e5,1e5,1e5]
N0_A = sp.array([[m+n,mu],[mu,m+1.0-n]])
PrintAndWrite('\n    iter\tn\t\tm\t\tmu',logfname)
t = time()
niter = 0
while any([sp.fabs(n-nold) > eps_hf,sp.fabs(m-mold) > eps_hf,sp.fabs(mu-muold) > eps_hf]):
	Nold_A = N0_A
	[nold,mold,muold] = Occupation(Nold_A)
	G0_iw = GFzero(N0_A,int(0.8*NMats),NMats,6)
	N0_A = np.real(G0_iw.density())
	[n,m,mu] = Occupation(N0_A)
	N0_A = alpha_hf*N0_A+(1.0-alpha_hf)*Nold_A
	niter += 1
	PrintAndWrite('  {0: 3d}\t{1: .8f}\t{2: .8f}\t{3: .8f}'.format(niter,n,m,mu),logfname)
	if niter > nmax_iter:
		PrintAndWrite('HF: no convergence after '+str(nmax_iter)+' iterations, exit.',logfname)
		exit(1)
PrintAndWrite('\n  HF cycle converged, dn = {0: .8f}, dm = {1: .8f}, dmu = {2: .8f}'\
.format(n-nold,m-mold,mu-muold),logfname)

PrintAndWrite('  ...done in {0: 3d} seconds.'.format(int(time()-t)),logfname)

WriteG_iw(G0_iw,'gwhf')

NZero_A = np.real(G0_iw.density())
PrintAndWrite('\n  Occupation matrix from G_HF(iw):',logfname)
WriteMatrix(NZero_A)

G0_real = PadeContinuation(G0_iw)
WriteG_real(G0_real,'grhf')

PrintAndWrite('\n{0: .5f}\t{1: .5f}\t{2: .5f}\t{3: .5f}\t{4: .8f}\t{5: .8f}\t{6: .8f}\t:HF'\
.format(U,eps,h,P,n,m,mu),logfname)

###########################################################
## second order correction ################################

N_A = sp.copy(NZero_A)
niter = 0
Gint_iw = G0_iw.copy()
PrintAndWrite('Calculating the two-particle bubble...',logfname)
Chi = TwoParticleBubbleFFT(G0_iw,G0_iw)
PrintAndWrite('Calculating the kernel of the Schwinger-Dyson equation...',logfname)
Psi = GfImFreq(indices = [0], beta = beta, n_points = NMats, statistic = 'Boson')
## no idea how to do this correctly in new TRIQS, 
Psi.data[:,0,0] = U**2*(Chi['up','up'].data[:]+Chi['up','dn'].data[:])
PrintAndWrite('Calculating the 2ndPT dynamic self-energy...',logfname)
Sigma = SelfEnergy(G0_iw,Psi)
Gint_iw << inverse(inverse(G0_iw)-Sigma)	## Dyson equation

## write auxilliary files:
if WriteAux:
	PrintAndWrite('\nWriting auxilliary files...',logfname)
	WriteG_iw(Chi,'chiw')
	Chi_real = PadeContinuation(Chi)
	WriteG_real(Chi_real,'chir')
	WriteG_iw(Psi,'psiw')
	Psi_real = PadeContinuation(Psi)
	WriteG_real(Psi_real,'psir')
	WriteG_iw(Sigma,'sw')
	Sigma_real = PadeContinuation(Sigma)
	WriteG_real(Sigma_real,'sr')

## write the non-corrected Green function and parameters:
WriteG_iw(Gint_iw,'gw_nc')
N_A = np.real(Gint_iw.density())
PrintAndWrite('\n  Occupation matrix from G(iw):',logfname)
WriteMatrix(N_A)
[n,m,mu] = Occupation(N_A)

Gint_real = PadeContinuation(Gint_iw)
WriteG_real(Gint_real,'gr_nc')
PrintAndWrite('\n{0: .5f}\t{1: .5f}\t{2: .5f}\t{3: .5f}\t{4: .8f}\t{5: .8f}\t{6: .8f}\t:2PT0'\
.format(U,eps,h,P,n,m,mu),logfname)

###########################################################
## correcting the HF part of the self-energy only:
PrintAndWrite('\nCorrecting the Hartree-Fock self-energy:',logfname)
[nold,mold,muold] = [1e5,1e5,1e5]
PrintAndWrite('\n    iter\tn\tm\t\tmu',logfname)
niter = 0
while any([sp.fabs(n-nold) > eps_2nd,sp.fabs(m-mold) > eps_2nd,sp.fabs(mu-muold) > eps_2nd]):
	Nold_A = N_A
	[nold,mold,muold] = Occupation(N_A)
	G0_iw = GFzero(N_A,int(0.8*NMats),NMats,6)
	Gint_iw << inverse(inverse(G0_iw)-Sigma)	## Dyson equation
	N_A = np.real(Gint_iw.density())
	N_A = alpha_2nd*N_A+(1.0-alpha_2nd)*Nold_A
	[n,m,mu] = Occupation(N_A)
	niter += 1
	PrintAndWrite('  {0: 3d}\t{1: .8f}\t{2: .8f}\t{3: .8f}'.format(niter,n,m,mu),logfname)
	if niter > nmax_iter:
		PrintAndWrite('2ndPT: no convergence after '+str(nmax_iter)+' iterations, exit.',logfname)
		exit(1)
PrintAndWrite('\n  cycle converged, dn = {0: .8f}, dm = {1: .8f}, dmu = {2: .8f}'\
.format(n-nold,m-mold,mu-muold),logfname)
	
## write the corrected Green function and parameters:
WriteG_iw(Gint_iw,'gw')
N_A = np.real(Gint_iw.density())
PrintAndWrite('\n  Occupation matrix from G(iw):',logfname)
WriteMatrix(N_A)
[n,m,mu] = Occupation(N_A)

Gint_real = PadeContinuation(Gint_iw)
WriteG_real(Gint_real,'gr')

PrintAndWrite('\n{0: .5f}\t{1: .5f}\t{2: .5f}\t{3: .5f}\t{4: .8f}\t{5: .8f}\t{6: .8f}\t:2PT'\
.format(U,eps,h,P,n,m,mu),logfname)

PrintAndWrite(str(argv[0])+' END, '+ctime(),logfname)

## 2PT_mats.py END

