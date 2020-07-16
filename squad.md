# Description of the parameter file *squad.in*

#### [params] section  
*Ndot* - number of quantum dots, only 1 and 2 are implemented  
*iw_cut* - cutoff in Matsubara frequencies, default: 200.0  
*BW* - half bandwidth, default: 100.0  
*eps_hf* - convergence criterium for Hartree-Fock densities, default: 1e-4  
*eps_2nd* - convergence criterium for final densities, default: 1e-4  
*alpha_hf* - mixing parameter for Hartree-Fock calculation, default: 0.5  
*alpha_2nd* - mixing parameter for final density, default: 0.5  
*magzero* - initial value of the magnetization, useful while searching for polarized solutions, default: 0.0  
*param* - independent parameter for output, e.g. U, eps, Phi...  

*Delta* - superconducting gap, default: 1.0 (energy unit)  
*GammaL* - coupling to the left electrode  
*GammaR* - coupling to the right electrode  
*GammaN*  - coupling to the normal electrode (optional)  

*WriteAuxFiles* - 0/1 switch  - write auxiliary files (bubbles, self-energy etc...)  

*run_FSC* - run fully self-consistent cycle (obsolete, not implemented in the updated version)  

#### [pade] section:  
*emax* - maximum of the real frequency axis, default: 10  
*NRealPoints* - number of real frequency slices, default: 20000  
*NMatsubara* - number of Matsubara frequencies used in the continuation, default: 50  
*izero* - small imaginary part of the real frequency to guarantee analytic properties, default: 1e-3  

