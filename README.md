SQUAD-MATS
==========
#### Description:
**SQUAD-MATS** is a python2 code to calculate properties of a single-level quantum dot connected to 
two BCS superconducting (sc) leads within the second-order perturbation theory (2ndPT) described in  
[Žonda2015, Žonda2016]. It is based on TRIQS 1.4 libraries [Parcollet2015, triqs.github.io] 
and uses Matsubara (imaginary-frequency) formalism. The spectral functions are obtained 
using simple Pade approximation. The results are being tested against the real-frequency 
implementation (github.com/pokornyv/SQUAD).  

This code is in a very early stage of development and just yet being tested for simplest setups.  

SQUAD-MATS is a free software distributed under the GPL license.  

#### Project homepage:
github.com/pokornyv/SQUAD-MATS  

#### Usage:
- set parameters in *squad.in*  
- run `python 2ndPT_mats.py <float U> <float β> <float ε> <float B>`  

#### Output files:
- *2ndPT.log* - log file  
- *gwhf*, *grhf* - Hartree-Fock Green function in Matsubaras/real frequencies  
- *gw*, *gr* - 2ndPT Green function in Matsubaras/real frequencies  
- *sw*, *sr* - 2ndPT self-energy (without the static part) in Matsubaras/real frequencies  
- *chiw*, *chir* - two-particle bubble in Matsubaras/real frequencies (for development only)  
- *psiw*, *psir* - kernel of the Schwinger-Dyson equation in Matsubaras/real frequencies (for development only)  

*grhf* and *gr* files can be labeled by the value of the parameter set in *squad.in* (*param*), e.g. *grU2.0* or *grB0.1* etc.
Useful for generating maps using *spectrum_map.py*  

Output of the HF calculation and 2ndPT calculation can be extracted from *2ndPT.log*: using 
`grep ":HF_OUT" 2ndPT.log` and `grep ":2nd_OUT" 2ndPT.log` (linux only).  

#### List of files:
- *2ndPT_mats* - main code  
- *matslib.py* - library of functions  
- *squad.in* - parameter file  
- *params.py* - script reading parameter file *squad.in*  
- *README.md* - this file  
- *spectrum_map.py* - script to generate input for color maps of the spectral function  

#### TODO list:
- [x] add magnetic field dependence  
- [x] parameters should be read from external file  
- [ ] test the code out of half-filling with magnetic field  
- [ ] add HDF5 support  

#### References:
- M. Žonda, V. Pokorný, V. Janiš and T. Novotný, Phys. Rev. B **93**, 024523 (2016).  
- M. Žonda, V. Pokorný, V. Janiš and T. Novotný, Sci. Rep. **5**, 8821 (2015).  
- O. Parcollet, M. Ferrero, T. Ayral, H. Hafermann, I. Krivenko, L. Messio, and P. Seth, *Comput. Phys. Commun.* **196**, 398 (2015).  

