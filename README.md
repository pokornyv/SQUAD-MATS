SQUAD-MATS
==========
#### Description:
**SQUAD-MATS** is a python 2 code to calculate properties of a quantum dot (QD) system connected to 
two BCS superconducting (sc) leads within the second-order perturbation theory (2PT) described in
[Žonda2015, Žonda2016]. It is based on TRIQS 2.x libraries [Parcollet2015, [homepage](https://triqs.github.io)] 
and uses Matsubara (imaginary-frequency) formalism. 
The code can calculate properties for a single or a double QD [Pokorný2020].
The spectral functions are obtained using simple Pade approximation. The results for a single QD 
were tested against the real-frequency implementation 
([SQUAD](https://github.com/pokornyv/SQUAD "github.com/pokornyv/SQUAD")).  

This code is in a early stage of development and should not be used for generating production data
unless you are sure what you are doing. Note that Hartree-Fock method is known for providing false
symmetry-broken (magnetically polarized) solutions. Restricted Hartree-Fock is not yet implemented.  

SQUAD-MATS is a free software distributed under the GPL license.  

#### Project homepage:
https://github.com/pokornyv/SQUAD-MATS  

#### Usage:
- set parameters in *squad.in*  
- single dot: run `python 2PT_mats.py <U> <β> <ε> <h> <P>`
- double dot: run `python 2PT_ddot_mats.py <U> <β> <ε> <h> <P>`    
where *U* is the ineraction strength, *β* is the inverse temperature, *ε* is the
gate voltage (energy shift of levels w.r.t a half-filled dot), *h* is the magnetic field
and *P=φ/π* is the phase difference.  

#### Output files:
- *2PT.log* - log file  
- *gwhf*, *grhf* - Hartree-Fock Green function in Matsubaras/real frequencies  
- *gw*, *gr* - 2PT Green function in Matsubaras/real frequencies with corrected HF part  
- *gw_nc*, *gr_nc* - 2PT Green function in Matsubaras/real frequencies w/o the corrected HF part    
- *sw*, *sr* - 2PT self-energy (without the static part) in Matsubaras/real frequencies  
- *chiw*, *chir* - two-particle bubble in Matsubaras/real frequencies (for development only)  
- *psiw*, *psir* - kernel of the Schwinger-Dyson equation in Matsubaras/real frequencies (for development only)  

Output of the HF calculation and 2ndPT calculation can be extracted from *2PT.log*: using (linux only)  
single dot: `grep ":HF" 2PT.log`, `grep ":2PT0" 2PT.log` and `grep ":2PT" 2PT.log`.  
double dot: `grep ":HF_L" 2PT.log` and `grep ":2PT_L" 2PT.log` for left QD and
similarly `grep ":HF_R" 2PT.log` and `grep ":2PT_R" 2PT.log` for right QD.  

#### List of files:
- *2PT_mats.py* - main code for single QD  
- *2PT_ddot_mats.py* - main code for double QD  
- *matslib.py* - library of functions for single QD
- *matslib2.py* - library of functions for double QD  
- *iolib.py* - input/output library  
- *config.py* - read the parameter file *squad.in*, set default values and global variables  
- *squad.in* - parameter file  
- *README.md* - this file  
- *spectrum_map.py* - script to generate input for color maps of the spectral function for a single dot  
- *spectrum_map_ddot.py* - script to generate input for color maps of the spectral function for a double dot  

#### TODO list:
- [x] implement the double QD setup  
- [x] parameters should be read from external file  
- [x] rewrite code from TRIQS 1.x to TRIQS 2.x  
- [x] add magnetic field dependence  
- [ ] implement restricted Hartree-Fock (no magnetic polarization)
- [ ] test the symmetries out of half-filling with magnetic field  
- [ ] add full HDF5 support (write and read)  
- [ ] test the complex hybridization in case of non-zero *φ* in double-dot setup  

#### References:
- M. Žonda, V. Pokorný, V. Janiš and T. Novotný, Sci. Rep. **5**, 8821 (2015).  
- M. Žonda, V. Pokorný, V. Janiš and T. Novotný, Phys. Rev. B **93**, 024523 (2016).  
- O. Parcollet, M. Ferrero, T. Ayral, H. Hafermann, I. Krivenko, L. Messio, and P. Seth, 
*Comput. Phys. Commun.* **196**, 398 (2015).  
- V. Pokorný, M. Žonda, G. Loukeris, and T. Novotný, JPS Conf. Proc. **30**, 011002 (2020).  

