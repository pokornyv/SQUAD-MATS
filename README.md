SQUAD-MATS
==========
#### Description:
**SQUAD-MATS** is a python2 code to calculate properties of a single-level quantum dot connected to two BCS superconducting (sc) leads
within the second-order perturbation theory (2ndPT) described in  [Žonda2015, Žonda2016]. It is based on TRIQS 1.4 libraries 
[Parcollet2015] and uses Matsubara (imaginary-frequency) formalism. The spectral functions are obtained using simple Pade
approximation. The results are being tested against the real-frequency implementation (https://github.com/pokornyv/SQUAD).  

This code is in very early stage of development and just yet being tested for simplest setups.  

SQUAD-MATS is a free software distributed under the GPL license.  

#### Project homepage:
https://github.com/pokornyv/SQUAD-MATS  

#### List of files:
- *2ndPT_mats* - main code  
- *matslib.py* - library of functions  
- *README.md* - this file  

#### TODO list:
- [x] add magnetic field dependence  
- [ ] parameters should be read from external file  
- [ ] test the code out of half-filling with magnetic field  

#### References:
- M. Žonda, V. Pokorný, V. Janiš and T. Novotný, Phys. Rev. B **93**, 024523 (2016).  
- M. Žonda, V. Pokorný, V. Janiš and T. Novotný, Sci. Rep. **5**, 8821 (2015).  
- O. Parcollet, M. Ferrero, T. Ayral, H. Hafermann, I. Krivenko, L. Messio, and P. Seth, *Comput. Phys. Commun.* **196**, 398 (2015).  

