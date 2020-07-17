################################################################
# SQUAD-MATS - second-order PT solver for the sc quantum dot   #
# Copyright (C) 2018-2020  Vladislav Pokorny; pokornyv@fzu.cz  #
# homepage: github.com/pokornyv/SQUAD-MATS                     #
# config.py - main config file, read inputs and prepare globals#
# method described in                                          #
#    Sci. Rep. 5, 8821 (2015).                                 #
#    Phys. Rev. B 93, 024523 (2016).                           #
################################################################

import numpy as np
import scipy as sp
from scipy.integrate import simps
from os import listdir
from sys import argv,exit
from os import getcwd
from time import time,ctime
from itertools import product

from ConfigParser import SafeConfigParser

from pytriqs.gf import *
from pytriqs.gf.descriptors import Function
from pytriqs.archive import *

config = SafeConfigParser()

cfile = 'squad.in'

if cfile not in listdir('.'): 
	print('Parameter file '+cfile+' missing. Exit.')
	exit(1)

config.read(cfile)

###########################################################
## default values #########################################

Ndot        = 1
iw_cut      = 100.0	
BW          = 1000.0
eps_hf      = 1e-4
eps_2nd     = 1e-4
alpha_hf    = 0.2
alpha_2nd   = 0.2
magzero     = 0.0
param       = 'U'

Delta       = 1.0
GammaL      = 0.5
GammaR      = 0.5
GammaN      = 0.0
GammaNL     = 0.0
GammaNR     = 0.0
tLR         = 0.0

WriteAuxFiles = True
run_FSC       = False

pade_emax       = 10.0
pade_NReal      = 20000
pade_NMats      = 100
pade_izero      = 1e-3

## params part ############################################

if config.has_option('params','Ndot'):
	Ndot      = int(config.get('params','Ndot'))
if config.has_option('params','iw_cut'):
	iw_cut    = float(config.get('params','iw_cut'))
if config.has_option('params','BW'):
	BW        = float(config.get('params','BW'))
if config.has_option('params','eps_hf'):
	eps_hf    = float(config.get('params','eps_hf'))
if config.has_option('params','eps_2nd'):
	eps_2nd   = float(config.get('params','eps_2nd'))
if config.has_option('params','alpha_hf'):
	alpha_hf  = float(config.get('params','alpha_hf'))
if config.has_option('params','alpha_2nd'):
	alpha_2nd = float(config.get('params','alpha_2nd'))
if config.has_option('params','alpha_2nd'):
	magzero   = float(config.get('params','magzero'))
if config.has_option('params','param'):
	param     = str(config.get('params','param'))

if config.has_option('params','Delta'):
	Delta   = float(config.get('params','Delta'))
if config.has_option('params','GammaL'):
	GammaL  = float(config.get('params','GammaL'))
if config.has_option('params','GammaR'):
	GammaR  = float(config.get('params','GammaR'))
if config.has_option('params','GammaN'):
	GammaN  = float(config.get('params','GammaN'))
if config.has_option('params','GammaNL'):
	GammaNL = float(config.get('params','GammaNL'))
if config.has_option('params','GammaNR'):
	GammaNR = float(config.get('params','GammaNR'))
if config.has_option('params','tLR'):
	tLR     = float(config.get('params','tLR'))

if config.has_option('params','WriteAuxFiles'):
	WriteAux = bool(int(config.get('params','WriteAuxFiles')))
if config.has_option('params','run_FSC'):
	run_FSC  = bool(int(config.get('params','run_FSC')))

## pade part ##############################################
if config.has_option('pade','emax'):
	pade_emax   = float(config.get('pade','emax'))
if config.has_option('pade','NRealPoints'):
	pade_NReal  = int(config.get('pade','NRealPoints'))
if config.has_option('pade','NMatsubara'):
	pade_NMats  = int(config.get('pade','NMatsubara'))
if config.has_option('pade','izero'):
	pade_izero  = float(config.get('pade','izero'))

###########################################################
## reading parameters from command line ###################

U      = float(argv[1])
beta   = float(argv[2])
eps    = float(argv[3])
h      = float(argv[4])
P      = float(argv[5])

Phi = P*np.pi

if Ndot == 2:
	UR = UL = U
	epsL = epsR = eps
	U2 = 0.0 ## interdot coupling

###########################################################
## output files ###########################################

## HDF5 output
hdf5file  = 'squad_mats.h5'
WriteHDF5 = False

## log file
logfname = '2PT.log'

###########################################################
## define other variables #################################

nmax_iter = 1000  ## maximum number of iterations before breaking a cycle

stars  = '*'*60
hashes = '#'*60

## parameter for output file names
if   param == 'U':   par = U
elif param == 'eps': par = eps
elif param == 'h':   par = h
elif param == 'P':   par = P
elif param == 'Phi': par = P
else:                par = U	## default

if Ndot == 1: bands_T = ['up','dn']
else:         bands_T = ['Lup','Ldn','Rup','Rdn']
NBand = len(bands_T)
params_A = [beta,U,Delta,GammaL,GammaR,GammaN,eps,P,h]
NMats = int(0.5*(beta*iw_cut/np.pi)+1.0)

## config.py END ##

