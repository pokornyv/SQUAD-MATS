## second-order PT solver for the superconducting quantum dot
## Matsubara frequency formulation
## uses TRIQS 1.4
## Vladislav Pokorny; 2018; pokornyv@fzu.cz

from os import listdir
from ConfigParser import SafeConfigParser

config = SafeConfigParser()

if 'squad.in' not in listdir('.'): 
	print('Parameter file squad.in missing. Exit.')
	exit(1)

config.read('squad.in')

P = {}

## default values #########################################

P['iw_cut']      = 100.0	
P['BW']          = 1000.0
P['eps_hf']      = 1e-4
P['eps_2nd']     = 1e-4
P['alpha_hf']    = 0.2
P['alpha_2nd']   = 0.2
P['magzero']     = 0.0
P['param']       = 'U'

P['Delta']       = 1.0
P['GammaL']      = 0.5
P['GammaR']      = 0.5
P['GammaN']      = 0.0

P['WriteAuxFiles'] = 1

P['pade_emax']       = 10.0
P['pade_NReal']      = 20000
P['pade_NMats']      = 100
P['pade_izero']      = 1e-3

## params part ############################################
for pa in ['iw_cut','BW','eps_hf','eps_2nd','alpha_hf','alpha_2nd','magzero']:
	if config.has_option('params',pa):
		P[pa] = float(config.get('params',pa))

if config.has_option('params','param'):
	P['param'] = str(config.get('params','param'))

for pa in ['Delta','GammaL','GammaR','GammaN']:
	if config.has_option('params',pa):
		P[pa] = float(config.get('params',pa))

if config.has_option('params','WriteAuxFiles'):
	P['WriteAux'] = bool(int(config.get('params','WriteAuxFiles')))

## pade part ##############################################
if config.has_option('pade','emax'):
	P['pade_emax'] = float(config.get('pade','emax'))
if config.has_option('pade','NRealPoints'):
	P['pade_NReal'] = int(config.get('pade','NRealPoints'))
if config.has_option('pade','NMatsubara'):
	P['pade_NMats'] = int(config.get('pade','NMatsubara'))
if config.has_option('pade','izero'):
	P['pade_izero'] = float(config.get('pade','izero'))

## params.py end ##

