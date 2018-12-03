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
P['P']           = 0.5

P['pade_emax']       = 10.0
P['pade_NReal']      = 20000
P['pade_NMats']      = 50
P['pade_izero']      = 1e-3

## params part ############################################
if config.has_option('params','iw_cut'):
	P['iw_cut'] = float(config.get('params','iw_cut'))
if config.has_option('params','BW'):
	P['BW'] = float(config.get('params','BW'))
if config.has_option('params','eps_hf'):
	P['eps_hf'] = float(config.get('params','eps_hf'))
if config.has_option('params','eps_2nd'):
	P['eps_2nd'] = float(config.get('params','eps_2nd'))
if config.has_option('params','alpha_hf'):
	P['alpha_hf'] = float(config.get('params','alpha_hf'))
if config.has_option('params','alpha_2nd'):
	P['alpha_2nd'] = float(config.get('params','alpha_2nd'))
if config.has_option('params','magzero'):
	P['magzero'] = float(config.get('params','magzero'))
if config.has_option('params','param'):
	P['param'] = str(config.get('params','param'))

for pa in ['Delta','GammaL','GammaR','GammaN','P']:
	if config.has_option('params',pa):
		P[pa] = float(config.get('params',pa))

## pade part ##############################################
if config.has_option('pade','emax'):
	P['pade_emax'] = float(config.get('pade','emax'))
if config.has_option('pade','NReal'):
	P['pade_NReal'] = float(config.get('pade','NReal'))
if config.has_option('pade','NMats'):
	P['pade_NMats'] = float(config.get('pade','NMats'))
if config.has_option('pade','izero'):
	P['pade_izero'] = float(config.get('pade','izero'))

## params.py end ##

