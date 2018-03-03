import numpy as np
import emcee
from math import *
from random import random
from random import gauss
import sys
from orbit_utils import *
import pyorbits
import copy
import matplotlib.pyplot as plt

from snaptools import convert

def to_magellanic_coords(l, b, l0=278.5, lambda0=32.8610, epsilon=97.5):
    sin = lambda x: np.sin(np.radians(x))
    cos = lambda x: np.cos(np.radians(x))
    atan = np.arctan2
    asin = np.arcsin

    lam =  atan(sin(b)*sin(epsilon) + cos(b)*cos(epsilon)*sin(l - l0), (cos(b)*cos(l-l0)) )  # lambda

    term2 = sin(b)*cos(epsilon) - cos(b)*sin(epsilon)*sin(l - l0)
    Beta = asin(term2)  # Beta
    return np.degrees(lam) + lambda0, np.degrees(Beta)


mw_pos = np.array([0.0, 0.0, 0.0])
mw_vel = np.array([0.0, 0.0, 0.0])

lmc_pos = np.array([-1.0, -41.0, -28.0])
smc_pos = np.array([ 14.88893139, -38.08632213, -44.20286079])

lmc_vel = np.array([-78.9574267558, -227.486025797, 207.835618258])
smc_vel = np.array([19.0, -153.0, 153.0])


LMCMASS = 18
SMCMASS = 6

gadget_LMC = {"rad": 21.7,
          "mass": LMCMASS,
          "gamma": 1.0,
          "a2": 1.53,
          "m2": 0.02*LMCMASS,
          "b2": 0.34,
          "b1": 1.0,
          "m1": 0.0,
          "dynamical_friction": 1,
          "mass_growth": 0,
          "dyn_C_uneq": 1.6*3.0/21.7,
          "dyn_L_uneq": 0.0,
          "dyn_alpha_uneq": 1.5,#uneq_alpha,
          "dyn_C_eq": 0.17,
          "dyn_L_eq": 0.02,
          "dyn_alpha_eq": 1.0,
          "tidal_truncation": 0,
          "tidal_radius": 30.,
          "inplace": 0,
          "tracers": 0,
          "pos": lmc_pos,
          "vel": lmc_vel
        }

gadget_SMC = {"rad": 7.3,
          "mass": SMCMASS,
          "gamma": 1.0,
          "a2": 2.9,
          "m2": 0.05*SMCMASS,
          "b2": 0.65,
          "b1": 1.0,
          "m1": 0.0,
          "dynamical_friction": 0,
          "mass_growth": 0,
          "dyn_C_uneq": 1.6*3.0/6., #1.22
          "dyn_L_uneq": 0.0,
          "dyn_alpha_uneq": 1.5,#uneq_alpha,
          "dyn_C_eq": 0.17,
          "dyn_L_eq": 0.02,
          "dyn_alpha_eq": 1.0,
          "tidal_truncation": 0,
          "tidal_radius": 20.,
          "inplace": 0,
          "tracers": 0,
          "pos": smc_pos,
          "vel": smc_vel
        }

gadget_MW = {"rad": 29.2,
          "mass": 150.0,
          "gamma": 1.0,
          "a2": 3.5,
          "m2": 0.04*120,
          "b1": 0.7,
          "b2": 0.53,
          "m1": 0.01*120,
          "dynamical_friction": 0,
          "mass_growth": 0,
          "dyn_C_uneq": 1.6*3.0/29.2, #1.22
          "dyn_L_uneq": 0.0,
          "dyn_alpha_uneq": 1.5,#uneq_alpha,
          "dyn_C_eq": 0.17,
          "dyn_L_eq": 0.02,
          "dyn_alpha_eq": 1.0,
          "tidal_truncation": 0,
          "tidal_radius": np.nan,
          "inplace": 1,
          "tracers": 0,
          "pos": np.array([0.0, 0.0, 0.0]),
          "vel": np.array([0.0, 0.0, 0.0])
        }

params_stream = {"galaxies":{"MW": gadget_MW,
                             "LMC": gadget_LMC,
                             "SMC": gadget_SMC,
                             },
                  "pos_err": 0.1,  # degrees
                  "vel_err": 20.0,  # km/s rough avg. from 3rd epoch data
                  "tpast": -1.0,
                  "tfuture": 0.0,
                  "dt0": 1e-3,
                  "outputdir":"./output_stream/",
                  "dtout":0.05,
                  "variable_timesteps": 0,
                  "save_snapshot": 0,
                  "save_tracers": 0,
                  "output_tracers": 0}



results = pyorbits.run(params_stream)
print results[-1][0]['t']
print "lmc_pos = np.array([{:8.5f}, {:8.5f}, {:8.5f}])".format(*(np.array(results[-1][1]['pos']) -
                                                                  np.array(results[-1][0]['pos'])))

print "lmc_vel = np.array([{:8.5f}, {:8.5f}, {:8.5f}])".format(*(np.array(results[-1][1]['vel']) -
                                                                  np.array(results[-1][0]['vel'])))

print "smc_pos = np.array([{:8.5f}, {:8.5f}, {:8.5f}])".format(*(np.array(results[-1][2]['pos']) -
                                                                  np.array(results[-1][0]['pos'])))

print "smc_vel = np.array([{:8.5f}, {:8.5f}, {:8.5f}])".format(*(np.array(results[-1][2]['vel']) -
                                                                  np.array(results[-1][0]['vel'])))





if True:
	fig, axis = plt.subplots(1, figsize=(10, 5))

	lmc_wMW = np.zeros((len(results), 3))
	smc_wMW = np.zeros((len(results), 3))
	times = np.zeros((len(results)))

	for i, result in enumerate(results):
		lmc_wMW[i, :] = result[1]['pos']
		smc_wMW[i, :] = result[2]['pos'] 
		times[i] = result[0]['t']

	#lb_lmc, _ = convert.vcart_to_gal(lmc_wMW.astype('float_')*1000., np.zeros_like(lmc_wMW).astype('float_'))
	#lb_smc, _ = convert.vcart_to_gal(smc_wMW.astype('float_')*1000., np.zeros_like(smc_wMW).astype('float_'))

	#lam_lmc, Beta_lmc = to_magellanic_coords(lb_lmc[:, 0], lb_lmc[:, 1], l0=278.5,
	#                                         lambda0=32.8610, epsilon=97.5)
	#lam_smc, Beta_smc = to_magellanic_coords(lb_smc[:, 0], lb_smc[:, 1], l0=278.5,
	#                                         lambda0=32.8610, epsilon=97.5)

	#axis.plot(lam_lmc, Beta_lmc, color='black')
	#axis.plot(lam_smc, Beta_smc, color='black', linestyle='--')
	#axis.set_xlim([100, -150])
	#axis.set_ylim([-50, 50])

	#axis.set_xlabel('Magellanic Longitude', fontsize=15)
	#axis.set_ylabel('Magellanic Latitude', fontsize=15)

	axis.plot(times, np.linalg.norm(lmc_wMW, axis=1), color='green')
	axis.plot(times, np.linalg.norm(smc_wMW, axis=1), color='black')

	plt.tight_layout()
	plt.show()