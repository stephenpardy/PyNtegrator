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
mw_pos = np.array([0.0, 0.0, 0.0])
mw_vel = np.array([0.0, 0.0, 0.0])

lmc_pos = np.array([-1.0, -41.0, -28.0])
smc_pos = np.array([ 14.88893139, -38.08632213, -44.20286079])

lmc_vel = np.array([-78.9574267558, -227.486025797, 207.835618258])
smc_vel = np.array([0.263575776662, -163.219855122, 145.413236082])



LMCMASS = 15
SMCMASS = 2.5

gadget_LMC = {"rad": 19,#7.93195, #14.0013#21.4,
          "mass": LMCMASS,
          "gamma": 1.0,
          "a2": 1.53,
          "m2": 0.05*LMCMASS,
          "b2": 0.34,
          "b1": 1.0,
          "m1": 0.0,
          "dynamical_friction": 1,
          "mass_growth": 0,
          "dyn_C_uneq": 1.6*3.0/19., #1.22
          "dyn_L_uneq": 0.0,
          "dyn_alpha_uneq": 1.5,#uneq_alpha,
          "dyn_C_eq": 0.17,
          "dyn_L_eq": 0.02,
          "dyn_alpha_eq": 1.0,
          "tidal_truncation": 1,
          "tidal_radius": 30.,
          "inplace": 0,
          "tracers": 0,
          "pos": lmc_pos,
          "vel": lmc_vel
        }

gadget_SMC = {"rad": 6,#7.93195, #14.0013#21.4,
          "mass": SMCMASS,
          "gamma": 1.0,
          "a2": 2.9,
          "m2": SMCMASS*0.1,
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
          "tidal_truncation": 1,
          "tidal_radius": 20.,
          "inplace": 0,
          "tracers": 0,
          "pos": smc_pos,
          "vel": smc_vel
        }

gadget_MW = {"rad": 29.2,
          "mass": 130.0,
          "gamma": 1.0,
          "a2": 3.5,
          "m2": 0.04*120,
          "b1": 0.7,
          "b2": 0.53,
          "m1": 0.01*120,
          "dynamical_friction": 1,
          "mass_growth": 0,
          "dyn_C_uneq": 1.6*3.0/29.2, #1.22
          "dyn_L_uneq": 0.0,
          "dyn_alpha_uneq": 1.5,#uneq_alpha,
          "dyn_C_eq": 0.17,
          "dyn_L_eq": 0.02,
          "dyn_alpha_eq": 1.0,
          "tidal_truncation": 0,
          "tidal_radius": np.nan,
          "inplace": 0,
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
                  "tpast": -7.0,
                  "tfuture": 0.01,
                  "dt0": 1e-3,
                  "outputdir":"./output_stream/",
                  "dtout":0.05,
                  "variable_timesteps": 0,
                  "save_snapshot": 0,
                  "save_tracers": 0,
                  "output_tracers": 0}


stream_data = np.loadtxt('stream.pos.csv', delimiter=',')

#stream_data = np.loadtxt('stream.pos.csv', delimiter=',')

ln_likelihood = pyorbits.test_stream_2ndpass(params_stream, stream_data, False)
print ln_likelihood


# results = pyorbits.run(params)
# print "lmc_pos = np.array([{:8.5f}, {:8.5f}, {:8.5f}])".format(*(np.array(results[-1][1]['pos']) -
#                                                                  np.array(results[-1][0]['pos'])))

# print "lmc_vel = np.array([{:8.5f}, {:8.5f}, {:8.5f}])".format(*(np.array(results[-1][1]['vel']) -
#                                                                  np.array(results[-1][0]['vel'])))

# print "smc_pos = np.array([{:8.5f}, {:8.5f}, {:8.5f}])".format(*(np.array(results[-1][2]['pos']) -
#                                                                  np.array(results[-1][0]['pos'])))

# print "smc_vel = np.array([{:8.5f}, {:8.5f}, {:8.5f}])".format(*(np.array(results[-1][2]['vel']) -
#                                                                  np.array(results[-1][0]['vel'])))


# blue_dot, = plt.plot([0], "b.")
# black_dot, = plt.plot([0], "k.")
# plt.clf()

# fig, axes = plt.subplots(1, 2, figsize=(7, 3))

# axis = axes[0]

# axis.set_xlabel('Time [Gyr]')
# axis.set_ylabel('Distance from Milky Way [kpc]')

# #axis.legend([blue_dot, black_dot], ['LMC', 'SMC'])

# for result in results:
#     axis.plot(result[1]['t'], np.sqrt(np.sum((np.array(result[1]['pos']) -
#                              np.array(result[0]['pos']))**2)), 'b.')
#     axis.plot(result[1]['t'], np.sqrt(np.sum((np.array(result[2]['pos']) -
#                              np.array(result[0]['pos']))**2)), 'k.')

# axis = axes[1]


# axis.set_xlabel('Time [Gyr]')
# axis.set_ylabel('Distance from LMC [kpc]')

# for result in results:
#     axis.plot(result[1]['t'], np.sqrt(np.sum((np.array(result[2]['pos']) -
#                              np.array(result[1]['pos']))**2)), 'k.')

# plt.tight_layout()
# plt.show()

# #print ln_likelihood
