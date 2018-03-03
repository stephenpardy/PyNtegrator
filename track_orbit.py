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

# Obs:

lmc_pos = np.array([-1.0, -41.0, -28.0])
smc_pos = np.array([ 14.88893139, -38.08632213, -44.20286079])

lmc_vel = np.array([-78.9574267558, -227.486025797, 207.835618258])
smc_vel = np.array([0.263575776662, -163.219855122, 145.413236082])

# Besla + 2012:

lmc_pos = np.array([48, 198, -85])
lmc_vel = np.array([-17, -160, -29])

smc_pos = np.array([56, 193, -90])
smc_vel = np.array([-51, -289, 88])


#Besla + 2012

LMCMASS = 17.64
SMCMASS = 1.995

gadget_LMC = {"rad": 21.4,
          "mass": LMCMASS,
          "gamma": 1.0,
          "a2": 1.53,
          "m2": 0.36,
          "b2": 0.34,
          "b1": 1.0,
          "m1": 0.0,
          "dynamical_friction": 0,
          "mass_growth": 0,
          "dyn_C_uneq": 1.6*3.0/21.4, #1.22
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

gadget_SMC = {"rad": 7.3,
          "mass": SMCMASS,
          "gamma": 1.0,
          "a2": 2.9,
          "m2": 0.105,
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

gadget_MW = {"rad": 35.66,
          "mass": 150.0,
          "gamma": 1.0,
          "a2": 3.5,
          "m2": 0.0,
          "b1": 0.7,
          "b2": 0.53,
          "m1": 0.0,
          "dynamical_friction": 0,
          "mass_growth": 0,
          "dyn_C_uneq": 1.6*3.0/35.66, #1.22
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
                  "tpast": 0.0,
                  "tfuture": 0.98,
                  "dt0": 1e-3,
                  "outputdir":"./output_stream/",
                  "dtout":0.05,
                  "variable_timesteps": 1,
                  "save_snapshot": 0,
                  "save_tracers": 0,
                  "output_tracers": 0}




snaps = pyorbits.run(params_stream)

#obs
lmc_pos_obs = np.array([-1.0, -41.0, -28.0])
smc_pos_obs = np.array([ 14.88893139, -38.08632213, -44.20286079])

lmc_vel_obs = np.array([-78.9574267558, -227.486025797, 207.835618258])
smc_vel_obs = np.array([19.0, -153.0, 153.0])
#model 2 at R200
lmc_pos = np.array([48, 198, -85])
lmc_vel = np.array([-17, -160, -29])

smc_pos = np.array([56, 193, -90])
smc_vel = np.array([-51, -289, 88])
#model 2 now
lmc_pos_now = np.array([-1, -42, -26])
lmc_vel_now = np.array([-82, -263, 249])

smc_pos_now = np.array([6, -39, -35])
smc_vel_now = np.array([-66, -258, 198])


#lmc_positions = np.empty((len(snaps), 3))
#for i in xrange(len(snaps)):
#  print snaps[i][1]['pos'], lmc_pos, np.linalg.norm(snaps[i][1]['pos'] - lmc_pos)
#  lmc_positions[i, :] = snaps[i][1]['pos']

print snaps[-1][1]['pos'], lmc_pos_obs, lmc_pos_now
