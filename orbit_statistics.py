import numpy as np
from random import random
from random import gauss
import pyorbits
import copy
from pathos.multiprocessing import ProcessingPool as Pool
import sys
import os

mw_pos = np.array([0.0, 0.0, 0.0])
mw_vel = np.array([0.0, 0.0, 0.0])

lmc_pos = np.array([-1.0, -41.0, -28.0])
#lmc_vel = np.array([-57.0, -226.0, 221.0])  # New solar
lmc_vel = np.array([-77, -224, 227])  # vDM02

smc_vel = np.array([19.0, -153.0, 153.0])  # New Solar
smc_pos = np.array([ 14.88893139, -38.08632213, -44.20286079])

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

params_init = {"galaxies":{"MW": gadget_MW,
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



def run_pool(i):
    pid = os.getpid()
    params = copy.deepcopy(params_init)
    #Set randomly
    params['galaxies']['MW']['mass'] = 120+5*gauss(0.0, 1)
    params['galaxies']['LMC']['mass'] = 20+5*gauss(0.0, 1)
    params['galaxies']['SMC']['mass'] = 2+0.5*gauss(0.0, 1)

    params['galaxies']['LMC']['vel'][0] += 8*gauss(0.0, 1)
    params['galaxies']['LMC']['vel'][1] += 14*gauss(0.0, 1)
    params['galaxies']['LMC']['vel'][2] += 18*gauss(0.0, 1)
    params['galaxies']['SMC']['vel'][0] += 18*gauss(0.0, 1)
    params['galaxies']['SMC']['vel'][1] += 21*gauss(0.0, 1)
    params['galaxies']['SMC']['vel'][2] += 17*gauss(0.0, 1)
    try:
        results = pyorbits.run(params)
        mindist, maxdist, apos, peris, stripped = pyorbits.test_orbit_statistics(results,
                                                                        3, "MW", "LMC",
                                                                        min_time=-6.0)

        mindist_SMC, maxdist_SMC, apos_SMC, peris_SMC, stripped_SMC = pyorbits.test_orbit_statistics(results,
                                                                        3, "LMC", "SMC",
                                                                        min_time=-2.0)
        with open('parameter_space_%d.csv' % pid, 'a') as f:
            line = "{:g},"*18+"{:g}"
            line = line.format(params['galaxies']['MW']['mass'],
                               params['galaxies']['LMC']['mass'],
                               params['galaxies']['SMC']['mass'],
                               params['galaxies']['LMC']['vel'][0],
                               params['galaxies']['LMC']['vel'][1],
                               params['galaxies']['LMC']['vel'][2],
                               params['galaxies']['SMC']['vel'][0],
                               params['galaxies']['SMC']['vel'][1],
                               params['galaxies']['SMC']['vel'][2],
                               mindist, maxdist, apos, peris, stripped,
                               mindist_SMC, maxdist_SMC, apos_SMC, peris_SMC, stripped_SMC)
            f.write(line+"\n")

    except RuntimeError, e:
        with open('parameter_space_%d.csv' % pid, 'a') as f:
            line = "{:g},"*18+"{:g}"
            line = line.format(params['galaxies']['MW']['mass'],
                               params['galaxies']['LMC']['mass'],
                               params['galaxies']['SMC']['mass'],
                               params['galaxies']['LMC']['vel'][0],
                               params['galaxies']['LMC']['vel'][1],
                               params['galaxies']['LMC']['vel'][2],
                               params['galaxies']['SMC']['vel'][0],
                               params['galaxies']['SMC']['vel'][1],
                               params['galaxies']['SMC']['vel'][2],
                               np.nan, np.nan, np.nan, np.nan, np.nan,
                               np.nan, np.nan, np.nan, np.nan, np.nan)
            f.write(line+"\n")
    except KeyError, e:
        return None


def main():

    N_TRIALS = 10
    N_Processors = 4
    pool = Pool(N_Processors)
    name = 'future_parameter_space'

    print "Running"
    stats = pool.map(run_pool, xrange(N_TRIALS))


if __name__ == "__main__":
    main()

