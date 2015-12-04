import numpy as np
from random import random
from random import gauss
import pyorbits
import copy
from pathos.multiprocessing import ProcessingPool as Pool
import sys


def compute_stats(p):
    """
    Wraps pyorbits.orbit_statistics and checks for errors
    """
    stats = pyorbits.orbit_statistics(p, 'MW', 'LMC')
    print(stats)
    if stats is None:
      return np.nan, np.nan, np.nan, np.nan, np.nan

    return stats


def main():

    lmc_pos = np.array([-1.0, -41.0, -28.0])
    lmc_vel = np.array([-57.0, -226.0, 221.0])  # New solar

    smc_vel = np.array([19.0, -153.0, 153.0])  # New Solar
    smc_pos = np.array([ 14.88893139, -38.08632213, -44.20286079])

    gadget_MW = {"rad": 29.8,
              "mass": 120.0,
              "gamma": 1.0,
              "c": 9.56,
              "a2": 3.5,
              "m2": 4.8,
              "b1": 0.7,
              "b2": 0.53,
              "m1": 1.2,
              "type": 1,
              "dynamical_friction": 1,
              "dyn_C_uneq": np.nan,
              "dyn_L_uneq": 0.0,
              "dyn_alpha_uneq": 1.5,
              "dyn_C_eq": 0.17,
              "dyn_L_eq": 0.02,
              "dyn_alpha_eq": 1.0,
              "tidal_truncation": 0,
              "inplace": 0,
              "pos": np.array([0.0, 0.0, 0.0]),
              "vel": np.array([0.0, 0.0, 0.0])
            }

    mass = 9.90

    gadget_LMC = {"rad": 14.0013,#7.93195, #14.0013#21.4,
              "mass": mass,
              "gamma": 1.0,
              "c": 9.56,
              "a2": 1.43,
              "m2": 0.02*mass,
              "b2": 0.312,
              "b1": 1.0,
              "m1": 0.0,
              "type": 1,
              "dynamical_friction": 0,
              "tidal_truncation": 1,
              "inplace": 0,
              "pos": lmc_pos,
              "vel": lmc_vel
            }

    m31 = {"rad": 28.7,
           "mass": 150.0,
           "gamma": 1.0,
           "c": 12.0,
              "a2": 3.5,
              "m2": 7.0,
              "b1": 0.7,
              "b2": 0.53,
              "m1": 1.9,
              "type": 1,
              "dynamical_friction": 1,
              "dyn_C_uneq": 1.6*3.0/28.7,
              "dyn_L_uneq": 0.0,
              "dyn_alpha_uneq": 1.5,
              "dyn_C_eq": 0.17,
              "dyn_L_eq": 0.02,
              "dyn_alpha_eq": 1.0,
              "tidal_truncation": 0,
              "inplace": 0,
              "pos": np.array([-378.9, 612.7, -283.1]),
              "vel": np.array([66.1, -76.3, 45.1])
            }

    kal_smc = {"rad": 1.0,
               "mass": 0.3,
               "gamma": 1.0,
               "c": 1.0,
               "a2": 0.43405,
               "b2": 0.09586,
               "m2": 0.0,
               "b1": 1.0,
               "m1": 0.0,
               "type": 1,
               "dynamical_friction": 0,
               "tidal_truncation": 1,
               "inplace": 0,
               "pos": smc_pos,
               "vel": smc_vel
               }

    params_init = {"galaxies":{"MW": gadget_MW,
                               "M31": m31,
                               "LMC": gadget_LMC,
                               "SMC": kal_smc
                              },
                    "pos_err": 10.0,  # 0.1 kpc randomly choosen for now
                    "vel_err": 20.0,  # km/s rough avg. from 3rd epoch data
                    "tpast": 0.0,
                    "tfuture": 7.0,
                    "dt0": 1e-4,
                    "outputdir":"./output_testfric/",
                    "dtout":0.05,
                    "save_snapshot": 0}


    N_TRIALS = 40
    N_per_file = 1000
    N_Processors = 4
    pool = Pool(N_Processors)
    name = 'future_parameter_space'
    for ID in xrange(N_TRIALS//N_per_file+1):
        params = []
        print "Making parameters"
        for i in xrange(N_per_file):
            params.append(copy.deepcopy(params_init))
            #Set randomly
            params[i]['galaxies']['MW']['mass'] = 120+5*gauss(1.0, 1)
            params[i]['galaxies']['LMC']['mass'] = 9.9+5*gauss(-1.0, 1)

            params[i]['galaxies']['LMC']['vel'][0] += 13*gauss(0.0, 1)
            params[i]['galaxies']['LMC']['vel'][1] += 15*gauss(0.0, 1)
            params[i]['galaxies']['LMC']['vel'][2] += 19*gauss(0.0, 1)

            #Set by other params
            params[i]['galaxies']['MW']['dyn_C_uneq'] = 1.6*3.0/params[i]['galaxies']['MW']['rad']

        print "Running"
        stats = pool.map(compute_stats, params)

        # Print chain to file
        print("Writing out")
        with file('{:s}{:d}.output'.format(name, ID), 'w') as f:
          f.write("#min,max,apo,peri,mwmass,lmcmass,lmcvx,lmcvy,lmcvz\n")
          for i in xrange(N_per_file):
                if any(np.isnan(stats[i])):
                    continue
                f.write("{:g}, {:g}, {:g}, {:g}, {:g}, {:g}, {:g}, {:g}, {:g}, {:g}\n".format(stats[i][0],
                                                                      stats[i][1],
                                                                      stats[i][2],
                                                                      stats[i][3],
                                                                      stats[i][4],
                                                                      params[i]['galaxies']['MW']['mass'],
                                                                      params[i]['galaxies']['LMC']['mass'],
                                                                      params[i]['galaxies']['LMC']['vel'][0],
                                                                      params[i]['galaxies']['LMC']['vel'][1],
                                                                      params[i]['galaxies']['LMC']['vel'][2]))


if __name__ == "__main__":
    main()

