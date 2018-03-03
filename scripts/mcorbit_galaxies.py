import numpy as np
import emcee
from math import *
from random import random
from random import gauss
import sys
from orbit_utils import *
import pyorbits
import copy

mw_pos = np.array([0.0, 0.0, 0.0])
mw_vel = np.array([0.0, 0.0, 0.0])

lmc_pos = np.array([-1.0, -41.0, -28.0])
lmc_vel = np.array([-57.0, -226.0, 221.0])  # New solar

smc_vel = np.array([19.0, -153.0, 153.0])  # New Solar
smc_pos = np.array([ 14.88893139, -38.08632213, -44.20286079])

LMCMASS = 9.9

gadget_LMC = {"rad": 14.0013,#7.93195, #14.0013#21.4,
          "mass": LMCMASS,
          "gamma": 1.0,
          "a2": 1.43,
          "m2": 0.02*LMCMASS,
          "b2": 0.312,
          "b1": 1.0,
          "m1": 0.0,
          "dynamical_friction": 1,
          "dyn_C_uneq": 1.6*3.0/14., #1.22
          "dyn_L_uneq": 0.0,
          "dyn_alpha_uneq": 1.5,#uneq_alpha,
          "dyn_C_eq": 0.17,
          "dyn_L_eq": 0.02,
          "dyn_alpha_eq": 1.0,
          "tidal_truncation": 0,
          "inplace": 0,
          "tracers": 0,
          "pos": lmc_pos,
          "vel": lmc_vel
        }


gadget_MW = {"rad": 29.8,
          "mass": 120.0,
          "gamma": 1.0,
          "a2": 3.5,
          "m2": 4.8,
          "b1": 0.7,
          "b2": 0.53,
          "m1": 1.2,
          "dynamical_friction": 1,
          "dyn_C_uneq": 1.6*3.0/29.8, #1.22
          "dyn_L_uneq": 0.0,
          "dyn_alpha_uneq": 1.5,#uneq_alpha,
          "dyn_C_eq": 0.17,
          "dyn_L_eq": 0.02,
          "dyn_alpha_eq": 1.0,
          "tidal_truncation": 0,
          "inplace": 0,
          "tracers": 0,
          "pos": np.array([0.0, 0.0, 0.0]),
          "vel": np.array([0.0, 0.0, 0.0])
        }


gadget_SMC = {"rad": 2.99579,
              "mass": 1.9,
              "gamma": 1.0,
              "a2": 0.279,
              "m2": 0.1,
              "b2": 0.062,
              "b1": 1.0,
              "m1": 0.0,
              "dynamical_friction": 0,
              "tidal_truncation": 0,
              "inplace": 0,
              "tracers": 0,
              "pos": smc_pos,
              "vel": smc_vel
              }


params_stream = {"galaxies":{"MW": gadget_MW,
                             "LMC": gadget_LMC,
                             "SMC": gadget_SMC
                             },
                  "pos_err": 100.0,  # 0.1 kpc randomly choosen for now
                  "vel_err": 20.0,  # km/s rough avg. from 3rd epoch data
                  "tpast": -8.0,
                  "tfuture": 0.0,
                  "dt0": 1e-3,
                  "outputdir":"./output_stream/",
                  "dtout":0.05,
                  "save_snapshot": 0,
                  "save_tracers": 0,
                  "output_tracers": 0}

stream_data = np.loadtxt('stream.pos.csv', delimiter=',')


def Serial_orbit(name):
    pool = get_pool(mpi=False, threads=8)

    # Set up MCMC
    nwalkers = 32  # Number of chains
    # Number of parameters: M, mu_alphacosdelta, mu_delta, Mass_DM, qz, a3,
    # dsun, mloss_rate
    ndim = 5

    i = 0

    p0 = (np.random.rand(ndim * nwalkers).reshape((nwalkers, ndim)))

    # Initial guesses for parameters in each of the chains
    # from Pal 5 best fits
    while i < nwalkers:
        print("\nWalker: %i" % i)

        p0[i][0] = 120+5*gauss(1.0, 1)
        p0[i][1] = 9.9+5*gauss(-1.0, 1)

        p0[i][2] = lmc_vel[0] + 13*gauss(0.0, 1)
        p0[i][3] = lmc_vel[1] + 15*gauss(0.0, 1)
        p0[i][4] = lmc_vel[2] + 19*gauss(0.0, 1)


        fit = loglikelihood(p0[i])
        if fit > -2000.:
            i += 1

        sampler = emcee.EnsembleSampler(nwalkers,
                                        ndim,
                                        loglikelihood,
                                        args=[],
                                        pool=pool)

        # Burn in
    likepos, prob, state = sampler.run_mcmc(p0, 500)

    print("Burn-in complete")

    # Fill with median values instead!
    medians = 1.0 * np.arange(ndim)
    flatchain_temp = sampler.chain.copy()
    flatchain_temp.shape = -1, ndim

    for i in range(ndim):
        medians[i] = np.median(flatchain_temp[:, i])

    print(medians)

    sampler.reset()

    # MCMC
    counter = 0
    while counter < 3000:
        counter = counter + 1
        likepos, prob, state = sampler.run_mcmc(likepos,
                                                1,
                                                rstate0=state,
                                                lnprob0=prob)

        print("MCMC step %i complete, chain written to chain%s.dat"
              % (counter, name))

        # Print chain to file
        printcol(sampler.flatlnprobability,
                 sampler.flatchain[:, 0],
                 sampler.flatchain[:, 1],
                 sampler.flatchain[:, 2],
                 sampler.flatchain[:, 3],
                 sampler.flatchain[:, 4],
                 fout="chain%s.dat" % name)

        # some diagnostics
        # acceptance fraction
        acc_fra = sampler.acceptance_fraction
        print("\n\nAcceptance fraction (Ntemps x Nwalkers)\n", acc_fra)

def bad_parameters(p):
    if (p[0] < 80) or (p[0] > 140):
        return True
    if (p[1] < 8) or (p[1] > 23):
        return True
    if (p[2] < -100) or (p[2] > 0):
        return True
    if (p[3] < -270) or (p[3] > -200):
        return True
    if (p[4] < 180) or (p[4] > 270):
        return True

    return False


def loglikelihood(p):
    """
    Wraps pyorbits.orbit_statistics and checks for errors
    """
    params = copy.deepcopy(params_stream)
    params['galaxies']['MW']['mass'] = p[0]
    params['galaxies']['LMC']['mass'] = p[1]

    params['galaxies']['LMC']['vel'][0] = p[2]
    params['galaxies']['LMC']['vel'][1] = p[3]
    params['galaxies']['LMC']['vel'][2] = p[4]

    if bad_parameters(p):
        return np.NINF

    ln_likelihood = pyorbits.test_stream(params, stream_data, False)
    if ln_likelihood is None:
        return np.NINF
    else:
        return ln_likelihood


if __name__ == '__main__':
    Serial_orbit(sys.argv[1])
