import numpy as np
import emcee
from random import random
from random import gauss
import sys
from orbit_utils import *
import pyorbits
import copy
from time import sleep

mw_pos = np.array([0.0, 0.0, 0.0])
mw_vel = np.array([0.0, 0.0, 0.0])

lmc_pos = np.array([-1.0, -41.0, -28.0])
#lmc_vel = np.array([-57.0, -226.0, 221.0])  # New solar
lmc_vel = np.array([-77, -224, 227])  # vDM02

smc_vel = np.array([19.0, -153.0, 153.0])  # New Solar
smc_pos = np.array([ 14.88893139, -38.08632213, -44.20286079])



LMCMASS = 25

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
          "mass": 2.5,
          "gamma": 1.0,
          "a2": 2.9,
          "m2": 2.5*0.05,
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
          "mass": 120.0,
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


def Serial_orbit(name):
    pool = get_pool(mpi=True, threads=2)

    # Set up MCMC
    nwalkers = 32  # Number of chains
    ndim = 6

    i = 0

    p0 = (np.random.rand(ndim * nwalkers).reshape((nwalkers, ndim)))

    # Initial guesses for parameters in each of the chains
    # from Pal 5 best fits
    while i < nwalkers:
        print("\nWalker: %i" % i)

        p0[i][0] = lmc_vel[0] + 8*gauss(0.0, 1)
        p0[i][1] = lmc_vel[1] + 14*gauss(0.0, 1)
        p0[i][2] = lmc_vel[2] + 18*gauss(0.0, 1)

        p0[i][3] = smc_vel[0] + 18*gauss(0.0, 1)
        p0[i][4] = smc_vel[1] + 21*gauss(0.0, 1)
        p0[i][5] = smc_vel[2] + 17*gauss(0.0, 1)

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
    while counter < 5000:
        counter = counter + 1
        likepos, prob, state = sampler.run_mcmc(likepos,
                                                1,
                                                rstate0=state,
                                                lnprob0=prob)

        print("MCMC step %i complete, chain written to chain%s.dat"
              % (counter, name))

        # Print chain to file
        ioloop(sampler, name)

        # some diagnostics
        # acceptance fraction
        acc_fra = sampler.acceptance_fraction
        print("\n\nAcceptance fraction (Ntemps x Nwalkers)\n", acc_fra)

    pool.close()

def bad_parameters(p):
    init_vel = np.array([-77.0, -224.0, 227.0, 19.0, -153.0, 153.0])
    sigma_errors = [8, 14, 18, 18, 21, 17]

    for param, vel, err in zip(p, init_vel, sigma_errors):
        if np.abs(param - vel) > 3*err:
            return True

    return False


def loglikelihood(p):
    """
    Wraps pyorbits.orbit_statistics and checks for errors
    """
    params = copy.deepcopy(params_stream)

    print "Params: {}, {}, {}, {}, {}, {}".format(*p)
    params['galaxies']['LMC']['vel'][0] = p[0]
    params['galaxies']['LMC']['vel'][1] = p[1]
    params['galaxies']['LMC']['vel'][2] = p[2]

    params['galaxies']['SMC']['vel'][0] = p[3]
    params['galaxies']['SMC']['vel'][1] = p[4]
    params['galaxies']['SMC']['vel'][2] = p[5]

    if bad_parameters(p):
        return np.NINF

    ln_likelihood = pyorbits.test_stream_1stpass(params, stream_data, False)
    if ln_likelihood is None:
        return np.NINF
    else:
        return ln_likelihood


def ioloop(sampler, name, counter=3, wait=5):
    """
    Try to save to disk several times if the endpoint fails. Could happen if headnode is busy.
    kwargs:
        counter: number of times to try
        wait:    number of seconds to wait
    """
    try:
        printcol(sampler.flatlnprobability,
                 sampler.flatchain[:, 0],
                 sampler.flatchain[:, 1],
                 sampler.flatchain[:, 2],
                 sampler.flatchain[:, 3],
                 sampler.flatchain[:, 4],
                 sampler.flatchain[:, 5],
                 fout="chain%s.dat" % name)
    except IOError, e:
        if counter > 0:
            sleep(wait)
            ioloop(sample, name, counter=counter-1)
        else:
            raise IOError(e)


if __name__ == '__main__':
    Serial_orbit(sys.argv[1])
