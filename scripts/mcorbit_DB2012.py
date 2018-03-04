import numpy as np
import emcee
from random import random
from random import gauss
import sys
from orbit_utils import *
import orbit_probe
import copy
from time import sleep
from functools import partial

mw_pos = np.array([0.0, 0.0, 0.0])
mw_vel = np.array([0.0, 0.0, 0.0])

lmc_pos = np.array([-1.0, -41.0, -28.0])
#lmc_vel = np.array([-57.0, -226.0, 221.0])  # New solar
smc_pos = np.array([ 14.88893139, -38.08632213, -44.20286079])


lmc_vel = np.array([-77, -224, 227])  # vDM02
lmc_vel_err = np.array([8, 14, 18])

smc_vel = np.array([19.0, -153.0, 153.0])  # New Solar
smc_vel_err = np.array([18, 21, 17])

# This next set is drawn from the parameter space search

#lmc_vel = np.array([-78.000, -227.000, 209.000])
#lmc_vel_err = np.array([2., 2., 2.])  # keep close to original

#smc_vel = np.array([2.526, -163.461, 148.478])
#smc_vel_err = np.array([2., 2., 2.])  # keep close to original


LMCMASS = 1
SMCMASS = 0.3

gadget_LMC = {"rad": 3,#7.93195, #14.0013#21.4,
          "mass": LMCMASS,
          "gamma": 1.0,
          "a2": 1.53,
          "m2": 0.0,
          "b2": 0.34,
          "b1": 1.0,
          "m1": 0.0,
          "dynamical_friction": 0,
          "halo_type": 1,
          "mass_growth": 0,
          "dyn_C_uneq": 1.6*3.0/3., #1.22
          "dyn_L_uneq": 0.0,
          "dyn_alpha_uneq": 1.5,#uneq_alpha,
          "dyn_C_eq": 0.17,
          "dyn_L_eq": 0.02,
          "dyn_alpha_eq": 1.0,
          "tidal_truncation": 1,
          "tidal_radius": 50.,
          "inplace": 0,
          "tracers": 0,
          "pos": lmc_pos,
          "vel": lmc_vel
        }

gadget_SMC = {"rad": 2,#7.93195, #14.0013#21.4,
          "mass": SMCMASS,
          "gamma": 1.0,
          "a2": 2.9,
          "m2": 0.0,
          "b2": 0.65,
          "b1": 1.0,
          "m1": 0.0,
          "dynamical_friction": 0,
          "halo_type": 1,
          "mass_growth": 0,
          "dyn_C_uneq": 1.6*3.0/2., #1.22
          "dyn_L_uneq": 0.0,
          "dyn_alpha_uneq": 1.5,#uneq_alpha,
          "dyn_C_eq": 0.17,
          "dyn_L_eq": 0.02,
          "dyn_alpha_eq": 1.0,
          "tidal_truncation": 1,
          "tidal_radius": 40.,
          "inplace": 0,
          "tracers": 0,
          "pos": smc_pos,
          "vel": smc_vel
        }

MW_Mass = 173.0

gadget_MW = {"rad": 26.43,  # from rvir=175 and c=12
          "mass": MW_Mass,
          "gamma": 1.0,
          "a2": 3.5,
          "m2": 5,
          "b1": 0.7,
          "b2": 0.35,
          "m1": 0.5,
          "dynamical_friction": 1,
          "halo_type": 0,
          "mass_growth": 0,
          "dyn_C_uneq": 0.0, #1.22
          "dyn_L_uneq": np.log(3),  # set others to 0 to use this value
          "dyn_alpha_uneq": 0.0,#uneq_alpha,
          "dyn_C_eq": 0.0,
          "dyn_L_eq": 0.0,
          "dyn_alpha_eq": 0.0,
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
                  "tpast": -5.0,
                  "tfuture": 0.01,
                  "dt0": 1e-3,
                  "outputdir":"./output_stream/",
                  "dtout":0.05,
                  "variable_timesteps": 1,
                  "save_snapshot": 0,
                  "save_tracers": 0,
                  "output_tracers": 0}

stream_data = np.loadtxt('stream.pos.csv', delimiter=',')

def init_walker(params, lmc_vel=[0, 0, 0], smc_vel=[0, 0, 0]):
    fit = np.NINF
    while (fit < -2000):
        params[0] = lmc_vel[0] + lmc_vel_err[0]*gauss(0.0, 1)
        params[1] = lmc_vel[1] + lmc_vel_err[1]*gauss(0.0, 1)
        params[2] = lmc_vel[2] + lmc_vel_err[2]*gauss(0.0, 1)
        params[3] = smc_vel[0] + smc_vel_err[0]*gauss(0.0, 1)
        params[4] = smc_vel[1] + smc_vel_err[1]*gauss(0.0, 1)
        params[5] = smc_vel[2] + smc_vel_err[2]*gauss(0.0, 1)
        fit = loglikelihood(params)
        print "Params: {}, {}, {}, {}, {}, {}".format(*params), fit

    return params


def Serial_orbit(name):
    pool = get_pool(mpi=True, threads=2)

    # Set up MCMC
    nwalkers = 32  # Number of chains
    ndim = 6

    i = 0

    p0 = (np.random.rand(ndim * nwalkers).reshape((nwalkers, ndim)))

    # Initial guesses for parameters in each of the chains
    # from Pal 5 best fits
    p0 = np.array(pool.map(partial(init_walker,
                                   lmc_vel=lmc_vel,
                                   smc_vel=smc_vel), p0))

    print("Initial conditions set")

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
    while counter < 6000:
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
    init_vel = np.concatenate((lmc_vel, smc_vel))
    sigma_errors = np.concatenate((lmc_vel_err, smc_vel_err))

    for param, vel, err in zip(p, init_vel, sigma_errors):
        if np.abs(param - vel) > 3*err:
            return True

    return False


def loglikelihood(p):
    """
    Wraps pyorbits.orbit_statistics and checks for errors
    """
    params = copy.deepcopy(params_stream)

    params['galaxies']['LMC']['vel'][0] = p[0]
    params['galaxies']['LMC']['vel'][1] = p[1]
    params['galaxies']['LMC']['vel'][2] = p[2]

    params['galaxies']['SMC']['vel'][0] = p[3]
    params['galaxies']['SMC']['vel'][1] = p[4]
    params['galaxies']['SMC']['vel'][2] = p[5]

    if bad_parameters(p):
        return np.NINF

    ln_likelihood = orbit_probe.test_stream_DB2012(params, stream_data, False)
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
