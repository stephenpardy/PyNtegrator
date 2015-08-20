import _orbit
import numpy as np
import emcee
from math import *
from random import random
from random import gauss
import sys
from mcorbit_utils import *

n_OD = []
n_vr = []


def Serial_orbit():
    pool = get_pool(mpi=False, threads=2)

    if __name__ == '__main__':
        i = len(sys.argv)
        print("%i arguments passed to python:\n" % i)
        for j in np.arange(i):
            print(j, sys.argv[j])
        name = sys.argv[1]
        ICS = sys.argv[2]
        datafile = sys.argv[3]
        newspaper = 0  # suppress output to temporary file for performance
    else:
        name = ""
        newspaper = 1  # write each model to a temporary file for plotting

    n_OD = np.loadtxt(datafile+"_OD").astype('double')
    n_VR = np.loadtxt(datafile+"_VR").astype('double')

    # Set up MCMC
    nwalkers = 16  # Number of chains
    # Number of parameters: M, mu_alphacosdelta, mu_delta, Mass_DM, qz, a3,
    # dsun, mloss_rate
    ndim = 15

    i = 0

    p0 = (np.random.rand(ndim * nwalkers).reshape((nwalkers, ndim)))

    IC = np.loadtxt(ICS)

    # Initial guesses for parameters in each of the chains
    # from Pal 5 best fits
    while i < nwalkers:
        print("\nWalker: %i" % i)

        p0[i][0] = IC[0, 0] + IC[1, 0] * random()
#        p0[i][1] = IC[0, 1] + IC[1, 1] * random()
#        p0[i][2] = IC[0, 2] + IC[1, 2] * random()
        p0[i][1] = gauss(IC[0, 1], IC[1, 1])
        p0[i][2] = gauss(IC[0, 2], IC[1, 2])
        p0[i][3] = IC[0, 3] + IC[1, 3] * random()
        p0[i][4] = IC[0, 4] + IC[1, 4] * random()
        p0[i][5] = IC[0, 5] + IC[1, 5] * random()
        p0[i][6] = IC[0, 6] + IC[1, 6] * random()
        p0[i][7] = IC[0, 7] + IC[1, 7] * random()
        p0[i][8] = IC[0, 8] + IC[1, 8] * random()
        p0[i][9] = IC[0, 9] + IC[1, 9] * random()
        p0[i][10] = IC[0, 10] + IC[1, 10] * random()
        p0[i][11] = IC[0, 11] + IC[1, 11] * random()
        p0[i][12] = IC[0, 12] + IC[1, 12] * random()
        p0[i][13] = IC[0, 13] + IC[1, 13] * random()
        p0[i][14] = IC[0, 14] + IC[1, 14] * random()

        mass_cluster = p0[i][0]
        pm_mu_delta = p0[i][1]
        pm_mu_alphacosdelta = p0[i][2]
        mass_halo = p0[i][3]
        distance_cluster = p0[i][4]
        mass_loss_rate = p0[i][5]
        q_halo = p0[i][6]
        r_halo = p0[i][7]
        tpast = p0[i][8]
        sigma_x = p0[i][9]
        sigma_v = p0[i][10]
        sigma_vx = p0[i][11]
        sigma_mu = p0[i][12]
        rgalsun = p0[i][13]
        vLSR = p0[i][14]

        fit = _orbit.orbit(mass_cluster,
                           pm_mu_delta,
                           pm_mu_alphacosdelta,
                           mass_halo,
                           distance_cluster,
                           mass_loss_rate,
                           q_halo,
                           r_halo,
                           tpast,
                           rgalsun,
                           vLSR,
                           sigma_x,
                           sigma_v,
                           sigma_vx,
                           sigma_mu,
                           newspaper,
                           n_OD,
                           n_VR)
        if fit > -2000.:
            i += 1

    sampler = emcee.EnsembleSampler(nwalkers,
                                    ndim,
                                    loglikelihood,
                                    args=[],
                                    pool=pool)

        # Burn in
    pos, prob, state = sampler.run_mcmc(p0, 200)

    print("Burn-in complete")

    # Fill with median values instead!
    medians = 1.0 * np.arange(ndim)
    flatchain_temp = sampler.chain.copy()
    flatchain_temp.shape = -1, ndim

    for i in range(ndim):
        medians[i] = np.median(flatchain_temp[:, i])

    print(medians)

    #pos = (np.random.rand(ndim * nwalkers).reshape((nwalkers, ndim)))
        # Rescale random positions to 1% around the best position
    #for i in range(ndim):
    #    pos[:, i] = ((pos[:, i]-0.5)*0.01+1.)*medians[i]
    sampler.reset()

    # MCMC
    counter = 0
    while counter < 1024:
        counter = counter + 1
        pos, prob, state = sampler.run_mcmc(pos,
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
                 sampler.flatchain[:, 5],
                 sampler.flatchain[:, 6],
                 sampler.flatchain[:, 7],
                 sampler.flatchain[:, 8],
                 sampler.flatchain[:, 9],
                 sampler.flatchain[:, 10],
                 sampler.flatchain[:, 11],
                 sampler.flatchain[:, 12],
                 sampler.flatchain[:, 13],
                 sampler.flatchain[:, 14],
                 fout="chain%s.dat" % name)

        # some diagnostics
        # acceptance fraction
        acc_fra = sampler.acceptance_fraction
        print("\n\nAcceptance fraction (Ntemps x Nwalkers)\n", acc_fra)


def loglikelihood(x):

    global n_VR, n_OD

    if __name__ == '__main__':
        newspaper = 0  # suppress output to temporary file for performance
    else:
        newspaper = 1  # write each model to a temporary file for plotting

    Supersmall = -1.e50

    # Returns probability of input initial velocities
    mass_cluster = x[0]
    pm_mu_delta = x[1]
    pm_mu_alphacosdelta = x[2]
    mass_halo = x[3]
    distance_cluster = x[4]
    mass_loss_rate = x[5]
    q_halo = x[6]
    r_halo = x[7]
    tpast = x[8]
    sigma_x = x[9]
    sigma_v = x[10]
    sigma_vx = x[11]
    sigma_mu = x[12]
    rgalsun = x[13]
    vLSR = x[14]
    # Generate streakline model and determine log-likelihood
    #Upper and lower bounds
    if mass_cluster > 250000.0 or mass_cluster <= 25000.0:
        fit = Supersmall
    elif pm_mu_delta < -0.004 or pm_mu_delta > 0.005:
        fit = Supersmall
    elif pm_mu_alphacosdelta < -0.007 or pm_mu_alphacosdelta > 0.001:
        fit = Supersmall
    elif unbound(pm_mu_alphacosdelta, pm_mu_delta):
        fit = Supersmall
# leave unconstrained for the LM2010 halo
#    elif mass_halo > 5e12 or mass_halo < 3.e11:
#        fit = Supersmall
    elif distance_cluster < 10.0 or distance_cluster > 1883.0:
        fit = Supersmall
    elif mass_loss_rate < 0.0 or mass_loss_rate > 100.0:
        fit = Supersmall
    elif q_halo < 0.2 or q_halo > 1.8:
        fit = Supersmall
    elif r_halo > 100000.0 or r_halo < 1000.0:
        fit = Supersmall
    elif tpast > -500.0 or tpast < -10000.0:
        fit = Supersmall
    elif sigma_x > 100.0 or sigma_x < 0.0:
        fit = Supersmall
    elif sigma_v > 100.0 or sigma_v < 0.0:
        fit = Supersmall
    elif sigma_vx > 100.0 or sigma_vx < 0.0:
        fit = Supersmall
    elif rgalsun > 9000.0 or rgalsun < 7500.0:
        fit = Supersmall
    elif vLSR > 280.0 or vLSR < 200.0:
        fit = Supersmall
    elif sigma_mu > 100.0 or sigma_mu < 0.0:
        fit = Supersmall
    else:
        fit = _orbit.orbit(mass_cluster,
                           pm_mu_delta,
                           pm_mu_alphacosdelta,
                           mass_halo,
                           distance_cluster,
                           mass_loss_rate,
                           q_halo,
                           r_halo,
                           tpast,
                           rgalsun,
                           vLSR,
                           sigma_x,
                           sigma_v,
                           sigma_vx,
                           sigma_mu,
                           newspaper,
                           n_OD,
                           n_VR)

    if np.any(np.isnan(fit)):
        fit = Supersmall

    if not __name__ == '__main__' and fit > Supersmall + 1.0:
        plt.ion()
        draw_last_model()

    return fit


# Use a flat prior
def logp(x):
    return 0.0


if __name__ == '__main__':
    Serial_orbit()
