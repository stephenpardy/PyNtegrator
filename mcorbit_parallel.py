import _orbit
import numpy as np
#from subprocess import call
import emcee
from math import *
from random import random
from random import gauss
from myutils import *
import sys  # alternatively use argparse package
from mcorbit_data import *
from mcorbit_utils import *


#Parallel tempering MCMC
#Is this better??
def Parallel_orbit():

    pool = get_pool(mpi=True, threads=4)

    if __name__ == '__main__':
        i = len(sys.argv)
        print("%i arguments passed to python:\n" % i)
        for j in np.arange(i):
            print(j, sys.argv[j])
        name = sys.argv[1]
        ICS = sys.argv[2]
        newspaper = 0  # suppress output to temporary file for performance
    else:
        name = ""
        newspaper = 1  # write each model to a temporary file for plotting

    #Supersmall = -1.e50

    # Set up MCMC
    ntemps = 2      # Number of Temperatures
    nwalkers = 256  # Number of chains
                    # should be 2x free parameters (more the better!)
    # Number of parameters: M, mu_alphacosdelta, mu_delta, Mass_DM, qz, a3,
    # dsun, mloss_rate, dl, db
    ndim = 15

    j = -1

    p0 = (np.random.rand(ntemps * ndim * nwalkers)
          .reshape((ntemps, nwalkers, ndim)))

    IC = np.loadtxt(ICS)
    # Initial guesses for parameters in each of the chains
    while j < ntemps - 1:
        j += 1
        i = 0
        while i < nwalkers:
            print("\nTemperature: %i\tWalker: %i" % (j, i))
            p0[j][i][0] = IC[0, 0] + IC[1, 0] * random()
#            p0[j][i][1] = IC[0, 1] + IC[1, 1] * random()
#            p0[j][i][2] = IC[0, 2] + IC[1, 2] * random()
            p0[j][i][1] = gauss(IC[0, 1], IC[1, 1])
            p0[j][i][2] = gauss(IC[0, 2], IC[1, 2])
            p0[j][i][3] = IC[0, 3] + IC[1, 3] * random()
            p0[j][i][4] = IC[0, 4] + IC[1, 4] * random()
            p0[j][i][5] = IC[0, 5] + IC[1, 5] * random()
            p0[j][i][6] = IC[0, 6] + IC[1, 6] * random()
            p0[j][i][7] = IC[0, 7] + IC[1, 7] * random()
            p0[j][i][8] = IC[0, 8] + IC[1, 8] * random()
            p0[j][i][9] = IC[0, 9] + IC[1, 9] * random()
            p0[j][i][10] = IC[0, 10] + IC[1, 10] * random()
            p0[j][i][11] = IC[0, 11] + IC[1, 11] * random()
            p0[j][i][12] = IC[0, 12] + IC[1, 12] * random()
            p0[j][i][13] = IC[0, 13] + IC[1, 13] * random()
            p0[j][i][14] = IC[0, 14] + IC[1, 14] * random()

            mass_cluster = p0[j][i][0]
            pm_mu_delta = p0[j][i][1]
            pm_mu_alphacosdelta = p0[j][i][2]
            mass_halo = p0[j][i][3]
            distance_cluster = p0[j][i][4]
            mass_loss_rate = p0[j][i][5]
            q_halo = p0[j][i][6]
            r_halo = p0[j][i][7]
            tpast = p0[j][i][8]
            sigma_x = p0[j][i][9]
            sigma_v = p0[j][i][10]
            sigma_vx = p0[j][i][11]
            sigma_mu = p0[j][i][12]
            rgalsun = p0[j][i][13]
            vLSR = p0[j][i][14]
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
                               newspaper)
            if fit > -2000:
                i += 1

    sampler = emcee.PTSampler(ntemps,
                              nwalkers,
                              ndim,
                              loglikelihood,
                              logp,
                              pool=pool)

    # Burn in
    for pos, prob, state in sampler.sample(p0, iterations=100):
        pass

    print("Burn-in complete")

    # Fill with median values instead!
    #medians = 1.0 * np.arange(ndim)
    #flatchain_temp = sampler.chain.copy()
    #flatchain_temp.shape = ntemps, -1, ndim
#
    #for i in range(ndim):
    #    medians[i] = np.median(flatchain_temp[0, :, i])
#
    #print(medians)

    #pos = (np.random.rand(ntemps * ndim * nwalkers)
    #       .reshape((ntemps, nwalkers, ndim)))
        # Rescale random positions to 1% around the best position
    #for i in range(ndim):
    #    pos[:, :, i] = ((pos[:, :, i]-0.5)*0.01+1.)*medians[i]
    sampler.reset()

    # MCMC
    counter = 0

    while counter < 10000:
        counter = counter + 1
        for pos, prob, state in sampler.sample(pos,
                                               lnprob0=prob,
                                               lnlike0=state,
                                               iterations=1):
            pass

        print("MCMC step %i complete, chain written to chain%s.dat"
              % (counter, name))

        flatchain_temp = sampler.chain.copy()
        flatchain_temp.shape = ntemps, -1, ndim

        flatlnprobability_temp = sampler.lnprobability.copy()
        flatlnprobability_temp.shape = ntemps, -1

        # Print chains to files
        for i in range(ntemps):
            printcol(flatlnprobability_temp[i, :],
                     flatchain_temp[i, :, 0],
                     flatchain_temp[i, :, 1],
                     flatchain_temp[i, :, 2],
                     flatchain_temp[i, :, 3],
                     flatchain_temp[i, :, 4],
                     flatchain_temp[i, :, 5],
                     flatchain_temp[i, :, 6],
                     flatchain_temp[i, :, 7],
                     flatchain_temp[i, :, 8],
                     flatchain_temp[i, :, 9],
                     flatchain_temp[i, :, 10],
                     flatchain_temp[i, :, 11],
                     flatchain_temp[i, :, 12],
                     flatchain_temp[i, :, 13],
                     flatchain_temp[i, :, 14],
                     fout="chain%s_t%i.dat" % (name, i))

        # some diagnostics
        # betas
        betas = sampler.betas
        print("\n\nSequence of inverse temperatures in the ladder: ", betas)

        # temperature swap fractions
        swa_fra = sampler.tswap_acceptance_fraction
        print("\n\nAccepted temperature swap fractions for each temperature: ",
              swa_fra)

        # acceptance fraction
        acc_fra = sampler.acceptance_fraction
        print("\n\nAcceptance fraction (Ntemps x Nwalkers)\n", acc_fra)


def loglikelihood(x):

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
    elif pm_mu_alphacosdelta < -0.008 or pm_mu_alphacosdelta > 0.001:
        fit = Supersmall
    elif unbound(pm_mu_alphacosdelta, pm_mu_delta):
        fit = Supersmall
    elif mass_halo > 1.e13 or mass_halo < 1.e11:
        fit = Supersmall
    elif distance_cluster < 10.0 or distance_cluster > 20.0:
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
                           newspaper)

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
    Parallel_orbit()
