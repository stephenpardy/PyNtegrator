import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
#from subprocess import call
from math import *
from myutils import *
import triangle
from mcorbit_utils import *


def param_histo_Halo(name, folder='./'):
    mpl.rcParams.update({'font.size': 12})

    print "File name: chain%s.dat" % name

        # Plot histograms of halo parameters
    prob,\
        mass_cluster,\
        pm_mu_delta,\
        pm_mu_alphacosdelta,\
        mass_halo,\
        distance_cluster,\
        mass_loss_rate,\
        q_halo,\
        r_halo,\
        tpast,\
        sigma_x,\
        sigma_v,\
        sigma_vx,\
        sigma_mu,\
        rgalsun,\
        vLSR = np.loadtxt("chain%s.dat" % name, unpack=True)

    maxid = np.argmax(prob)
    Mh_b = mass_halo[maxid]
    rh_b = r_halo[maxid]
    q_b = q_halo[maxid]
    dM_b = mass_loss_rate[maxid]
    Mh_m = np.mean(mass_halo)
    rh_m = np.mean(r_halo)
    q_m = np.mean(q_halo)
    dM_m = np.mean(mass_loss_rate)
    Mh_md = np.median(mass_halo)
    rh_md = np.median(r_halo)
    q_md = np.median(q_halo)
    dM_md = np.median(mass_loss_rate)
    Mh_84 = np.percentile(mass_halo, 84)
    rh_84 = np.percentile(r_halo, 84)
    q_84 = np.percentile(q_halo, 84)
    dM_84 = np.percentile(mass_loss_rate, 84)
    Mh_16 = np.percentile(mass_halo, 16)
    rh_16 = np.percentile(r_halo, 16)
    q_16 = np.percentile(q_halo, 16)
    dM_16 = np.percentile(mass_loss_rate, 16)
    Mh_95 = np.percentile(mass_halo, 95)
    rh_95 = np.percentile(r_halo, 95)
    q_95 = np.percentile(q_halo, 95)
    dM_95 = np.percentile(mass_loss_rate, 95)
    Mh_05 = np.percentile(mass_halo, 5)
    rh_05 = np.percentile(r_halo, 5)
    q_05 = np.percentile(q_halo, 5)
    dM_05 = np.percentile(mass_loss_rate, 5)
    Mh_975 = np.percentile(mass_halo, 97.5)
    rh_975 = np.percentile(r_halo, 97.5)
    q_975 = np.percentile(q_halo, 97.5)
    dM_975 = np.percentile(mass_loss_rate, 97.5)
    Mh_025 = np.percentile(mass_halo, 2.5)
    rh_025 = np.percentile(r_halo, 2.5)
    q_025 = np.percentile(q_halo, 2.5)
    dM_025 = np.percentile(mass_loss_rate, 2.5)
    print("Best-fit: ", Mh_b, rh_b, q_b, dM_b)
    print("Mean: ", Mh_m, rh_m, q_m, dM_m)
    print("Median: ", Mh_md, rh_md, q_md, dM_md)
    print("up (68%): ", Mh_84, rh_84, q_84, dM_84)
    print("down (68%)", Mh_16, rh_16, q_16, dM_16)
    print("up (90%): ", Mh_95, rh_95, q_95, dM_95)
    print("down (90%)", Mh_05, rh_05, q_05, dM_05)
    print("up (95%): ", Mh_975, rh_975, q_975, dM_975)
    print("down (95%)", Mh_025, rh_025, q_025, dM_025)
    print("Mhalo = %e +%e -%e (+%e -%e)"
          % (Mh_md,
             Mh_84 - Mh_md,
             Mh_md - Mh_16,
             Mh_95 - Mh_md,
             Mh_md - Mh_05))
    print("Rhalo = %f +%f -%f (+%f -%f)"
          % (rh_md,
             rh_84 - rh_md,
             rh_md - rh_16,
             rh_95 - rh_md,
             rh_md - rh_05))
    print("qz = %f +%f -%f (+%f -%f)"
          % (q_md,
             q_84 - q_md,
             q_md - q_16,
             q_95 - q_md,
             q_md - q_05))
    print("dM = %f +%f -%f (+%f -%f)"
          % (dM_md,
             dM_84 - dM_md,
             dM_md - dM_16,
             dM_95 - dM_md,
             dM_md - dM_05))

    # Plot histograms
    plt.close()
    plt.figure()

    fig, axes = plt.subplots(2, 2, figsize=(10.0, 10.0))
    axes[0, 0].hist(mass_halo[prob > -1.e4],
                    bins=10 ** np.linspace(11, 13, 20),
                    color=mpl.colors.hex2color("#ABCEEA"),
                    histtype="stepfilled")
    #axes[0, 0].set_xlim(1.e11, 1.e13)
    axes[0, 0].set_xscale('log')
    axes[0, 0].set_xlabel("M$_{halo}$ (M$_\odot$)")
    axes[0, 0].set_ylabel("N")
    axes[0, 0].axvline(Mh_md, color='r', linewidth=1)
    axes[0, 0].axvline(Mh_m, color='b', linewidth=1)
    axes[0, 0].axvline(Mh_84, color='k', linewidth=1, linestyle='--')
    axes[0, 0].axvline(Mh_16, color='k', linewidth=1, linestyle='--')
    axes[0, 0].axvline(Mh_95, color='k', linewidth=1, linestyle=':')
    axes[0, 0].axvline(Mh_05, color='k', linewidth=1, linestyle=':')

    axes[1, 0].hist(r_halo[prob > -1.e4], 30,
                    color=mpl.colors.hex2color("#ABCEEA"),
                    histtype="stepfilled")
    #axes[1, 0].set_xlim(0, 50000)
    axes[1, 0].set_xlabel("R$_{halo}$ (pc)")
    axes[1, 0].set_ylabel("N")

    axes[1, 0].axvline(rh_md, color='r', linewidth=1)
    axes[1, 0].axvline(rh_m, color='b', linewidth=1)
    axes[1, 0].axvline(rh_84, color='k', linewidth=1, linestyle='--')
    axes[1, 0].axvline(rh_16, color='k', linewidth=1, linestyle='--')
    axes[1, 0].axvline(rh_95, color='k', linewidth=1, linestyle=':')
    axes[1, 0].axvline(rh_05, color='k', linewidth=1, linestyle=':')

    axes[1, 1].hist(q_halo[prob > -1.e4], 30,
                    color=mpl.colors.hex2color("#ABCEEA"),
                    histtype="stepfilled")
    #axes[1, 1].set_xlim(0.5, 2)
    axes[1, 1].set_xlabel("q$_z$")
    axes[1, 1].set_ylabel("N")
    axes[1, 1].axvline(q_md, color='r', linewidth=1)
    axes[1, 1].axvline(q_m, color='b', linewidth=1)
    axes[1, 1].axvline(q_84, color='k', linewidth=1, linestyle='--')
    axes[1, 1].axvline(q_16, color='k', linewidth=1, linestyle='--')
    axes[1, 1].axvline(q_95, color='k', linewidth=1, linestyle=':')
    axes[1, 1].axvline(q_05, color='k', linewidth=1, linestyle=':')

    try:
        axes[0, 1].hist(mass_loss_rate[prob > -1.e4],
                        bins=10 ** np.linspace(-2, 2, 20),
                        color=mpl.colors.hex2color("#ABCEEA"),
                        histtype="stepfilled")
    except:
        pass
    #axes[0, 1].set_xlim(1.e-2, 1.e1)
    axes[0, 1].set_xscale('log')
    axes[0, 1].set_xlabel("dM/dt (M$_\odot/Myr$)")
    axes[0, 1].set_ylabel("N")
    axes[0, 1].axvline(dM_md, color='r', linewidth=1)
    axes[0, 1].axvline(dM_m, color='b', linewidth=1)
    axes[0, 1].axvline(dM_84, color='k', linewidth=1, linestyle='--')
    axes[0, 1].axvline(dM_16, color='k', linewidth=1, linestyle='--')
    axes[0, 1].axvline(dM_95, color='k', linewidth=1, linestyle=':')
    axes[0, 1].axvline(dM_05, color='k', linewidth=1, linestyle=':')

    # Save figure
    plt.draw()
    filename = folder + 'chain%s_params_histo.halo' % name
    plt.savefig(filename + '.png', format='png', dpi=150, bbox_inches='tight')


def param_histo_acc(name,
                    streamname='stream',
                    gal_latitude=42.15,
                    gal_longitude=73.59,
                    folder='./'):
    mpl.rcParams.update({'font.size': 20})

    print "File name: chain%s.dat" % name

        # Plot histograms of acceleration parameters
    prob,\
        mass_cluster,\
        pm_mu_delta,\
        pm_mu_alphacosdelta,\
        mass_halo,\
        distance_cluster,\
        mass_loss_rate,\
        q_halo,\
        r_halo,\
        tpast,\
        sigma_x,\
        sigma_v,\
        sigma_vx,\
        sigma_mu,\
        rgalsun,\
        vLSR = np.loadtxt("chain%s.dat" % name, unpack=True)

    good_inds = [good_vals(vals) for vals in [prob,
                                              rgalsun,
                                              mass_halo,
                                              q_halo,
                                              r_halo,
                                              distance_cluster]]
    ind = good_inds[0]
    for i in xrange(1, len(good_inds)):
        ind = np.intersect1d(good_inds[i], ind)

    prob = prob[ind]
    rgalsun = rgalsun[ind]
    mass_halo = mass_halo[ind]
    q_halo = q_halo[ind]
    r_halo = r_halo[ind]
    distance_cluster = distance_cluster[ind]

    Vsun, R200, cNFW, M200, aPal5 = calculate_halo_properties(prob,
                                                              rgalsun,
                                                              mass_halo,
                                                              q_halo,
                                                              r_halo,
                                                              distance_cluster,
                                                              gal_latitude,
                                                              gal_longitude)

    maxid = np.argmax(prob)
    a_b = aPal5[maxid]
    V_b = Vsun[maxid]
    R_b = R200[maxid]
    M_b = M200[maxid]
    cNFW_b = cNFW[maxid]
    a_m = np.mean(aPal5)
    V_m = np.mean(Vsun)
    R_m = np.mean(R200)
    M_m = np.mean(M200)
    cNFW_m = np.mean(cNFW)
    a_md = np.median(aPal5)
    V_md = np.median(Vsun)
    R_md = np.median(R200)
    M_md = np.median(M200)
    cNFW_md = np.median(cNFW)
    a_84 = np.percentile(aPal5, 84)
    V_84 = np.percentile(Vsun, 84)
    R_84 = np.percentile(R200, 84)
    M_84 = np.percentile(M200, 84)
    cNFW_84 = np.percentile(cNFW, 84)
    a_16 = np.percentile(aPal5, 16)
    V_16 = np.percentile(Vsun, 16)
    R_16 = np.percentile(R200, 16)
    M_16 = np.percentile(M200, 16)
    cNFW_16 = np.percentile(cNFW, 16)
    a_95 = np.percentile(aPal5, 95)
    V_95 = np.percentile(Vsun, 95)
    R_95 = np.percentile(R200, 95)
    M_95 = np.percentile(M200, 95)
    cNFW_95 = np.percentile(cNFW, 95)
    a_05 = np.percentile(aPal5, 5)
    V_05 = np.percentile(Vsun, 5)
    R_05 = np.percentile(R200, 5)
    M_05 = np.percentile(M200, 5)
    cNFW_05 = np.percentile(cNFW, 5)
    a_975 = np.percentile(aPal5, 97.5)
    V_975 = np.percentile(Vsun, 97.5)
    R_975 = np.percentile(R200, 97.5)
    M_975 = np.percentile(M200, 97.5)
    cNFW_975 = np.percentile(cNFW, 97.5)
    a_025 = np.percentile(aPal5, 2.5)
    V_025 = np.percentile(Vsun, 2.5)
    R_025 = np.percentile(R200, 2.5)
    M_025 = np.percentile(M200, 2.5)
    cNFW_025 = np.percentile(cNFW, 2.5)
    print("Best-fit: ", a_b, V_b, R_b, M_b, cNFW_b)
    print("Mean: ", a_m, V_m, R_m, M_m, cNFW_m)
    print("Median: ", a_md, V_md, R_md, M_md, cNFW_md)
    print("up (68%): ", a_84, V_84, R_84, M_84, cNFW_84)
    print("down (68%): ", a_16, V_16, R_16, M_16, cNFW_16)
    print("up (90%): ", a_95, V_95, R_95, M_95, cNFW_95)
    print("down (90%): ", a_05, V_05, R_05, M_05, cNFW_05)
    print("up (95%): ", a_975, V_975, R_975, M_975, cNFW_975)
    print("down (95%): ", a_025, V_025, R_025, M_025, cNFW_025)
    print("a = %f +%f -%f (+%f -%f)"
          % (a_md,
             a_84 - a_md,
             a_md - a_16,
             a_95 - a_md,
             a_md - a_05))
    print("V = %f +%f -%f (+%f -%f)"
          % (V_md,
             V_84 - V_md,
             V_md - V_16,
             V_95 - V_md,
             V_md - V_05))
    print("R = %f +%f -%f (+%f -%f)"
          % (R_md,
             R_84 - R_md,
             R_md - R_16,
             R_95 - R_md,
             R_md - R_05))
    print("M = %f +%f -%f (+%f -%f)"
          % (M_md,
             M_84 - M_md,
             M_md - M_16,
             M_95 - M_md,
             M_md - M_05))
    print("c = %f +%f -%f (+%f -%f)"
          % (cNFW_md,
             cNFW_84 - cNFW_md,
             cNFW_md - cNFW_16,
             cNFW_95 - cNFW_md,
             cNFW_md - cNFW_05))

    # Plot single histograms
    plt.close()
    plt.figure()

    ax1 = plt.axes()
    print(len(aPal5), len(prob), len(prob > -1.e4))

    ax1.hist(aPal5[prob > -1.e4], bins=10 ** np.linspace(-1, 1, 20),
             color=mpl.colors.hex2color("#ABCEEA"), histtype="stepfilled")
    #ax1.set_xlim(0.1, 10)
    ax1.set_xscale('log')
    ax1.set_xlabel("a$_{" + streamname + "}$ (pc/Myr$^2$)")
    ax1.set_ylabel("N")
    ax1.axvline(a_md, color='r', linewidth=2)
    ax1.axvline(a_84, color='k', linewidth=2, linestyle='--')
    ax1.axvline(a_16, color='k', linewidth=2, linestyle='--')

    # Save figure
    plt.draw()
    filename = folder + 'chain%s_a%s' % (name, streamname)
    plt.savefig(filename + '.png', format='png', dpi=150, bbox_inches='tight')

    plt.close()
    plt.figure()

    plt.hist(Vsun[prob > -1.e4], 30,
             color=mpl.colors.hex2color("#ABCEEA"), histtype="stepfilled")
    plt.xlabel("V$_{C}$ (km/s)")
    plt.ylabel("N")
    plt.axvline(V_md, color='r', linewidth=2)
    plt.axvline(V_84, color='k', linewidth=2, linestyle='--')
    plt.axvline(V_16, color='k', linewidth=2, linestyle='--')

    # Save figure
    plt.draw()
    filename = folder + 'chain%s_vsun' % name
    plt.savefig(filename + '.png', format='png', dpi=150, bbox_inches='tight')

    plt.close()
    plt.figure()

    plt.hist(cNFW[prob > -1.e4], 30,
             color=mpl.colors.hex2color("#ABCEEA"), histtype="stepfilled")
    plt.xlabel("c")
    plt.ylabel("N")
    plt.axvline(cNFW_md, color='r', linewidth=2)
    plt.axvline(cNFW_84, color='k', linewidth=2, linestyle='--')
    plt.axvline(cNFW_16, color='k', linewidth=2, linestyle='--')

    # Save figure
    plt.draw()
    filename = folder + 'chain%s_c' % name
    plt.savefig(filename + '.png', format='png', dpi=150, bbox_inches='tight')

    mpl.rcParams.update({'font.size': 12})

    # Plot histograms
    plt.close()
    plt.figure()

    fig, axes = plt.subplots(2, 2, figsize=(10.0, 10.0))
    axes[0, 0].hist(aPal5[prob > -1.e4], bins=10 ** np.linspace(-1, 1, 20),
                    color=mpl.colors.hex2color("#ABCEEA"),
                    histtype="stepfilled")
    #axes[0, 0].set_xlim(0.1, 10)
    axes[0, 0].set_xscale('log')
    axes[0, 0].set_xlabel("a$_{" + streamname + "}$ (pc/Myr$^2$)")
    axes[0, 0].set_ylabel("N")
    axes[0, 0].axvline(a_md, color='r', linewidth=1)
    axes[0, 0].axvline(a_m, color='b', linewidth=1)
    axes[0, 0].axvline(a_84, color='k', linewidth=1, linestyle='--')
    axes[0, 0].axvline(a_16, color='k', linewidth=1, linestyle='--')
    axes[0, 0].axvline(a_95, color='k', linewidth=1, linestyle=':')
    axes[0, 0].axvline(a_05, color='k', linewidth=1, linestyle=':')

    axes[1, 0].hist(Vsun[prob > -1.e4], 30,
                    color=mpl.colors.hex2color("#ABCEEA"),
                    histtype="stepfilled")
    axes[1, 0].set_xlabel("V$_{C}$ (km/s)")
    axes[1, 0].set_ylabel("N")
    axes[1, 0].axvline(V_md, color='r', linewidth=1)
    axes[1, 0].axvline(V_m, color='b', linewidth=1)
    axes[1, 0].axvline(V_84, color='k', linewidth=1, linestyle='--')
    axes[1, 0].axvline(V_16, color='k', linewidth=1, linestyle='--')
    axes[1, 0].axvline(V_95, color='k', linewidth=1, linestyle=':')
    axes[1, 0].axvline(V_05, color='k', linewidth=1, linestyle=':')

    axes[1, 1].hist(R200[prob > -1.e4], 30,
                    color=mpl.colors.hex2color("#ABCEEA"),
                    histtype="stepfilled")
    #axes[1, 1].set_xlim(50000, 250000)
    axes[1, 1].set_xlabel("R$_200$")
    axes[1, 1].set_ylabel("N")
    axes[1, 1].axvline(R_md, color='r', linewidth=1)
    axes[1, 1].axvline(R_m, color='b', linewidth=1)
    axes[1, 1].axvline(R_84, color='k', linewidth=1, linestyle='--')
    axes[1, 1].axvline(R_16, color='k', linewidth=1, linestyle='--')
    axes[1, 1].axvline(R_95, color='k', linewidth=1, linestyle=':')
    axes[1, 1].axvline(R_05, color='k', linewidth=1, linestyle=':')

    axes[0, 1].hist(M200[prob > -1.e4], bins=10 ** np.linspace(11, 13, 20),
                    color=mpl.colors.hex2color("#ABCEEA"),
                    histtype="stepfilled")
    #axes[0, 1].set_xlim(1.e11, 1.e13)
    axes[0, 1].set_xscale('log')
    axes[0, 1].set_xlabel("M$_{200}$ (M$_\odot$)")
    axes[0, 1].set_ylabel("N")
    axes[0, 1].set_ylabel("N")
    axes[0, 1].axvline(M_md, color='r', linewidth=1)
    axes[0, 1].axvline(M_m, color='b', linewidth=1)
    axes[0, 1].axvline(M_84, color='k', linewidth=1, linestyle='--')
    axes[0, 1].axvline(M_16, color='k', linewidth=1, linestyle='--')
    axes[0, 1].axvline(M_95, color='k', linewidth=1, linestyle=':')
    axes[0, 1].axvline(M_05, color='k', linewidth=1, linestyle=':')

    # Save figure
    plt.draw()
    filename = folder + 'chain%s_params_histo.acc' % name
    plt.savefig(filename + '.png', format='png', dpi=150, bbox_inches='tight')


def param_histo_stream(name, streamname='stream', folder='./'):
    mpl.rcParams.update({'font.size': 12})

    print "File name: chain%s.dat" % name

        # Plot histograms of stream parameters
    prob,\
        mass_cluster,\
        pm_mu_delta,\
        pm_mu_alphacosdelta,\
        mass_halo,\
        distance_cluster,\
        mass_loss_rate,\
        q_halo,\
        r_halo,\
        tpast,\
        sigma_x,\
        sigma_v,\
        sigma_vx,\
        sigma_mu,\
        rgalsun,\
        vLSR = np.loadtxt("chain%s.dat" % name, unpack=True)
    thin = 1
    prob = prob[::thin]
    mass_cluster = mass_cluster[::thin]
    pm_mu_delta = pm_mu_delta[::thin]
    pm_mu_alphacosdelta = pm_mu_alphacosdelta[::thin]
    distance_cluster = distance_cluster[::thin]

    #step = np.arange(len(prob[prob > -1.e4]))

    maxid = np.argmax(prob)
    Mcl_b = 0.001 * mass_cluster[maxid]
    mud_b = 1000 * pm_mu_delta[maxid]
    mua_b = 1000 * pm_mu_alphacosdelta[maxid]
    d_b = distance_cluster[maxid]
    Mcl_m = 0.001 * np.mean(mass_cluster)
    mud_m = 1000 * np.mean(pm_mu_delta)
    mua_m = 1000 * np.mean(pm_mu_alphacosdelta)
    d_m = np.mean(distance_cluster)
    Mcl_md = 0.001 * np.median(mass_cluster)
    mud_md = 1000 * np.median(pm_mu_delta)
    mua_md = 1000 * np.median(pm_mu_alphacosdelta)
    d_md = np.median(distance_cluster)
    Mcl_84 = 0.001 * np.percentile(mass_cluster, 84)
    mud_84 = 1000 * np.percentile(pm_mu_delta, 84)
    mua_84 = 1000 * np.percentile(pm_mu_alphacosdelta, 84)
    d_84 = np.percentile(distance_cluster, 84)
    Mcl_16 = 0.001 * np.percentile(mass_cluster, 16)
    mud_16 = 1000 * np.percentile(pm_mu_delta, 16)
    mua_16 = 1000 * np.percentile(pm_mu_alphacosdelta, 16)
    d_16 = np.percentile(distance_cluster, 16)
    Mcl_95 = 0.001 * np.percentile(mass_cluster, 95)
    mud_95 = 1000 * np.percentile(pm_mu_delta, 95)
    mua_95 = 1000 * np.percentile(pm_mu_alphacosdelta, 95)
    d_95 = np.percentile(distance_cluster, 95)
    Mcl_05 = 0.001 * np.percentile(mass_cluster, 5)
    mud_05 = 1000 * np.percentile(pm_mu_delta, 5)
    mua_05 = 1000 * np.percentile(pm_mu_alphacosdelta, 5)
    d_05 = np.percentile(distance_cluster, 5)
    Mcl_975 = 0.001 * np.percentile(mass_cluster, 97.5)
    mud_975 = 1000 * np.percentile(pm_mu_delta, 97.5)
    mua_975 = 1000 * np.percentile(pm_mu_alphacosdelta, 97.5)
    d_975 = np.percentile(distance_cluster, 97.5)
    Mcl_025 = 0.001 * np.percentile(mass_cluster, 2.5)
    mud_025 = 1000 * np.percentile(pm_mu_delta, 2.5)
    mua_025 = 1000 * np.percentile(pm_mu_alphacosdelta, 2.5)
    d_025 = np.percentile(distance_cluster, 2.5)
    print("Best-fit: ", Mcl_b, mud_b, mua_b, d_b)
    print("Mean: ", Mcl_m, mud_m, mua_m, d_m)
    print("Median: ", Mcl_md, mud_md, mua_md, d_md)
    print("up (68%): ", Mcl_84, mud_84, mua_84, d_84)
    print("down (68%): ", Mcl_16, mud_16, mua_16, d_16)
    print("up (90%): ", Mcl_95, mud_95, mua_95, d_95)
    print("down (90%): ", Mcl_05, mud_05, mua_05, d_05)
    print("up (95%): ", Mcl_975, mud_975, mua_975, d_975)
    print("down (95%): ", Mcl_025, mud_025, mua_025, d_025)
    print("M" + streamname + " = %f +%f -%f (+%f -%f)"
          % (Mcl_md,
             Mcl_84 - Mcl_md,
             Mcl_md - Mcl_16,
             Mcl_95 - Mcl_md,
             Mcl_md - Mcl_05))
    print("mu_delta = %f +%f -%f (+%f -%f)"
          % (mud_md,
             mud_84 - mud_md,
             mud_md - mud_16,
             mud_95 - mud_md,
             mud_md - mud_05))
    print("mu_alpha*cos(delta) = %f +%f -%f (+%f -%f)"
          % (mua_md,
             mua_84 - mua_md,
             mua_md - mua_16,
             mua_95 - mua_md,
             mua_md - mua_05))
    print("d = %f +%f -%f (+%f -%f)"
          % (d_md,
             d_84 - d_md,
             d_md - d_16,
             d_95 - d_md,
             d_md - d_05))

    # Plot histograms
    plt.close()
    plt.figure()

    fig, axes = plt.subplots(2, 2, figsize=(10.0, 10.0))
    axes[0, 0].hist(mass_cluster[prob > -1.e4] / 1000.0, 30,
                    color=mpl.colors.hex2color("#ABCEEA"),
                    histtype="stepfilled")
    #axes[0, 0].set_xlim(0, 60)
    axes[0, 0].set_xlabel("M$_{cl}$ ($10^3$ M$_\odot$)")
    axes[0, 0].set_ylabel("N")
    axes[0, 0].axvline(Mcl_md, color='r', linewidth=1)
    axes[0, 0].axvline(Mcl_m, color='b', linewidth=1)
    axes[0, 0].axvline(Mcl_84, color='k', linewidth=1, linestyle='--')
    axes[0, 0].axvline(Mcl_16, color='k', linewidth=1, linestyle='--')
    axes[0, 0].axvline(Mcl_95, color='k', linewidth=1, linestyle=':')
    axes[0, 0].axvline(Mcl_05, color='k', linewidth=1, linestyle=':')

    axes[1, 0].hist(1000.0 * pm_mu_delta[prob > -1.e4], 30,
                    color=mpl.colors.hex2color("#ABCEEA"),
                    histtype="stepfilled")
    #axes[1, 0].set_xlim(-3, -1.5)
    axes[1, 0].set_xlabel("$\mu_\delta$ (mas/yr)")
    axes[1, 0].set_ylabel("N")
    axes[1, 0].axvline(mud_md, color='r', linewidth=1)
    axes[1, 0].axvline(mud_m, color='b', linewidth=1)
    axes[1, 0].axvline(mud_84, color='k', linewidth=1, linestyle='--')
    axes[1, 0].axvline(mud_16, color='k', linewidth=1, linestyle='--')
    axes[1, 0].axvline(mud_95, color='k', linewidth=1, linestyle=':')
    axes[1, 0].axvline(mud_05, color='k', linewidth=1, linestyle=':')

    axes[1, 1].hist(1000.0 * pm_mu_alphacosdelta[prob > -1.e4], 30,
                    color=mpl.colors.hex2color("#ABCEEA"),
                    histtype="stepfilled")
    #axes[1, 1].set_xlim(-3, -1.5)
    axes[1, 1].set_xlabel("$\mu_{\\alpha} \cos{\delta}$ (mas/yr)")
    axes[1, 1].set_ylabel("N")
    axes[1, 1].axvline(mua_md, color='r', linewidth=1)
    axes[1, 1].axvline(mua_m, color='b', linewidth=1)
    axes[1, 1].axvline(mua_84, color='k', linewidth=1, linestyle='--')
    axes[1, 1].axvline(mua_16, color='k', linewidth=1, linestyle='--')
    axes[1, 1].axvline(mua_95, color='k', linewidth=1, linestyle=':')
    axes[1, 1].axvline(mua_05, color='k', linewidth=1, linestyle=':')

    axes[0, 1].hist(distance_cluster[prob > -1.e4],
                    30,
                    color=mpl.colors.hex2color("#ABCEEA"),
                    histtype="stepfilled")
    #axes[0, 1].set_xlim(20, 27)
    axes[0, 1].set_xlabel("d (kpc)")
    axes[0, 1].set_ylabel("N")
    axes[0, 1].axvline(d_md, color='r', linewidth=1)
    axes[0, 1].axvline(d_m, color='b', linewidth=1)
    axes[0, 1].axvline(d_84, color='k', linewidth=1, linestyle='--')
    axes[0, 1].axvline(d_16, color='k', linewidth=1, linestyle='--')
    axes[0, 1].axvline(d_95, color='k', linewidth=1, linestyle=':')
    axes[0, 1].axvline(d_05, color='k', linewidth=1, linestyle=':')

    # Save figure
    plt.draw()
    filename = folder + 'chain%s_params_histo.%s' % (name, streamname)
    plt.savefig(filename + '.png', format='png', dpi=150, bbox_inches='tight')


def plot_convergence(name, streamname='stream', folder='./'):
    mpl.rcParams.update({'font.size': 20})

    print "File name: chain%s.dat" % name

    prob,\
        mass_cluster,\
        pm_mu_delta,\
        pm_mu_alphacosdelta,\
        mass_halo,\
        distance_cluster,\
        mass_loss_rate,\
        q_halo,\
        r_halo,\
        tpast,\
        sigma_x,\
        sigma_v,\
        sigma_vx,\
        sigma_mu,\
        rgalsun,\
        vLSR = np.loadtxt("chain%s.dat" % name, unpack=True)

    thin = 1
    prob = prob[::thin]
    prob_m = np.mean(prob)
    prob_md = np.median(prob)
    prob_84 = np.percentile(prob, 84)
    prob_16 = np.percentile(prob, 16)

    step = np.arange(len(prob[prob > -1.e4]))

    plt.close()
    plt.figure()

    plt.plot(step, prob[prob > -1.e4], 'k-')
    # plt.xlim(0,300000)
    plt.xlabel("MCMC step")
    plt.ylabel("$log(L)$")
    plt.axhline(prob_md, color='r', linewidth=1)
    plt.axhline(prob_m, color='b', linewidth=1)
    plt.axhline(prob_84, color='k', linewidth=1, linestyle='--')
    plt.axhline(prob_16, color='k', linewidth=1, linestyle='--')

    # Save figure
    plt.draw()
    filename = folder + 'chain%s_P.%s' % (name, streamname)
    plt.savefig(filename + '.png', format='png', dpi=150, bbox_inches='tight')


def plot_model(stream_data,
               VRdata,
               ODdata,
               modelname,
               name='',
               streamname='',
               folder='./',
               radec=False,
               clusterpos=[211.36370833, 28.53444444],
               errors=False):

    ltemp = stream_data[:, 0]
    btemp = stream_data[:, 1]
    vrtemp = stream_data[:, 3]
    RAtemp = stream_data[:, 4]
    DECtemp = stream_data[:, 5]

    lcor = ltemp

    for i in xrange(len(ltemp)):
        if ltemp[i] < 180.0:
            lcor[i] = ltemp[i] * cos(btemp[i] * np.pi / 180.0)
        else:
            lcor[i] = (ltemp[i] - 360.0) * cos(btemp[i] * np.pi / 180.0)

    ODx = range(len(ODdata))
    ODy = range(len(ODdata))

    if radec:
        ODx = ODdata[:, 2]  # RA
        ODy = ODdata[:, 3]  # Dec

    else:
        for i in xrange(len(ODdata)):
            if ODdata[i][0] < 180.0:
                ODx[i] = ODdata[i][0] * cos(ODdata[i][1] * np.pi / 180.0)
            else:
                ODx[i] = (ODdata[i][0] - 360.0) * cos(ODdata[i][1]
                                                      * np.pi / 180.0)
            ODy[i] = ODdata[i][1]

    VRx = range(len(VRdata))
    VRy = range(len(VRdata))

    if radec:
        VRx = VRdata[:, 5]
    else:
        VRx = VRdata[:, 4] * np.cos(VRdata[:, 3] * np.pi / 180.0)

    VRy = VRdata[:, 0]
    vrerr = VRdata[:, 1]

    plt.close()
    fig, ax = plt.subplots(2, 2, figsize=(10.0, 5))
    fig.subplots_adjust(hspace=0.5)

    xmin = np.min(VRx)
    xmax = np.max(VRx)
    ymin = np.min(VRy)
    ymax = np.max(VRy)

    xlabel = ("l cos(b) [deg]", "RA [deg]")[radec]

    xpos = (lcor, RAtemp)[radec]

    ax[0, 0].scatter(xpos, vrtemp, alpha=.1, edgecolors='none')
    ax[0, 0].scatter(VRx, VRy, facecolor='none', edgecolors='red', alpha=0.5)
    if errors:
        ax[0, 0].errorbar(VRx, VRy, 0.5, vrerr, fmt=None,
                          ecolor='red', alpha=0.5)
    ax[0, 0].set_xlim(220, 205)
    ax[0, 0].set_ylim(100, 120)
    ax[0, 0].set_xlabel(xlabel)
    ax[0, 0].set_ylabel("V$_R$ [km/s]")

    ylabel = ("b [deg]", "DEC [deg]")[radec]

    xmin = (10, 225)[radec]
    xmax = (20, 190)[radec]
    ymin = (50, 20)[radec]
    ymax = (90, 38)[radec]
    ypos = (btemp, DECtemp)[radec]

    ax[0, 1].scatter(xpos, ypos, alpha=.1, edgecolors='none')
    ax[0, 1].scatter(ODx, ODy, facecolor='none', edgecolors='green', alpha=0.5)
    if errors:
        ax[0, 1].errorbar(ODx, ODy, 1.0, 1.0, fmt=None,
                          ecolor='green', alpha=0.5)
    ax[0, 1].scatter(ODx[0:3], ODy[0:3],
                     facecolor='none',
                     edgecolors='black',
                     alpha=0.5)
    ax[0, 1].plot(clusterpos[0], clusterpos[1], 'r^')
    ax[0, 1].set_xlim(xmin, xmax)
    ax[0, 1].set_ylim(ymin, ymax)
    ax[0, 1].set_xlabel(xlabel)
    ax[0, 1].set_ylabel(ylabel)

    VRy = (VRdata[:, 3], VRdata[:, 6])[radec]
    VRx = VRdata[:, 0]

    xmin = np.min(VRx)
    xmax = np.max(VRx)
    ymin = np.min(VRy)
    ymax = np.max(VRy)

    ax[1, 1].scatter(vrtemp, ypos, alpha=.1, edgecolors='none')
    ax[1, 1].scatter(VRx, VRy, facecolor='none', edgecolors='red', alpha=0.5)
    if errors:
        ax[1, 1].errorbar(VRx, VRy, vrerr, 0.5, fmt=None,
                          ecolor='red', alpha=0.5)
    ax[1, 1].set_xlim(100, 120)
    ax[1, 1].set_ylim(24, 32)
    ax[1, 1].set_ylabel(ylabel)
    ax[1, 1].set_xlabel("V$_R$ [km/s]")

    VRx = (VRdata[:, 4] * np.cos(VRdata[:, 3] * np.pi / 180.0),
           VRdata[:, 5])[radec]
    VRy = (VRdata[:, 3], VRdata[:, 6])[radec]

    xmin = np.min(VRx) * 0.1
    xmax = np.max(VRx) * 10.0
    ymin = np.min(VRy) * 0.1
    ymax = np.max(VRy) * 10.0

    ax[1, 0].scatter(xpos, ypos, alpha=.1, edgecolors='none')
    ax[1, 0].scatter(ODx, ODy, facecolor='none', edgecolors='green', alpha=0.5)
    if errors:
        ax[1, 0].errorbar(ODx, ODy, 1.0, 1.0, fmt=None,
                          ecolor='green', alpha=0.5)
    ax[1, 0].scatter(VRx, VRy, facecolor='none', edgecolors='red', alpha=0.5)
    if errors:
        ax[1, 0].errorbar(VRx, VRy, 0.5, 0.5, fmt=None,
                          ecolor='red', alpha=0.5)
    ax[1, 0].set_xlim(220, 205)
    ax[1, 0].set_ylim(24, 32)
    ax[1, 0].set_xlabel(xlabel)
    ax[1, 0].set_ylabel(ylabel)

    # Save or show figure
    #plt.draw()
    if modelname == 'Last':
        plt.draw()
    elif modelname == 'Temp':
        filename = folder + 'chain%s_model.%s' % (name, streamname)
        plt.savefig(filename + '.png',
                    format='png',
                    dpi=150,
                    bbox_inches='tight')
    else:
        filename = folder + 'chain%s_model%s.%s'\
            % (name, modelname, streamname)
        plt.savefig(filename + '.png',
                    format='png',
                    dpi=150,
                    bbox_inches='tight')


def plot_mstream(name, streamname='stream', folder='./'):
    mpl.rcParams.update({'font.size': 20})

    print "File name: chain%s.dat" % name

        # Plot histograms of Pal 5 parameters
    prob,\
        mass_cluster,\
        pm_mu_delta,\
        pm_mu_alphacosdelta,\
        mass_halo,\
        distance_cluster,\
        mass_loss_rate,\
        q_halo,\
        r_halo,\
        tpast,\
        sigma_x,\
        sigma_v,\
        sigma_vx,\
        sigma_mu,\
        rgalsun,\
        vLSR = np.loadtxt("chain%s.dat" % name, unpack=True)

    Mcl_md = 0.001 * np.median(mass_cluster)
    Mcl_84 = 0.001 * np.percentile(mass_cluster, 84)
    Mcl_16 = 0.001 * np.percentile(mass_cluster, 16)
    Mcl_95 = 0.001 * np.percentile(mass_cluster, 95)
    Mcl_05 = 0.001 * np.percentile(mass_cluster, 5)

    print("MPal5 = %f +%f -%f (+%f -%f)"
          % (Mcl_md,
             Mcl_84 - Mcl_md,
             Mcl_md - Mcl_16,
             Mcl_95 - Mcl_md,
             Mcl_md - Mcl_05))

    # Plot histograms
    plt.close()
    plt.figure()

    plt.hist(mass_cluster[prob > -1.e4] / 1000.0, 30,
             color=mpl.colors.hex2color("#ABCEEA"), histtype="stepfilled")
    #plt.xlim(0, 60)
    plt.xlabel("M$_{" + streamname + "}$ ($10^3$M$_\odot$)")
    plt.ylabel("N")
    plt.axvline(Mcl_md, color='r', linewidth=2)
    plt.axvline(Mcl_84, color='k', linewidth=2, linestyle='--')
    plt.axvline(Mcl_16, color='k', linewidth=2, linestyle='--')

    # Save figure
    plt.draw()
    filename = folder + 'chain%s_mpal5' % name
    plt.savefig(filename + '.png', format='png', dpi=150, bbox_inches='tight')


def plot_mudelta(name, folder='./'):
    mpl.rcParams.update({'font.size': 20})

    print "File name: chain%s.dat" % name

        # Plot histograms of Pal 5 parameters
    prob,\
        mass_cluster,\
        pm_mu_delta,\
        pm_mu_alphacosdelta,\
        mass_halo,\
        distance_cluster,\
        mass_loss_rate,\
        q_halo,\
        r_halo,\
        tpast,\
        sigma_x,\
        sigma_v,\
        sigma_vx,\
        sigma_mu,\
        rgalsun,\
        vLSR = np.loadtxt("chain%s.dat" % name, unpack=True)

    mud_md = 1000 * np.median(pm_mu_delta)
    mud_84 = 1000 * np.percentile(pm_mu_delta, 84)
    mud_16 = 1000 * np.percentile(pm_mu_delta, 16)
    mud_95 = 1000 * np.percentile(pm_mu_delta, 95)
    mud_05 = 1000 * np.percentile(pm_mu_delta, 5)

    print("mu_delta = %f +%f -%f (+%f -%f)"
          % (mud_md,
             mud_84 - mud_md,
             mud_md - mud_16,
             mud_95 - mud_md,
             mud_md - mud_05))

    # Plot histograms
    plt.close()
    plt.figure()

    plt.hist(1000.0 * pm_mu_delta[prob > -1.e4], 30,
             color=mpl.colors.hex2color("#ABCEEA"), histtype="stepfilled")
    #plt.xlim(-3, -1.5)
    plt.xlabel("$\mu_\delta$ (mas/yr)")
    plt.ylabel("N")
    plt.axvline(mud_md, color='r', linewidth=2)
    plt.axvline(mud_84, color='k', linewidth=2, linestyle='--')
    plt.axvline(mud_16, color='k', linewidth=2, linestyle='--')

    # Save figure
    plt.draw()
    filename = folder + 'chain%s_mudelta' % name
    plt.savefig(filename + '.png', format='png', dpi=150, bbox_inches='tight')


def plot_mualpha(name, folder='./'):
    mpl.rcParams.update({'font.size': 20})

    print "File name: chain%s.dat" % name

        # Plot histograms of Pal 5 parameters
    prob,\
        mass_cluster,\
        pm_mu_delta,\
        pm_mu_alphacosdelta,\
        mass_halo,\
        distance_cluster,\
        mass_loss_rate,\
        q_halo,\
        r_halo,\
        tpast,\
        sigma_x,\
        sigma_v,\
        sigma_vx,\
        sigma_mu,\
        rgalsun,\
        vLSR = np.loadtxt("chain%s.dat" % name, unpack=True)

    mua_md = 1000 * np.median(pm_mu_alphacosdelta)
    mua_84 = 1000 * np.percentile(pm_mu_alphacosdelta, 84)
    mua_16 = 1000 * np.percentile(pm_mu_alphacosdelta, 16)
    mua_95 = 1000 * np.percentile(pm_mu_alphacosdelta, 95)
    mua_05 = 1000 * np.percentile(pm_mu_alphacosdelta, 5)

    print("mu_alpha*cos(delta) = %f +%f -%f (+%f -%f)"
          % (mua_md,
             mua_84 - mua_md,
             mua_md - mua_16,
             mua_95 - mua_md,
             mua_md - mua_05))

    # Plot histograms
    plt.close()
    plt.figure()

    plt.hist(1000.0 * pm_mu_alphacosdelta[prob > -1.e4], 30,
             color=mpl.colors.hex2color("#ABCEEA"), histtype="stepfilled")
    #plt.xlim(-3, -1.5)
    plt.xlabel("$\mu_{\\alpha} \cos{\delta}$ (mas/yr)")
    plt.ylabel("N")
    plt.axvline(mua_md, color='r', linewidth=2)
    plt.axvline(mua_84, color='k', linewidth=2, linestyle='--')
    plt.axvline(mua_16, color='k', linewidth=2, linestyle='--')

    # Save figure
    plt.draw()
    filename = folder + 'chain%s_mualpha' % name
    plt.savefig(filename + '.png', format='png', dpi=150, bbox_inches='tight')


def plot_distance(name, folder='./'):
    mpl.rcParams.update({'font.size': 20})

    print "File name: chain%s.dat" % name

        # Plot histograms of Pal 5 parameters
    prob,\
        mass_cluster,\
        pm_mu_delta,\
        pm_mu_alphacosdelta,\
        mass_halo,\
        distance_cluster,\
        mass_loss_rate,\
        q_halo,\
        r_halo,\
        tpast,\
        sigma_x,\
        sigma_v,\
        sigma_vx,\
        sigma_mu,\
        rgalsun,\
        vLSR = np.loadtxt("chain%s.dat" % name, unpack=True)

    d_md = np.median(distance_cluster)
    d_84 = np.percentile(distance_cluster, 84)
    d_16 = np.percentile(distance_cluster, 16)
    d_95 = np.percentile(distance_cluster, 95)
    d_05 = np.percentile(distance_cluster, 5)

    print("d = %f +%f -%f (+%f -%f)"
          % (d_md,
             d_84 - d_md,
             d_md - d_16,
             d_95 - d_md,
             d_md - d_05))

    # Plot histograms
    plt.close()
    plt.figure()

    plt.hist(distance_cluster[prob > -1.e4],
             30,
             color=mpl.colors.hex2color("#ABCEEA"),
             histtype="stepfilled")
    #plt.xlim(20, 27)
    plt.xlabel("d (kpc)")
    plt.ylabel("N")
    plt.axvline(d_md, color='r', linewidth=2)
    plt.axvline(d_84, color='k', linewidth=2, linestyle='--')
    plt.axvline(d_16, color='k', linewidth=2, linestyle='--')

    # Save figure
    plt.draw()
    filename = folder + 'chain%s_distance' % name
    plt.savefig(filename + '.png', format='png', dpi=150, bbox_inches='tight')


def plot_mhalo(name, folder='./'):
    mpl.rcParams.update({'font.size': 20})

    print "File name: chain%s.dat" % name

        # Plot histograms of halo parameters
    prob,\
        mass_cluster,\
        pm_mu_delta,\
        pm_mu_alphacosdelta,\
        mass_halo,\
        distance_cluster,\
        mass_loss_rate,\
        q_halo,\
        r_halo,\
        tpast,\
        sigma_x,\
        sigma_v,\
        sigma_vx,\
        sigma_mu,\
        rgalsun,\
        vLSR = np.loadtxt("chain%s.dat" % name, unpack=True)

    Mh_md = np.median(mass_halo)
    Mh_84 = np.percentile(mass_halo, 84)
    Mh_16 = np.percentile(mass_halo, 16)
    Mh_95 = np.percentile(mass_halo, 95)
    Mh_05 = np.percentile(mass_halo, 5)

    print("Mhalo = %e +%e -%e (+%e -%e)"
          % (Mh_md,
             Mh_84 - Mh_md,
             Mh_md - Mh_16,
             Mh_95 - Mh_md,
             Mh_md - Mh_05))

    # Plot histograms
    plt.close()
    plt.figure()

    ax1 = plt.axes()
    plt.hist(mass_halo[prob > -1.e4], bins=10 ** np.linspace(11, 13, 20),
             color=mpl.colors.hex2color("#ABCEEA"), histtype="stepfilled")
    #plt.xlim(1.e11, 1.e13)
    ax1.set_xscale('log')
    plt.xlabel("M$_{halo}$ (M$_\odot$)")
    plt.ylabel("N")
    plt.axvline(Mh_md, color='r', linewidth=2)
    plt.axvline(Mh_84, color='k', linewidth=2, linestyle='--')
    plt.axvline(Mh_16, color='k', linewidth=2, linestyle='--')
    # Save figure
    plt.draw()
    filename = folder + 'chain%s_mhalo' % name
    plt.savefig(filename + '.png', format='png', dpi=150, bbox_inches='tight')


def plot_rhalo(name, folder='./'):
    mpl.rcParams.update({'font.size': 20})

    print "File name: chain%s.dat" % name

        # Plot histograms of halo parameters
    prob,\
        mass_cluster,\
        pm_mu_delta,\
        pm_mu_alphacosdelta,\
        mass_halo,\
        distance_cluster,\
        mass_loss_rate,\
        q_halo,\
        r_halo,\
        tpast,\
        sigma_x,\
        sigma_v,\
        sigma_vx,\
        sigma_mu,\
        rgalsun,\
        vLSR = np.loadtxt("chain%s.dat" % name, unpack=True)

    rh_md = np.median(r_halo)
    rh_84 = np.percentile(r_halo, 84)
    rh_16 = np.percentile(r_halo, 16)
    rh_95 = np.percentile(r_halo, 95)
    rh_05 = np.percentile(r_halo, 5)

    print("Rhalo = %f +%f -%f (+%f -%f)"
          % (rh_md,
             rh_84 - rh_md,
             rh_md - rh_16,
             rh_95 - rh_md,
             rh_md - rh_05))

    # Plot histograms
    plt.close()
    plt.figure()

    plt.hist(r_halo[prob > -1.e4], 30,
             color=mpl.colors.hex2color("#ABCEEA"), histtype="stepfilled")
    plt.xlim(0, 50000)
    plt.xlabel("R$_{halo}$ (pc)")
    plt.ylabel("N")
    plt.axvline(rh_md, color='r', linewidth=2)
    plt.axvline(rh_84, color='k', linewidth=2, linestyle='--')
    plt.axvline(rh_16, color='k', linewidth=2, linestyle='--')
    # Save figure
    plt.draw()
    filename = folder + 'chain%s_rhalo' % name
    plt.savefig(filename + '.png', format='png', dpi=150, bbox_inches='tight')


def plot_qhalo(name, folder='./'):
    mpl.rcParams.update({'font.size': 20})

    print "File name: chain%s.dat" % name

        # Plot histograms of halo parameters
    prob,\
        mass_cluster,\
        pm_mu_delta,\
        pm_mu_alphacosdelta,\
        mass_halo,\
        distance_cluster,\
        mass_loss_rate,\
        q_halo,\
        r_halo,\
        tpast,\
        sigma_x,\
        sigma_v,\
        sigma_vx,\
        sigma_mu,\
        rgalsun,\
        vLSR = np.loadtxt("chain%s.dat" % name, unpack=True)

    q_md = np.median(q_halo)
    q_84 = np.percentile(q_halo, 84)
    q_16 = np.percentile(q_halo, 16)
    q_95 = np.percentile(q_halo, 95)
    q_05 = np.percentile(q_halo, 5)

    print("qz = %f +%f -%f (+%f -%f)"
          % (q_md,
             q_84 - q_md,
             q_md - q_16,
             q_95 - q_md,
             q_md - q_05))

    # Plot histograms
    plt.close()
    plt.figure()

    plt.hist(q_halo[prob > -1.e4], 30,
             color=mpl.colors.hex2color("#ABCEEA"), histtype="stepfilled")
    plt.xlim(0.5, 2)
    plt.xlabel("q$_z$")
    plt.ylabel("N")
    plt.axvline(q_md, color='r', linewidth=2)
    plt.axvline(q_84, color='k', linewidth=2, linestyle='--')
    plt.axvline(q_16, color='k', linewidth=2, linestyle='--')

    # Save figure
    plt.draw()
    filename = folder + 'chain%s_qhalo' % name
    plt.savefig(filename + '.png', format='png', dpi=150, bbox_inches='tight')


def plot_dmdt(name, folder='./'):
    mpl.rcParams.update({'font.size': 20})

    print "File name: chain%s.dat" % name

        # Plot histograms of halo parameters
    prob,\
        mass_cluster,\
        pm_mu_delta,\
        pm_mu_alphacosdelta,\
        mass_halo,\
        distance_cluster,\
        mass_loss_rate,\
        q_halo,\
        r_halo,\
        tpast,\
        sigma_x,\
        sigma_v,\
        sigma_vx,\
        sigma_mu,\
        rgalsun,\
        vLSR = np.loadtxt("chain%s.dat" % name, unpack=True)

    dM_md = np.median(mass_loss_rate)
    dM_84 = np.percentile(mass_loss_rate, 84)
    dM_16 = np.percentile(mass_loss_rate, 16)
    dM_95 = np.percentile(mass_loss_rate, 95)
    dM_05 = np.percentile(mass_loss_rate, 5)
    print("dM = %f +%f -%f (+%f -%f)"
          % (dM_md,
             dM_84 - dM_md,
             dM_md - dM_16,
             dM_95 - dM_md,
             dM_md - dM_05))

    # Plot histograms
    plt.close()
    plt.figure()

    plt.hist(mass_loss_rate[prob > -1.e4],
             30,
             color=mpl.colors.hex2color("#ABCEEA"),
             histtype="stepfilled")
    plt.xlim(0, 15)
    plt.xlabel("dM/dt (M$_\odot/Myr$)")
    plt.ylabel("N")
    plt.axvline(dM_md, color='r', linewidth=2)
    plt.axvline(dM_84, color='k', linewidth=2, linestyle='--')
    plt.axvline(dM_16, color='k', linewidth=2, linestyle='--')

    # Save figure
    plt.draw()
    filename = folder + 'chain%s_dmdt' % name
    plt.savefig(filename + '.png', format='png', dpi=150, bbox_inches='tight')


def plot_vLSR(name, folder='./'):
    mpl.rcParams.update({'font.size': 20})

    print "File name: chain%s.dat" % name

        # Plot histograms of halo parameters
    prob,\
        mass_cluster,\
        pm_mu_delta,\
        pm_mu_alphacosdelta,\
        mass_halo,\
        distance_cluster,\
        mass_loss_rate,\
        q_halo,\
        r_halo,\
        tpast,\
        sigma_x,\
        sigma_v,\
        sigma_vx,\
        sigma_mu,\
        rgalsun,\
        vLSR = np.loadtxt("chain%s.dat" % name, unpack=True)

    vLSR_md = np.median(vLSR)
    vLSR_84 = np.percentile(vLSR, 84)
    vLSR_16 = np.percentile(vLSR, 16)
    vLSR_95 = np.percentile(vLSR, 95)
    vLSR_05 = np.percentile(vLSR, 5)

    print("vLSR = %f +%f -%f (+%f -%f)"
          % (vLSR_md,
             vLSR_84 - vLSR_md,
             vLSR_md - vLSR_16,
             vLSR_95 - vLSR_md,
             vLSR_md - vLSR_05))

    # Plot histograms
    plt.close()
    plt.figure()

    plt.hist(vLSR[prob > -1.e4], 30,
             color=mpl.colors.hex2color("#ABCEEA"), histtype="stepfilled")
    plt.xlim(200, 280)
    plt.xlabel("V$_{LSR}$ (km/s)")
    plt.ylabel("N")
    plt.axvline(vLSR_md, color='r', linewidth=2)
    plt.axvline(vLSR_84, color='k', linewidth=2, linestyle='--')
    plt.axvline(vLSR_16, color='k', linewidth=2, linestyle='--')

    # Save figure
    plt.draw()
    filename = folder + 'chain%s_vLSR' % name
    plt.savefig(filename + '.png', format='png', dpi=150, bbox_inches='tight')


def plot_rgalsun(name, folder='./'):
    mpl.rcParams.update({'font.size': 20})

    print "File name: chain%s.dat" % name

        # Plot histograms of halo parameters
    prob,\
        mass_cluster,\
        pm_mu_delta,\
        pm_mu_alphacosdelta,\
        mass_halo,\
        distance_cluster,\
        mass_loss_rate,\
        q_halo,\
        r_halo,\
        tpast,\
        sigma_x,\
        sigma_v,\
        sigma_vx,\
        sigma_mu,\
        rgalsun,\
        vLSR = np.loadtxt("chain%s.dat" % name, unpack=True)

    rgalsun_md = np.median(rgalsun)
    rgalsun_84 = np.percentile(rgalsun, 84)
    rgalsun_16 = np.percentile(rgalsun, 16)
    rgalsun_95 = np.percentile(rgalsun, 95)
    rgalsun_05 = np.percentile(rgalsun, 5)

    print("rgalsun = %f +%f -%f (+%f -%f)"
          % (rgalsun_md,
             rgalsun_84 - rgalsun_md,
             rgalsun_md - rgalsun_16,
             rgalsun_95 - rgalsun_md,
             rgalsun_md - rgalsun_05))

    # Plot histograms
    plt.close()
    plt.figure()

    plt.hist(rgalsun[prob > -1.e4], 30,
             color=mpl.colors.hex2color("#ABCEEA"), histtype="stepfilled")
    plt.xlim(7500, 9000)
    plt.xlabel("R$_{\odot}$ (km/s)")
    plt.ylabel("N")
    plt.axvline(rgalsun_md, color='r', linewidth=2)
    plt.axvline(rgalsun_84, color='k', linewidth=2, linestyle='--')
    plt.axvline(rgalsun_16, color='k', linewidth=2, linestyle='--')

    # Save figure
    plt.draw()
    filename = folder + 'chain%s_rgalsun' % name
    plt.savefig(filename + '.png', format='png', dpi=150, bbox_inches='tight')


def plot_triangle(name,
                  streamname='stream',
                  gal_latitude=42.15,
                  gal_longitude=73.59,
                  folder='./'):
    # Setup
    txt = 11
    mpl.rcParams['axes.labelsize'] = txt
    mpl.rcParams['xtick.labelsize'] = txt
    mpl.rcParams['ytick.labelsize'] = txt

    # Choice of stream type
    prob,\
        mass_cluster,\
        pm_mu_delta,\
        pm_mu_alphacosdelta,\
        mass_halo,\
        distance_cluster,\
        mass_loss_rate,\
        q_halo,\
        r_halo,\
        tpast,\
        sigma_x,\
        sigma_v,\
        sigma_vx,\
        sigma_mu,\
        rgalsun,\
        vLSR = np.loadtxt("chain%s.dat" % name, unpack=True)

        # Plot histograms of acceleration parameters

    good_inds = [good_vals(vals) for vals in [prob,
                                              rgalsun,
                                              mass_halo,
                                              q_halo,
                                              r_halo,
                                              distance_cluster]]
    ind = good_inds[0]
    for i in xrange(1, len(good_inds)):
        ind = np.intersect1d(good_inds[i], ind)

    prob = prob[ind]
    rgalsun = rgalsun[ind]
    mass_halo = mass_halo[ind]
    q_halo = q_halo[ind]
    r_halo = r_halo[ind]
    distance_cluster = distance_cluster[ind]

    distance_cluster = distance_cluster * 1000.0
    Vsun,\
        R200,\
        cNFW,\
        M200,\
        aPal5 = calculate_halo_properties(prob,
                                          rgalsun,
                                          mass_halo,
                                          q_halo,
                                          r_halo,
                                          distance_cluster,
                                          gal_latitude,
                                          gal_longitude)

    samples = np.array([mass_cluster,
                        pm_mu_delta,
                        pm_mu_alphacosdelta,
                        distance_cluster,
                        mass_loss_rate,
                        mass_halo,
                        r_halo,
                        q_halo,
                        tpast,
                        rgalsun,
                        vLSR,
                        sigma_x,
                        sigma_v,
                        sigma_vx,
                        aPal5,
                        Vsun,
                        M200,
                        R200,
                        cNFW]).T

    plt.close()
    triangle.corner(samples,
                    labels=["M$_{" + streamname + "}$ [M$_\odot$]",
                            "$\mu_\delta$ [arcsec/yr]",
                            "$\mu_\\alpha\cos(\delta)$ [arcsec/yr]",
                            "d$_{" + streamname + "}$ [kpc]",
                            "dM/dt [M$_\odot$]",
                            "M$_{Halo}$ [M$_\odot$]",
                            "R$_{Halo}$ [pc]",
                            "q$_z$",
                            "t$_{past}$ [Myr]",
                            "R$_{sun}$ [pc]",
                            "v$_{LSR}$ [km/s]",
                            "$\Sigma_x$",
                            "$\Sigma_v$",
                            "$\Sigma_{vx}$",
                            "a$_{" + streamname + "}$ [pc/Myr$^2$]",
                            "V$_{sun}$ [km/s]",
                            "M$_{200}$ [M$_\odot$]",
                            "R$_{200}$ [pc]", "c"],
                    quantiles=[0.16, 0.84],
                    plot_datapoints=False,
                    dpi=150,
                    verbose=False,
                    truth_color='g',
                    plot_contours=True)

    filename = folder + 'chain%s_triangle' % name
    plt.savefig(filename + '.png', format='png', dpi=150, bbox_inches='tight')


def plot_triangle_orbital(name, streamname='stream', folder='./'):
    # Setup
    txt = 11
    mpl.rcParams['axes.labelsize'] = txt
    mpl.rcParams['xtick.labelsize'] = txt
    mpl.rcParams['ytick.labelsize'] = txt

    # Choice of stream type
    prob,\
        mass_cluster,\
        pm_mu_delta,\
        pm_mu_alphacosdelta,\
        mass_halo,\
        distance_cluster,\
        mass_loss_rate,\
        q_halo,\
        r_halo,\
        tpast,\
        sigma_x,\
        sigma_v,\
        sigma_vx,\
        sigma_mu,\
        rgalsun,\
        vLSR = np.loadtxt("chain%s.dat" % name, unpack=True)

        # Plot histograms of acceleration parameters

    good_inds = [good_vals(vals) for vals in [prob,
                                              rgalsun,
                                              mass_halo,
                                              q_halo,
                                              r_halo,
                                              distance_cluster,
                                              pm_mu_delta,
                                              pm_mu_alphacosdelta,
                                              vLSR]]
    ind = good_inds[0]
    for i in xrange(1, len(good_inds)):
        ind = np.intersect1d(good_inds[i], ind)

    prob = prob[ind]
    rgalsun = rgalsun[ind]
    mass_halo = mass_halo[ind]
    q_halo = q_halo[ind]
    r_halo = r_halo[ind]
    distance_cluster = distance_cluster[ind]
    pm_mu_delta = pm_mu_delta[ind]
    pm_mu_alphacosdelta = pm_mu_alphacosdelta[ind]
    vLSR = vLSR[ind]

    distance_cluster = distance_cluster * 1000.0

    samples = np.array([pm_mu_delta,
                        distance_cluster,
                        rgalsun,
                        vLSR,
                        pm_mu_alphacosdelta]).T

    extent = [(np.percentile(sample, 16),
               np.percentile(sample, 84)) for sample in samples.T]
    for i in xrange(len(extent)):
        if extent[i][0] == extent[i][1]:
            extent[i] = (extent[i][0] - 1, extent[i][0] + 1)

    fig, axes = plt.subplots(5, 5, figsize=[10, 10])

    triangle.corner(samples,
                    labels=["$\mu_\delta$ [arcsec/yr]",
                            "d$_{" + streamname + "}$ [kpc]",
                            "R$_{sun}$ [pc]",
                            "v$_{LSR}$ [km/s]",
                            "$\mu_\\alpha\cos(\delta)$ [arcsec/yr]"],
                    extents=extent,
                    quantiles=[0.16, 0.84],
                    plot_datapoints=False,
                    dpi=150,
                    verbose=False,
                    truth_color='g',
                    plot_contours=True,
                    fig=fig)

    axes[4, 0].errorbar(0.0008, -0.00465,
                        xerr=0.00082, yerr=0.00082,
                        color='red')

    filename = folder + 'chain%s_triangle_orbital' % name
    plt.savefig(filename + '.png', format='png', dpi=150, bbox_inches='tight')
    plt.close()


def plot_triangle_halo(name, folder='./'):
    # Setup
    txt = 11
    mpl.rcParams['axes.labelsize'] = txt
    mpl.rcParams['xtick.labelsize'] = txt
    mpl.rcParams['ytick.labelsize'] = txt

    # Choice of stream type
    prob,\
        mass_cluster,\
        pm_mu_delta,\
        pm_mu_alphacosdelta,\
        mass_halo,\
        distance_cluster,\
        mass_loss_rate,\
        q_halo,\
        r_halo,\
        tpast,\
        sigma_x,\
        sigma_v,\
        sigma_vx,\
        sigma_mu,\
        rgalsun,\
        vLSR = np.loadtxt("chain%s.dat" % name, unpack=True)

        # Plot histograms of acceleration parameters

    good_inds = [good_vals(vals) for vals in [prob,
                                              rgalsun,
                                              mass_halo,
                                              q_halo,
                                              r_halo,
                                              distance_cluster]]
    ind = good_inds[0]
    for i in xrange(1, len(good_inds)):
        ind = np.intersect1d(good_inds[i], ind)

    prob = prob[ind]
    rgalsun = rgalsun[ind]
    mass_halo = mass_halo[ind]
    q_halo = q_halo[ind]
    r_halo = r_halo[ind]
    distance_cluster = distance_cluster[ind]

    samples = np.array([mass_halo, q_halo, r_halo]).T

    extent = [(np.percentile(sample, 16),
               np.percentile(sample, 84)) for sample in samples.T]
    for i in xrange(len(extent)):
        if extent[i][0] == extent[i][1]:
            extent[i] = (extent[i][0] - 1, extent[i][0] + 1)

    plt.close()
    fig, axes = plt.subplots(3, 3, figsize=[10, 10])

    triangle.corner(samples,
                    labels=["M$_{Halo}$ [M$_\odot$]",
                            "q$_z$", "R$_{Halo}$ [pc]"],
                    extents=extent,
                    quantiles=[0.16, 0.84],
                    plot_datapoints=False,
                    dpi=150,
                    verbose=False,
                    truth_color='g',
                    plot_contours=True,
                    fig=fig)
    axes[1, 0].plot([1.1e+12, 1.1e+12],
                    [extent[1][0],
                     extent[1][1]],
                    linewidth=5,
                    color='red',
                    alpha=0.5)
    axes[2, 0].plot([1.1e+12, 1.1e+12],
                    [extent[2][0],
                     extent[2][1]],
                    linewidth=5,
                    color='red',
                    alpha=0.5)
    filename = folder + 'chain%s_triangle_halo' % name
    plt.savefig(filename + '.png', format='png', dpi=150, bbox_inches='tight')


def plot_triangle_NGC5466_nomass(name, folder='./'):
    # Setup
    txt = 11
    mpl.rcParams['axes.labelsize'] = txt
    mpl.rcParams['xtick.labelsize'] = txt
    mpl.rcParams['ytick.labelsize'] = txt

    # Choice of stream type
    prob,\
        mass_cluster,\
        pm_mu_delta,\
        pm_mu_alphacosdelta,\
        mass_halo,\
        distance_cluster,\
        mass_loss_rate,\
        q_halo,\
        r_halo,\
        tpast,\
        sigma_x,\
        sigma_v,\
        sigma_vx,\
        sigma_mu,\
        rgalsun,\
        vLSR = np.loadtxt("chain%s.dat" % name, unpack=True)

        # Plot histograms of acceleration parameters

    good_inds = [good_vals(vals) for vals in [pm_mu_alphacosdelta,
                                              pm_mu_delta,
                                              q_halo]]
    ind = good_inds[0]
    for i in xrange(1, len(good_inds)):
        ind = np.intersect1d(good_inds[i], ind)

    q_halo = q_halo[ind]
    pm_mu_alphacosdelta = pm_mu_alphacosdelta[ind]
    pm_mu_delta = pm_mu_delta[ind]

    samples = np.array([pm_mu_delta, q_halo, pm_mu_alphacosdelta]).T
    #samples = np.array([pm_mu_delta, q_halo, pm_mu_alphacosdelta]).T
    extent = [(np.percentile(sample, 5),
               np.percentile(sample, 95)) for sample in samples.T]
    for i in xrange(len(extent)):
        if extent[i][0] == extent[i][1]:
            extent[i] = (extent[i][0] - 1, extent[i][0] + 1)

    plt.close()
    fig, axes = plt.subplots(3, 3, figsize=[10, 10])
    #fig, axes = plt.subplots(3, 3, figsize=[10, 10])
    triangle.corner(samples,
                    labels=["$\mu_\delta$ [arcsec/yr]",
                            "q$_z$",
                            "$\mu_\\alpha\cos(\delta)$ [arcsec/yr]"],
                    extents=extent,
                    quantiles=[0.16, 0.84],
                    plot_datapoints=False,
                    dpi=150,
                    verbose=False,
                    truth_color='g',
                    plot_contours=True,
                    fig=fig)

    circ = plt.Circle((0.0008, -0.00465), radius=0.003,
                      alpha=0.5, color='red', fill=False)
    axes[2, 0].add_artist(circ)

    axes[2, 0].errorbar(0.0008, -0.00465,
                        xerr=0.00082, yerr=0.00082,
                        color='red')

    filename = folder + 'chain%s_triangle_NGC5466' % name
    plt.savefig(filename + '.png', format='png', dpi=150, bbox_inches='tight')


def plot_triangle_NGC5466(name, folder='./'):
    # Setup
    txt = 11
    mpl.rcParams['axes.labelsize'] = txt
    mpl.rcParams['xtick.labelsize'] = txt
    mpl.rcParams['ytick.labelsize'] = txt

    # Choice of stream type
    prob,\
        mass_cluster,\
        pm_mu_delta,\
        pm_mu_alphacosdelta,\
        mass_halo,\
        distance_cluster,\
        mass_loss_rate,\
        q_halo,\
        r_halo,\
        tpast,\
        sigma_x,\
        sigma_v,\
        sigma_vx,\
        sigma_mu,\
        rgalsun,\
        vLSR = np.loadtxt("chain%s.dat" % name, unpack=True)

        # Plot histograms of acceleration parameters

    good_inds = [good_vals(vals) for vals in [pm_mu_alphacosdelta,
                                              pm_mu_delta,
                                              q_halo,
                                              tpast,
                                              mass_cluster]]
    ind = good_inds[0]
    for i in xrange(1, len(good_inds)):
        ind = np.intersect1d(good_inds[i], ind)

    q_halo = q_halo[ind]
    mass_cluster = mass_cluster[ind]
    pm_mu_alphacosdelta = pm_mu_alphacosdelta[ind]
    pm_mu_delta = pm_mu_delta[ind]
    tpast = tpast[ind]

    samples = np.array([pm_mu_delta, mass_cluster,
                        q_halo, tpast, pm_mu_alphacosdelta]).T
    #samples = np.array([pm_mu_delta, q_halo, pm_mu_alphacosdelta]).T
    extent = [(np.percentile(sample, 5),
               np.percentile(sample, 95)) for sample in samples.T]
    for i in xrange(len(extent)):
        if extent[i][0] == extent[i][1]:
            extent[i] = (extent[i][0] - 1, extent[i][0] + 1)

    plt.close()
    fig, axes = plt.subplots(5, 5, figsize=[10, 10])
    #fig, axes = plt.subplots(3, 3, figsize=[10, 10])
    triangle.corner(samples,
                    labels=["$\mu_\delta$ [arcsec/yr]",
                            "$M_{NGC5466} [M_\odot]$",
                            "q$_z$",
                            "t$_{past}$ [Myr]",
                            "$\mu_\\alpha\cos(\delta)$ [arcsec/yr]"],
                    extents=extent,
                    quantiles=[0.16, 0.84],
                    plot_datapoints=False,
                    dpi=150,
                    verbose=False,
                    truth_color='g',
                    plot_contours=True,
                    fig=fig)

    circ = plt.Circle((0.0008, -0.00465), radius=0.003,
                      alpha=0.5, color='red', fill=False)
    axes[4, 0].add_artist(circ)

    axes[4, 0].errorbar(0.0008, -0.00465,
                        xerr=0.00082, yerr=0.00082,
                        color='red')

    filename = folder + 'chain%s_triangle_NGC5466' % name
    plt.savefig(filename + '.png', format='png', dpi=150, bbox_inches='tight')


def plot_triangle_LM2010(name, folder='./'):
    # Setup
    txt = 11
    mpl.rcParams['axes.labelsize'] = txt
    mpl.rcParams['xtick.labelsize'] = txt
    mpl.rcParams['ytick.labelsize'] = txt

    # Choice of stream type
    prob,\
        mass_cluster,\
        pm_mu_delta,\
        pm_mu_alphacosdelta,\
        mass_halo,\
        distance_cluster,\
        mass_loss_rate,\
        q_halo,\
        r_halo,\
        tpast,\
        sigma_x,\
        sigma_v,\
        sigma_vx,\
        sigma_mu,\
        rgalsun,\
        vLSR = np.loadtxt("chain%s.dat" % name, unpack=True)

        # Plot histograms of acceleration parameters

    good_inds = [good_vals(vals) for vals in [pm_mu_alphacosdelta,
                                              pm_mu_delta,
                                              mass_cluster,
                                              q_halo,
                                              tpast,
                                              mass_halo]]
    ind = good_inds[0]
    for i in xrange(1, len(good_inds)):
        ind = np.intersect1d(good_inds[i], ind)

    pm_mu_alphacosdelta = pm_mu_alphacosdelta[ind]
    pm_mu_delta = pm_mu_delta[ind]
    mass_cluster = mass_cluster[ind]
    mass_halo = mass_halo[ind]
    tpast = tpast[ind]
    q_halo = q_halo[ind]

    samples = np.array([pm_mu_delta, mass_cluster, mass_halo,
                        tpast, q_halo, pm_mu_alphacosdelta]).T

    extent = [(np.percentile(sample, 16),
               np.percentile(sample, 84)) for sample in samples.T]
    for i in xrange(len(extent)):
        if extent[i][0] == extent[i][1]:
            extent[i] = (extent[i][0] - 1, extent[i][0] + 1)

    plt.close()
    fig, axes = plt.subplots(6, 6, figsize=[20, 20])

    triangle.corner(samples,
                    labels=["$\mu_\delta$ [arcsec/yr]",
                            "M$_{NGC5466}$ [M$_{\odot}$]",
                            "Vhalo [km/s]",
                            "t$_{past}$ [Myr]",
                            "q$_z$",
                            "$\mu_\\alpha\cos(\delta)$ [arcsec/yr]"],
                    extents=extent,
                    quantiles=[0.16, 0.84],
                    plot_datapoints=False,
                    dpi=150,
                    verbose=False,
                    truth_color='g',
                    plot_contours=True,
                    fig=fig)

    circ = plt.Circle((0.0008, -0.00465), radius=0.003,
                      alpha=0.5, color='red', fill=False)
    axes[5, 0].add_artist(circ)

    axes[5, 0].errorbar(0.0008, -0.00465,
                        xerr=0.00082, yerr=0.00082,
                        color='red')

    filename = folder + 'chain%s_triangle_NGC5466_dist' % name
    plt.savefig(filename + '.png', format='png', dpi=150, bbox_inches='tight')


def plot_triangle_cluster(name, streamname='stream', folder='./'):
    # Setup
    txt = 11
    mpl.rcParams['axes.labelsize'] = txt
    mpl.rcParams['xtick.labelsize'] = txt
    mpl.rcParams['ytick.labelsize'] = txt

    # Choice of stream type
    prob,\
        mass_cluster,\
        pm_mu_delta,\
        pm_mu_alphacosdelta,\
        mass_halo,\
        distance_cluster,\
        mass_loss_rate,\
        q_halo,\
        r_halo,\
        tpast,\
        sigma_x,\
        sigma_v,\
        sigma_vx,\
        sigma_mu,\
        rgalsun,\
        vLSR = np.loadtxt("chain%s.dat" % name, unpack=True)

    good_inds = [good_vals(vals) for vals in [prob,
                                              rgalsun,
                                              mass_halo,
                                              mass_cluster,
                                              q_halo,
                                              r_halo,
                                              distance_cluster,
                                              tpast,
                                              mass_loss_rate]]
    ind = good_inds[0]
    for i in xrange(1, len(good_inds)):
        ind = np.intersect1d(good_inds[i], ind)

    prob = prob[ind]
    rgalsun = rgalsun[ind]
    mass_halo = mass_halo[ind]
    q_halo = q_halo[ind]
    r_halo = r_halo[ind]
    distance_cluster = distance_cluster[ind]
    mass_cluster = mass_cluster[ind]
    tpast = tpast[ind]
    mass_loss_rate = mass_loss_rate[ind]

    samples = np.array([mass_cluster, mass_loss_rate, tpast]).T

    extent = [(np.percentile(sample, 16), np.percentile(sample, 84))
              for sample in samples.T]
    for i in xrange(len(extent)):
        if extent[i][0] == extent[i][1]:
            extent[i] = (extent[i][0] - 1, extent[i][0] + 1)

    fig, axes = plt.subplots(3, 3, figsize=[10, 10])

    triangle.corner(samples,
                    labels=["M$_{" + streamname + "}$ [M$_\odot$]",
                            "dM/dt [M$_\odot$]",
                            "t$_{past}$ [Myr]"],
                    extents=extent,
                    quantiles=[0.16, 0.84],
                    plot_datapoints=False,
                    dpi=150, verbose=False,
                    truth_color='g',
                    plot_contours=True,
                    fig=fig)

    filename = folder + 'chain%s_triangle_cluster' % name
    plt.savefig(filename + '.png', format='png', dpi=150, bbox_inches='tight')
    plt.close()


def plot_triangle_derived(name,
                          streamname='stream',
                          gal_latitude=42.15,
                          gal_longitude=73.59,
                          folder='./'):
    # Setup
    txt = 11
    mpl.rcParams['axes.labelsize'] = txt
    mpl.rcParams['xtick.labelsize'] = txt
    mpl.rcParams['ytick.labelsize'] = txt

    # Choice of stream type
    prob,\
        mass_cluster,\
        pm_mu_delta,\
        pm_mu_alphacosdelta,\
        mass_halo,\
        distance_cluster,\
        mass_loss_rate,\
        q_halo,\
        r_halo,\
        tpast,\
        sigma_x,\
        sigma_v,\
        sigma_vx,\
        sigma_mu,\
        rgalsun,\
        vLSR = np.loadtxt("chain%s.dat" % name, unpack=True)

    good_inds = [good_vals(vals) for vals in [prob,
                                              rgalsun,
                                              mass_halo,
                                              q_halo,
                                              r_halo,
                                              distance_cluster]]
    ind = good_inds[0]
    for i in xrange(1, len(good_inds)):
        ind = np.intersect1d(good_inds[i], ind)

    prob = prob[ind]
    rgalsun = rgalsun[ind]
    mass_halo = mass_halo[ind]
    q_halo = q_halo[ind]
    r_halo = r_halo[ind]
    distance_cluster = distance_cluster[ind]

    distance_cluster = distance_cluster * 1000.0
    Vsun,\
        R200,\
        cNFW,\
        M200,\
        aPal5 = calculate_halo_properties(prob,
                                          rgalsun,
                                          mass_halo,
                                          q_halo,
                                          r_halo,
                                          distance_cluster,
                                          gal_latitude,
                                          gal_longitude)

    samples = np.array([M200, cNFW, aPal5, Vsun, R200]).T

    extent = [(np.percentile(sample, 16), np.percentile(sample, 84))
              for sample in samples.T]
    for i in xrange(len(extent)):
        if extent[i][0] == extent[i][1]:
            extent[i] = (extent[i][0] - 1, extent[i][0] + 1)

    plt.close()
    triangle.corner(samples,
                    labels=["M$_{200}$ [M$_\odot$]",
                            "c", "a$_{" + streamname + "}$ [pc/Myr$^2$]",
                            "V$_{sun}$ [km/s]", "R$_{200}$ [pc]"],
                    extents=extent,
                    quantiles=[0.16, 0.84],
                    plot_datapoints=False,
                    dpi=150, verbose=False,
                    truth_color='g', plot_contours=True)

    filename = folder + 'chain%s_triangle_derived' % name
    plt.savefig(filename + '.png', format='png', dpi=150, bbox_inches='tight')
