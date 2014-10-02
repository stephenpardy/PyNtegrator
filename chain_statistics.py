import numpy as np
pma = np.zeros((4, 2))
pmd = np.zeros((4, 2))
qs = np.zeros((4, 2))
for j, i in enumerate([50, 100, 150, 200]):
    name = 'fixed_m%d' % i
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
    pma[j, 0] = np.median(pm_mu_alphacosdelta)
    pma[j, 1] = (np.percentile(pm_mu_alphacosdelta, 14)
                 - np.percentile(pm_mu_alphacosdelta, 84))
    pmd[j, 0] = np.median(pm_mu_delta)
    pmd[j, 1] = (np.percentile(pm_mu_delta, 14)
                 - np.percentile(pm_mu_delta, 84))
    qs[j, 0] = np.median(q_halo)
    qs[j, 1] = (np.percentile(q_halo, 14)
                - np.percentile(q_halo, 84))


print(pma)
print(pmd)
print(qs)

def mad(data, axis=None):
    return np.mean(np.absolute(data - np.mean(data, axis)), axis)
