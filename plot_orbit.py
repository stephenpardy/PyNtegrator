import numpy as np
import matplotlib.pyplot as plt
from os import listdir

folder = './'


def get_snap(f):
    x, y, z, v = np.loadtxt(f, skiprows=1, delimiter=',')
    return x, y, z

dir_contents = listdir('./')
files = [f for f in dir_contents if 'cluster' in f]

fig, axes = plt.subplots(2, figsize=[10, 20], sharex=True)

pericenter = 1e+6
apocenter = 0.0

for f in files:
    x, y, z = get_snap(f)
    dist = np.sqrt(x**2 + y**2 + z**2)
    if dist < pericenter:
        pericenter = dist
    if dist > apocenter:
        apocenter = dist

    axes[0].plot(x, y, 'b.')
    #axes[0].set_xlabel("x [kpc]")
    axes[0].set_ylabel("y [kpc]", fontsize='larger')

    axes[1].plot(x, z, 'b.')
    axes[1].set_xlabel("x [kpc]", fontsize='larger')
    axes[1].set_ylabel("z [kpc]", fontsize='larger')


print(apocenter, pericenter)
filename = folder + 'NGC5466_orbit'
plt.savefig(filename + '.png', format='png', dpi=150, bbox_inches='tight')
plt.clf()
plt.close()
