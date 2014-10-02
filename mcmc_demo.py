import numpy as np
import matplotlib.pyplot as plt
from myutils import *
import triangle

folder = '/Users/spardy/Research/projects/ISIMA/presentation/'

fig, axes = plt.subplots(2, 2, figsize=[10, 10])

parmA = np.random.randn(1000) + np.random.randn(1000) * 0.1
parmB = np.random.randn(1000) + np.random.randn(1000) * 0.1

samples = np.array([parmA, parmB]).T

triangle.corner(samples,
                labels=["Parameter A",
                        "Parameter B"],
                quantiles=[0.16, 0.84],
                plot_datapoints=False,
                dpi=150,
                verbose=False,
                truth_color='g',
                plot_contours=True,
                fig=fig)

filename = folder + 'MCMC_demo'
plt.savefig(filename + '.png', format='png', dpi=150, bbox_inches='tight')
plt.clf()
plt.close()
