import matplotlib as mpl
mpl.use("agg")
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

from myutils import *
import triangle
from mcorbit_utils import *

fig, axes = plt.subplots(2, 2)
#ax = plt.axes(xlim=(0, 2), ylim=(-2, 2))
#line, = ax.plot([], [], lw=2)

name = 'weight2_alldata_July28'

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
    rgalsun,\
    vLSR = np.loadtxt("chain%s.dat" % name, unpack=True)
    # Plot histograms of acceleration parameters
good_inds = [good_vals(vals) for vals in [q_halo, r_halo]]
ind = good_inds[0]
for i in xrange(1, len(good_inds)):
    ind = np.intersect1d(good_inds[i], ind)
q_halo = q_halo[ind]
r_halo = r_halo[ind]
samples = np.array([q_halo, r_halo]).T
extent = [(np.min(sample), np.max(sample)) for sample in samples.T]
for i in xrange(len(extent)):
    if extent[i][0] == extent[i][1]:
        extent[i] = (extent[i][0] - 1, extent[i][0] + 1)

for i in xrange(1000):
    fig, axes = plt.subplots(2, 2)

    samples = np.array([q_halo[0:256 * (i + 1)],
                        r_halo[0:256 * (i + 1)]]).T
    triangle.corner(samples,
                    labels=["q$_z$", "R$_{Halo}$ [pc]"],
                    extents=extent,
                    plot_datapoints=False,
                    dpi=150,
                    vebose=False,
                    truth_color='g',
                    plot_contours=True,
                    fig=fig)
    plt.savefig('triangle' + str(i).zfill(4) + '.png',
                format='png',
                dpi=150,
                bbox_inches='tight')
    plt.clf()
    plt.close()


def init():
    #line.set_data([], [])
    global q_halo, r_halo, extent
    samples = np.array([q_halo[0:100], r_halo[0:100]]).T
    triangle.corner(samples,
                    labels=["q$_z$", "R$_{Halo}$ [pc]"],
                    extents=extent,
                    plot_datapoints=False,
                    dpi=150,
                    vebose=False,
                    truth_color='g',
                    plot_contours=True,
                    fig=fig)
    #return line,


def animate(i):
    #x = np.linspace(0, 2, 1000)
    #y = np.sin(2 * np.pi * (x - 0.01 * i))
    #line.set_data(x, y)
    #return line,
    samples = np.array([q_halo[0:100 * (i + 1)],
                        r_halo[0:100 * (i + 1)]]).T
    triangle.corner(samples,
                    labels=["q$_z$", "R$_{Halo}$ [pc]"],
                    extents=extent,
                    plot_datapoints=False,
                    dpi=150,
                    vebose=False,
                    truth_color='g',
                    plot_contours=True,
                    fig=fig)

# call the animator.  blit=True means only re-draw the parts that have changed.
#anim = animation.FuncAnimation(fig, animate, init_func=init,
#                               frames=200, interval=20, blit=False)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
#anim.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
