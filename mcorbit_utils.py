import numpy as np
from math import *
import sys
import scipy.optimize
import inspect


def printcol(*arg, **kwarg):
    # Print vectors in columns
    # Use: printcol <vec1> <vec2> .. <vecn> (<fout='path to file'>)
    # Default: fout=sys.stdout
    #

    # Set output
    if kwarg:
        f = open(kwarg['fout'], 'w')
    else:
        f = sys.stdout

    # Get variable names
    frame = inspect.currentframe()
    frame2 = inspect.getouterframes(frame)[1]
    string = inspect.getframeinfo(frame2[0]).code_context[0].strip()
    args = string[string.find('(') + 1:-1].split(',')

    names = []
    for i in args:
        if i.find('=') != -1:
            names.append(i.split('=')[1].strip())
        else:
            names.append(i)

    Ncol = len(arg)
    Nrow = np.zeros(Ncol)

    for i in range(Ncol):
        Nrow[i] = len(arg[i])

    Nmax = int(np.max(Nrow))

    # Print
    print>>f, ("#"),
    for i in range(len(names)):
        print>>f, ("%s\t" % names[i]),
    print>>f, ("\n#\n"),

    for i in range(Nmax):
        for j in range(Ncol):
            if i < Nrow[j]:
                print>>f, ('%g\t' % arg[j][i]),
            else:
                print>>f, ('\t'),
        print>>f, ('\n'),

    # print Ncol, Nrow
    if kwarg:
        f.close()


def load_OD(filename):
    # n[][0] = ra*
    # n[][1] = dec
    # n[][3] = significance/weight
    # n[][4] = FWHM
    (ra, dec, weights, FWHM) = np.loadtxt(filename, unpack=True, delimiter=",")
#{224, 21.0, 0.00000,  2.0000,  1.0000}

def load_VR(filename):
    (vr, dvr, ra, dec, rastar, mura, mudel, mura_e, mudel_e) = np.loadtxt(filename, unpack=True, delimiter=",")
    # Helper to load radial velocity and propermotion data
      # n_vr[][0] = vr
      # n_vr[][1] = dvr
      # n_vr[][2] = l
      # n_vr[][3] = b
      # n_vr[][4] = l*
      # n_vr[][5] = #
      # n_vr[][6] = mura
      # n_vr[][7] = mudel
      # n_vr[][8] = mura_e
      # n_vr[][9] = mudel_e


def get_pool(mpi=False, threads=None):
    """ Get a pool object to pass to emcee for parallel processing.
        If mpi is False and threads is None, pool is None.

        Parameters
        ----------
        mpi : bool
        Use MPI or not. If specified, ignores the threads kwarg.
        threads : int (optional)
        If mpi is False and threads is specified, use a Python
        multiprocessing pool with the specified number of threads.
        """

    if mpi:
        from emcee.utils import MPIPool

        # Initialize the MPI pool
        pool = MPIPool()

        # Make sure the thread we're running on is the master
        if not pool.is_master():
            pool.wait()
            sys.exit(0)
        print("Running with MPI...")

    elif threads > 1:
        import multiprocessing
        print("Running with multiprocessing on " + str(threads) + " cores...")
        pool = multiprocessing.Pool(threads)

    else:
        print("Running serial...")
        pool = None

    return pool


def unbound(pmra, pmdec):
    pmra *= 1E3
    pmdec *= 1E3
    vrad = 106.7
    return np.sqrt((251.74 +
                    56.32 * pmdec +
                    48.0127 * pmra +
                    0.189495 * vrad) ** 2 +
                   (11.1 -
                    50.0918 * pmdec +
                    54.0621 * pmra +
                    0.209465 * vrad) ** 2 +
                   (7.25 -
                    0.187537 * pmdec -
                    21.2892 * pmra +
                    0.959279 * vrad) ** 2) > 600.


def good_vals(vals):
    """
    Find very bad values
    """
    large_val = 1e+36
    w = np.where(np.abs(vals) < large_val)[0]
    if len(w) > 0:
        return w
    else:
        return 0


# Custom function for NGC5466 properties
def calculate_halo_properties(prob,
                              rgalsun,
                              mass_halo,
                              q_halo,
                              r_halo,
                              distance_cluster,
                              gal_latitude,
                              gal_longitude):

    step = np.arange(len(prob))

        # Plot histograms of acceleration parameters
    Vsun = 1.0 * np.arange(len(prob))
    aStream = 1.0 * np.arange(len(prob))
    R200 = 1.0 * np.arange(len(prob))
    M200 = 1.0 * np.arange(len(prob))
    cNFW = 1.0 * np.arange(len(prob))

    distance_cluster = distance_cluster * 1000.0

    for i in step:
        x = -rgalsun[i]
        y = 0.0
        z = 0.0
        ax = 0.0
        ay = 0.0
        az = 0.0
        G = 0.0043009211  # gravitational constant in [km^2/s^2/Msun*pc]
        # Law, Majewski & Johnston (2009) potential constants
        b1_LMJ = 700.0  # [pc]
        M1_LMJ = 3.4e10  # [solar masses]
        a2_LMJ = 6500.0  # [pc]
        b2_LMJ = 260.0  # [pc]
        M2_LMJ = 1.0e11  # [solar masses]

        # Hernquist bulge
        r = sqrt(x * x + y * y + z * z)

        a1x = -G * M1_LMJ / ((r + b1_LMJ) * (r + b1_LMJ)) * x / r
        a1y = -G * M1_LMJ / ((r + b1_LMJ) * (r + b1_LMJ)) * y / r
        a1z = -G * M1_LMJ / ((r + b1_LMJ) * (r + b1_LMJ)) * z / r

        # Miyamato disk
        r2 = sqrt(x * x + y * y + (a2_LMJ + sqrt(z * z + b2_LMJ * b2_LMJ))
                  * (a2_LMJ + sqrt(z * z + b2_LMJ * b2_LMJ)))

        a2x = -G * M2_LMJ / (r2 * r2 * r2) * x
        a2y = -G * M2_LMJ / (r2 * r2 * r2) * y
        a2z = -G * M2_LMJ / (r2 * r2 * r2) *\
            (a2_LMJ + sqrt(z * z + b2_LMJ * b2_LMJ))\
            / sqrt(z * z + b2_LMJ * b2_LMJ) * z

        # NFW Halo
        r3 = sqrt(x * x + y * y + z / q_halo[i] * z / q_halo[i])

        a3x = -G * mass_halo[i] / r3 * (log(1.0 + r3 / r_halo[i])
                                        / r3 - 1.0 / (r_halo[i] + r3)) * x / r3
        a3y = -G * mass_halo[i] / r3 * (log(1.0 + r3 / r_halo[i])
                                        / r3 - 1.0 / (r_halo[i] + r3)) * y / r3
        a3z = -G * mass_halo[i] / r3 * (log(1.0 + r3 / r_halo[i])
                                        / r3 - 1.0 / (r_halo[i] + r3))\
            * z / (q_halo[i] * q_halo[i] * r3)

        ax = a1x + a2x + a3x
        ay = a1y + a2y + a3y
        az = a1z + a2z + a3z

        Vsun[i] = sqrt(r * sqrt(ax * ax + ay * ay + az * az))

        xsun = -rgalsun[i]

        # Cluster acceleration
        brad = gal_latitude / 360.0 * 2.0 * np.pi
        # Cluster galactic latitude [rad]
        lrad = gal_longitude / 360.0 * 2.0 * np.pi
        # Cluster galactic longitude [deg]

        z = sin(brad) * distance_cluster[i]
        dxy = sqrt(distance_cluster[i] * distance_cluster[i] - z * z)
        x = cos(lrad) * dxy + xsun
        y = sin(lrad) * dxy

        # Hernquist bulge
        r = sqrt(x * x + y * y + z * z)

        a1x = -G * M1_LMJ / ((r + b1_LMJ) * (r + b1_LMJ)) * x / r
        a1y = -G * M1_LMJ / ((r + b1_LMJ) * (r + b1_LMJ)) * y / r
        a1z = -G * M1_LMJ / ((r + b1_LMJ) * (r + b1_LMJ)) * z / r

        # Miyamato disk
        r2 = sqrt(x * x + y * y + (a2_LMJ + sqrt(z * z + b2_LMJ * b2_LMJ))
                  * (a2_LMJ + sqrt(z * z + b2_LMJ * b2_LMJ)))

        a2x = -G * M2_LMJ / (r2 * r2 * r2) * x
        a2y = -G * M2_LMJ / (r2 * r2 * r2) * y
        a2z = -G * M2_LMJ / (r2 * r2 * r2) *\
            (a2_LMJ + sqrt(z * z + b2_LMJ * b2_LMJ))\
            / sqrt(z * z + b2_LMJ * b2_LMJ) * z

        # NFW Halo
        r3 = sqrt(x * x + y * y + z / q_halo[i] * z / q_halo[i])

        a3x = -G * mass_halo[i] / r3 * (log(1.0 + r3 / r_halo[i])
                                        / r3 - 1.0 / (r_halo[i] + r3)) * x / r3
        a3y = -G * mass_halo[i] / r3 * (log(1.0 + r3 / r_halo[i])
                                        / r3 - 1.0 / (r_halo[i] + r3)) * y / r3
        a3z = -G * mass_halo[i] / r3 * (log(1.0 + r3 / r_halo[i])
                                        / r3 - 1.0 / (r_halo[i] + r3))\
            * z / (q_halo[i] * q_halo[i] * r3)

        ax = a1x + a2x + a3x
        ay = a1y + a2y + a3y
        az = a1z + a2z + a3z

        aStream[i] = 1.0 * sqrt(ax * ax + ay * ay + az * az)

        k = 1.3e-7  # rho_crit in Msun/pc^3

        def f(y):
            z = 3.0*mass_halo[i]/(8.0*np.pi*y**3)*(r_halo[i]
                                                   / (r_halo[i] + y)
                                                   + log((r_halo[i] + y)
                                                         / r_halo[i])
                                                   - 1.0) - 200.0*k
            return z

        cNFW[i] = R200[i] / r_halo[i]
        M200[i] = mass_halo[i] * (r_halo[i] / (r_halo[i] + R200[i])
                                  + log(r_halo[i] + R200[i])
                                  - 1.0 - log(r_halo[i]))

    return Vsun, R200, cNFW, M200, aStream
