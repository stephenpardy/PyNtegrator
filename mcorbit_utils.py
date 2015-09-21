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
