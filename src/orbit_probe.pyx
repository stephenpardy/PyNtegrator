import pyorbits
import numpy as np
cimport numpy as np

def test_orbit(dict input_parameters,
               str gal_name,
               str gal_name2,
               np.ndarray[double, ndim=2, mode="c"] data_pos):
    """
    Test against a full orbit (i.e. Gadget or another simulation)
    Args:
        dict input_parameters
            See above
        str gal_name
            galaxy to compare with data
        np.array data_pos
            2d array with x, y, z for each snapshot for data
        np.array data_vel
            2d array with vx, vy, vz for each snapshot for data
    """
    cdef list results
    cdef double ln_likelihood
    cdef int nsnaps, g
    try:
        results = pyorbits.run(input_parameters)
    except RuntimeError, e:
        print({"Runtime Error: {:s}".format(e)})
        return None
    except KeyError, e:
        print("Problem in initial conditions. Missing parameter {:s}".format(e))
        return None

    cdef int ngals = len(input_parameters['galaxies'])

    nsnaps = len(results)
    g = np.where([results[0][i]['name'] == gal_name for i in xrange(ngals)])[0]
    g2 = np.where([results[0][i]['name'] == gal_name2 for i in xrange(ngals)])[0]

    model_pos = np.zeros((nsnaps, 3))
    for i in xrange(nsnaps):
        for j in xrange(3):
            model_pos[i, j] = results[i][g]['pos'][j]-results[i][g2]['pos'][j]

    ln_likelihood = pyorbits.likelihood2(model_pos,
                                data_pos,
                                input_parameters['pos_err'])


    return ln_likelihood


def test_location(dict input_parameters):
    """
    Test against the current location of galaxies.
    Args:
        dict input_parameters
        See above
    """
    cdef list results
    try:
        results = pyorbits.run(input_parameters)
    except RuntimeError, e:
        print({"Runtime Error: {:s}".format(e)})
        return None
    except KeyError, e:
        print("Problem in initial conditions. Missing parameter {:s}".format(e))
        return None

    cdef int ngals = len(input_parameters['galaxies'])
    cdef double ln_likelihood

    model_pos = np.empty((ngals, 3))
    model_vel = np.empty((ngals, 3))
    for i, gal in enumerate(results[-1]):
        model_pos[i, :] = gal['pos']
        model_vel[i, :] = gal['vel']

    data_pos = np.empty((ngals, 3))
    data_vel = np.empty((ngals, 3))
    for i, (galname, gal) in enumerate(input_parameters['galaxies'].iteritems()):
        data_pos[i, :] = gal['pos']
        data_vel[i, :] = gal['vel']


    ln_likelihood = pyorbits.likelihood(ngals,
                               model_pos,
                               model_vel,
                               data_pos,
                               data_vel,
                               input_parameters['pos_err'],
                               input_parameters['vel_err'])
    return ln_likelihood


def orbit_statistics(dict input_parameters,
                     str gal_name,
                     str gal_name2):

    cdef list results

    try:
        results = pyorbits.run(input_parameters)
    except RuntimeError, e:
        print("Runtime Error: {:s}".format(e))
        return None
    except KeyError, e:
        print("Problem in initial conditions. Missing parameter {:s}".format(e))
        return None

    return pyorbits.test_orbit_statistics(results, len(input_parameters['galaxies']), gal_name, gal_name2)


#Gives 2/3 passes
def test_stream(dict input_parameters,
                np.ndarray[double, ndim=2, mode="c"] input_position,
                bint use_tracers):

    from astropy.coordinates import SkyCoord # High-level coordinates
    import astropy.coordinates as coord

    cdef list results
    try:
        if use_tracers:
            results, tracers = pyorbits.run(input_parameters)
        else:
            results = pyorbits.run(input_parameters)
    except RuntimeError, e:
        print("Runtime Error: {:s}".format(e))
        return None
    except KeyError, e:
        print("Problem in initial conditions. Missing parameter {:s}".format(e))
        return None

    cdef float rvir = 200.0

    mindist, maxdist, apos, peris, stripped = pyorbits.test_orbit_statistics(results,
                                                                    len(input_parameters['galaxies']),
                                                                    "MW", "LMC",
                                                                    min_time=-6.0)

    # Must be a second passage model
    if (peris < 2):
        print "Pericenters: ", peris
        return np.NINF

    #Check that the maxdistance is greater than rvir, this ensures that we have entered the system with the last 6Gyr
    if (maxdist <= rvir):
        print "LMC Maxdist: {:5.5f}, with {:d} pericenters".format(maxdist, peris)
        return np.NINF

    mindist, maxdist, apos, peris, stripped = pyorbits.test_orbit_statistics(results,
                                                                    len(input_parameters['galaxies']),
                                                                    "LMC", "SMC",
                                                                    min_time=-2.0)


    #Check that the two clouds have had a recent close encounter
    if (mindist >= 10):
        print "SMC mindist: ", mindist
        return np.NINF

    cdef double ln_likelihood = 0.0
    cdef int ngals = len(input_parameters['galaxies'])

    #Test final position/velocity against starting position/velocity

    model_pos = np.empty((ngals, 3))
    model_vel = np.empty((ngals, 3))
    for i, gal in enumerate(results[-1]):
        model_pos[i, :] = gal['pos']
        model_vel[i, :] = gal['vel']

    data_pos = np.empty((ngals, 3))
    data_vel = np.empty((ngals, 3))
    for i, (galname, gal) in enumerate(input_parameters['galaxies'].iteritems()):
        data_pos[i, :] = gal['pos']
        data_vel[i, :] = gal['vel']

    ln_likelihood += pyorbits.likelihood(ngals,
                                model_pos,
                                model_vel,
                                data_pos,
                                data_vel,
                                2,  # kpc pos err,
                                5)  # kms vel err


    #Test that the orbit roughly follows the current location of the MS

    g = np.where([results[0][i]['name'] == "MW" for i in xrange(ngals)])[0]
    g2 = np.where([results[0][i]['name'] == "LMC" for i in xrange(ngals)])[0]
    model_pos = np.zeros((len(results), 3))
    for r in results:
        c = SkyCoord(w=r[g2]['pos'][0] - r[g]['pos'][0],
                     u=r[g2]['pos'][1] - r[g]['pos'][1],
                     v=r[g2]['pos'][2] - r[g]['pos'][2], unit='kpc',
                     frame='galactic', representation='cartesian')
        coords = c.transform_to(coord.Galactic)
        model_pos[i, 0] = coords.l.value
        model_pos[i, 1] = coords.b.value

    ln_likelihood += pyorbits.likelihood2(model_pos,
                                 input_position,
                                 0.1)  # degrees pos error

    return ln_likelihood


def test_stream_1stpass(dict input_parameters,
                np.ndarray[double, ndim=2, mode="c"] input_position,
                bint use_tracers):

    from astropy.coordinates import SkyCoord # High-level coordinates
    import astropy.coordinates as coord

    cdef list results
    try:
        if use_tracers:
            results, tracers = pyorbits.run(input_parameters)
        else:
            results = pyorbits.run(input_parameters)
    except RuntimeError, e:
        print("Runtime Error: {:s}".format(e))
        return None
    except KeyError, e:
        print("Problem in initial conditions. Missing parameter {:s}".format(e))
        return None

    cdef float rvir = 300.0

    mindist, maxdist, apos, peris, stripped = pyorbits.test_orbit_statistics(results,
                                                                    len(input_parameters['galaxies']),
                                                                    "MW", "LMC",
                                                                    min_time=-6.0)

    # Must be a first passage model
    if (peris > 1):
        print "Pericenters: ", peris
        return np.NINF

    #Check that the maxdistance is greater than rvir, this ensures that we have entered the system with the last 6Gyr
    if (maxdist <= rvir):
        print "LMC Maxdist: ", maxdist
        return np.NINF

    mindist, maxdist, apos, peris, stripped = pyorbits.test_orbit_statistics(results,
                                                                    len(input_parameters['galaxies']),
                                                                    "LMC", "SMC",
                                                                    min_time=-2.0)


    #Check that the two clouds have had a recent close encounter
    if (mindist >= 10):
        print "SMC mindist: ", mindist
        return np.NINF

    cdef double ln_likelihood = 0.0
    cdef int ngals = len(input_parameters['galaxies'])

    #Test final position/velocity against starting position/velocity

    model_pos = np.empty((ngals, 3))
    model_vel = np.empty((ngals, 3))
    for i, gal in enumerate(results[-1]):
        model_pos[i, :] = gal['pos']
        model_vel[i, :] = gal['vel']

    data_pos = np.empty((ngals, 3))
    data_vel = np.empty((ngals, 3))
    for i, (galname, gal) in enumerate(input_parameters['galaxies'].iteritems()):
        data_pos[i, :] = gal['pos']
        data_vel[i, :] = gal['vel']

    ln_likelihood += pyorbits.likelihood(ngals,
                                model_pos,
                                model_vel,
                                data_pos,
                                data_vel,
                                2,  # kpc pos err,
                                5)  # kms vel err


    #Test that the orbit roughly follows the current location of the MS

    g = np.where([results[0][i]['name'] == "MW" for i in xrange(ngals)])[0]
    g2 = np.where([results[0][i]['name'] == "LMC" for i in xrange(ngals)])[0]
    model_pos = np.zeros((len(results), 3))
    for r in results:
        c = SkyCoord(w=r[g2]['pos'][0] - r[g]['pos'][0],
                     u=r[g2]['pos'][1] - r[g]['pos'][1],
                     v=r[g2]['pos'][2] - r[g]['pos'][2], unit='kpc',
                     frame='galactic', representation='cartesian')
        coords = c.transform_to(coord.Galactic)
        model_pos[i, 0] = coords.l.value
        model_pos[i, 1] = coords.b.value

    ln_likelihood += pyorbits.likelihood2(model_pos,
                                 input_position,
                                 0.1)  # degrees pos error

    return ln_likelihood


def test_stream_2ndpass(dict input_parameters,
                np.ndarray[double, ndim=2, mode="c"] input_position,
                bint use_tracers):

    from astropy.coordinates import SkyCoord # High-level coordinates
    import astropy.coordinates as coord

    cdef list results
    try:
        if use_tracers:
            results, tracers = pyorbits.run(input_parameters)
        else:
            results = pyorbits.run(input_parameters)
    except RuntimeError, e:
        print("Runtime Error: {:s}".format(e))
        return None
    except KeyError, e:
        print("Problem in initial conditions. Missing parameter {:s}".format(e))
        return None

    cdef float rvir = 200.0

    mindist, maxdist, apos, peris, stripped = pyorbits.test_orbit_statistics(results,
                                                                    len(input_parameters['galaxies']),
                                                                    "MW", "LMC",
                                                                    min_time=-6.0)

    # Must be a second passage model (I think that some 2nd pass models are getting labeled as 1 pass)
    if not (0 < peris < 3):
        print "Pericenters: ", peris
        return np.NINF

    #Check that the maxdistance is greater than rvir, this ensures that we have entered the system with the last 6Gyr
    if (maxdist <= rvir):
        print "LMC Maxdist: ", maxdist
        return np.NINF

    mindist, maxdist, apos, peris, stripped = pyorbits.test_orbit_statistics(results,
                                                                    len(input_parameters['galaxies']),
                                                                    "LMC", "SMC",
                                                                    min_time=-2.0)


    #Check that the two clouds have had a recent close encounter
    if (mindist >= 15):
        print "SMC mindist: ", mindist
        return np.NINF

    cdef double ln_likelihood = 0.0
    cdef int ngals = len(input_parameters['galaxies'])

    #Test final position/velocity against starting position/velocity

    model_pos = np.empty((ngals, 3))
    model_vel = np.empty((ngals, 3))
    for i, gal in enumerate(results[-1]):
        model_pos[i, :] = gal['pos']
        model_vel[i, :] = gal['vel']




    data_pos = np.empty((ngals, 3))
    data_vel = np.empty((ngals, 3))
    for i, (galname, gal) in enumerate(input_parameters['galaxies'].iteritems()):
        data_pos[i, :] = gal['pos']
        data_vel[i, :] = gal['vel']

    ln_likelihood += pyorbits.likelihood(ngals,
                                model_pos,
                                model_vel,
                                data_pos,
                                data_vel,
                                2,  # kpc pos err,
                                5)  # kms vel err

    #Test that the orbit roughly follows the current location of the MS

    g = np.where([results[0][i]['name'] == "MW" for i in xrange(ngals)])[0]
    g2 = np.where([results[0][i]['name'] == "LMC" for i in xrange(ngals)])[0]
    model_pos = np.zeros((len(results), 3))
    for r in results:
        c = SkyCoord(w=r[g2]['pos'][0] - r[g]['pos'][0],
                     u=r[g2]['pos'][1] - r[g]['pos'][1],
                     v=r[g2]['pos'][2] - r[g]['pos'][2], unit='kpc',
                     frame='galactic', representation='cartesian')
        coords = c.transform_to(coord.Galactic)
        model_pos[i, 0] = coords.l.value
        model_pos[i, 1] = coords.b.value

    ln_likelihood += pyorbits.likelihood2(model_pos,
                                 input_position,
                                 0.1)  # degrees pos error


    return ln_likelihood


def test_stream_DB2012(dict input_parameters,
                np.ndarray[double, ndim=2, mode="c"] input_position,
                bint use_tracers):

    from astropy.coordinates import SkyCoord # High-level coordinates
    import astropy.coordinates as coord

    cdef list results
    try:
        if use_tracers:
            results, tracers = pyorbits.run(input_parameters)
        else:
            results = pyorbits.run(input_parameters)
    except RuntimeError, e:
        print("Runtime Error: {:s}".format(e))
        return None
    except KeyError, e:
        print("Problem in initial conditions. Missing parameter {:s}".format(e))
        return None

    mindist, maxdist, apos, peris, stripped = pyorbits.test_orbit_statistics(results,
                                                                    len(input_parameters['galaxies']),
                                                                    "MW", "LMC",
                                                                    min_time=-5.0)

    # Must be more than one passage
    if (peris < 2):
        print "Pericenters: ", peris
        return np.NINF

    #Test LMC-SMC conditions in DB2012
    cdef int nsnaps, g, g2, s
    cdef int ngals = len(input_parameters['galaxies'])

    nsnaps = len(results)
    g = np.where([results[0][i]['name'] == 'LMC' for i in xrange(ngals)])[0]
    g2 = np.where([results[0][i]['name'] == 'SMC' for i in xrange(ngals)])[0]
    old_dist = None
    direction = None

    cdef int apocenters = 0
    cdef int pericenters = 0

    # Roll these together with the loop below
    dist = np.array([np.sqrt((results[i][g]['pos'][0]-results[i][g2]['pos'][0])**2 +
                    (results[i][g]['pos'][1]-results[i][g2]['pos'][1])**2 +
                    (results[i][g]['pos'][2]-results[i][g2]['pos'][2])**2) for i in xrange(nsnaps)])

    stripped = np.array([results[i][g2]['stripped'] for i in xrange(nsnaps)])
    times = np.array([results[i][g2]['t'] for i in xrange(nsnaps)])

    idx = (times < 0.0)*(times > -3.0)
    a = dist[idx]
    #np.r_ is an indexing trick that lets you build arrays
    # this gets all locations where the array is at a minimum
    minimums = np.r_[True, a[1:] < a[:-1]] & np.r_[a[:-1] < a[1:], True]
    maximums = np.r_[True, a[1:] > a[:-1]] & np.r_[a[:-1] > a[1:], True]

    if a.mean() > 50:
        print "Average LMC-SMC separation:", a.mean()
        return np.NINF

    if sum(minimums) != 2:
        #not two close interactions
        print "LMC-SMC pericenters: ", sum(minimums)
        return np.NINF

    interaction_times = times[idx][minimums]

    if not (-2.25 < interaction_times[0] < -1.75):
        print "First interaction at: ", interaction_times[0]
        return np.NINF

    if not (-0.65 < interaction_times[1] < -0.15):
        print "Second interaction at:", interaction_times[1]
        return np.NINF

    if a.min() < 5:
        print "LMC-SMC minimum distance:", a.min()
        return np.NINF

    cdef double ln_likelihood = 0.0

    #Test final position/velocity against starting position/velocity

    model_pos = np.empty((ngals, 3))
    model_vel = np.empty((ngals, 3))
    for i, gal in enumerate(results[-1]):
        model_pos[i, :] = gal['pos']
        model_vel[i, :] = gal['vel']


    data_pos = np.empty((ngals, 3))
    data_vel = np.empty((ngals, 3))
    for i, (galname, gal) in enumerate(input_parameters['galaxies'].iteritems()):
        data_pos[i, :] = gal['pos']
        data_vel[i, :] = gal['vel']

    ln_likelihood += pyorbits.likelihood(ngals,
                                model_pos,
                                model_vel,
                                data_pos,
                                data_vel,
                                2,  # kpc pos err,
                                5)  # kms vel err

    return ln_likelihood


def test_stream_nosmc(dict input_parameters,
                      np.ndarray[double, ndim=2, mode="c"] input_position,
                      bint use_tracers):

    from astropy.coordinates import SkyCoord # High-level coordinates
    import astropy.coordinates as coord

    cdef list results
    try:
        if use_tracers:
            results, tracers = pyorbits.run(input_parameters)
        else:
            results = pyorbits.run(input_parameters)
    except RuntimeError, e:
        print("Runtime Error: {:s}".format(e))
        return None
    except KeyError, e:
        print("Problem in initial conditions. Missing parameter {:s}".format(e))
        return None

    cdef float rvir = 300.0

    mindist, maxdist, apos, peris, stripped = pyorbits.test_orbit_statistics(results,
                                                                    len(input_parameters['galaxies']),
                                                                    "MW", "LMC",
                                                                    min_time=-6.0)

    # Must be a second passage model
    if (peris < 2):
        print "Pericenters: ", peris
        return np.NINF

    #Check that the maxdistance is greater than rvir, this ensures that we have entered the system with the last 6Gyr
    if (maxdist <= rvir):
        print "LMC Maxdist: ", maxdist
        return np.NINF


    cdef double ln_likelihood = 0.0
    cdef int ngals = len(input_parameters['galaxies'])

    #Test final position/velocity against starting position/velocity

    model_pos = np.empty((ngals, 3))
    model_vel = np.empty((ngals, 3))
    for i, gal in enumerate(results[-1]):
        model_pos[i, :] = gal['pos']
        model_vel[i, :] = gal['vel']

    data_pos = np.empty((ngals, 3))
    data_vel = np.empty((ngals, 3))
    for i, (galname, gal) in enumerate(input_parameters['galaxies'].iteritems()):
        data_pos[i, :] = gal['pos']
        data_vel[i, :] = gal['vel']

    ln_likelihood += pyorbits.likelihood(ngals,
                                model_pos,
                                model_vel,
                                data_pos,
                                data_vel,
                                2,  # kpc pos err,
                                5)  # kms vel err


    #Test that the orbit roughly follows the current location of the MS

    g = np.where([results[0][i]['name'] == "MW" for i in xrange(ngals)])[0]
    g2 = np.where([results[0][i]['name'] == "LMC" for i in xrange(ngals)])[0]
    model_pos = np.zeros((len(results), 3))
    for r in results:
        c = SkyCoord(w=r[g2]['pos'][0] - r[g]['pos'][0],
                     u=r[g2]['pos'][1] - r[g]['pos'][1],
                     v=r[g2]['pos'][2] - r[g]['pos'][2], unit='kpc',
                     frame='galactic', representation='cartesian')
        coords = c.transform_to(coord.Galactic)
        model_pos[i, 0] = coords.l.value
        model_pos[i, 1] = coords.b.value

    ln_likelihood += pyorbits.likelihood2(model_pos,
                                 input_position,
                                 0.1)  # degrees pos error

    return ln_likelihood
