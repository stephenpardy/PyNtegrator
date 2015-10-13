#encoding: utf-8
#cython: profile=False

import numpy as np
cimport numpy as np
import random
from libc.stdlib cimport malloc, free
import copy


cdef extern from 'orbit.c':
    int orbit(int ngals,
              Params parameters,
              Gal *gal,
              Snapshot **output_snapshots)

cdef extern from *:

    ctypedef struct Tracer:
        int nparticles
        double *pos
        double *vel

    ctypedef struct Gal:
        double pos[3]
        double vel[3]
        double post[3]
        double velt[3]
        int ID
        double mhalo
        double mtidal
        double r_halo
        double gamma
        double a2_LMJ
        double b2_LMJ
        double M2_LMJ
        double M1_LMJ
        double b1_LMJ
        int dyn_fric
        double dyn_C_eq
        double dyn_L_eq
        double dyn_alpha_eq
        double dyn_C_uneq
        double dyn_L_uneq
        double dyn_alpha_uneq
        int tidal_trunc
        int stripped
        double rt
        int inplace
        Tracer test_particles
        char *name

    ctypedef struct Params:
        double tpast
        double tfuture
        double dt0
        double dtout
        int ngals
        char *outputdir
        int snapshot
        int write_tracers

    ctypedef struct Snapshot:
        char *name
        int stripped
        double pos[3]
        double vel[3]
        double t


def run(dict input_parameters):
    cdef str param
    cdef dict galaxy
    cdef int i, n, j
    TRACERS = False
    cdef Params parameters
    cdef Snapshot **output_snapshots

    if 'galaxies' not in input_parameters.keys():
        raise ValueError("Must define galaxies to integrate.")

    cdef int ngals = len(input_parameters['galaxies'])
    cdef Gal *gal = <Gal *> malloc(ngals*sizeof(Gal))

    try:
        for n, (gal_name, galaxy) in enumerate(input_parameters['galaxies'].iteritems()):
            gal[n].name = gal_name
            gal[n].mhalo = galaxy['mass']
            gal[n].mtidal = galaxy['mass']
            gal[n].r_halo = galaxy['rad']
            gal[n].gamma = galaxy['gamma']
            gal[n].a2_LMJ = galaxy['a2']
            gal[n].b2_LMJ = galaxy['b2']
            gal[n].M2_LMJ = galaxy['m2']
            gal[n].M1_LMJ = galaxy['m1']
            gal[n].b1_LMJ = galaxy['b1']
            gal[n].dyn_fric = galaxy['dynamical_friction']
            gal[n].inplace = galaxy['inplace']
            if galaxy['dynamical_friction'] == 1:
                gal[n].dyn_C_eq = galaxy['dyn_C_eq']
                gal[n].dyn_L_eq = galaxy['dyn_L_eq']
                gal[n].dyn_alpha_eq = galaxy['dyn_alpha_eq']
                gal[n].dyn_C_uneq = galaxy['dyn_C_uneq']
                gal[n].dyn_L_uneq = galaxy['dyn_L_uneq']
                gal[n].dyn_alpha_uneq = galaxy['dyn_alpha_uneq']

            gal[n].tidal_trunc = galaxy['tidal_truncation']
            gal[n].rt = np.nan
            gal[n].stripped = 0
            for i in range(3):
                gal[n].pos[i] = galaxy['pos'][i]
                gal[n].vel[i] = galaxy['vel'][i]
                gal[n].post[i] = galaxy['pos'][i]
                gal[n].velt[i] = galaxy['vel'][i]

            # TESTING tracers
            if galaxy['tracers'] == 1:
                TRACERS = True
                gal[n].test_particles.nparticles = 1000
                gal[n].test_particles.pos = <double *>malloc(sizeof(double) * 3* 1000)
                gal[n].test_particles.vel = <double *>malloc(sizeof(double) * 3* 1000)
            else:
                gal[n].test_particles.nparticles = 0

    except KeyError, e:
        print('Missing parameter from galaxy %s' % gal[n].name)
        raise KeyError, e
    # Read parameters
    parameters.tpast = input_parameters["tpast"]
    parameters.dtout = input_parameters["dtout"]
    parameters.tfuture = input_parameters["tfuture"]
    parameters.dt0 = input_parameters['dt0']
    parameters.ngals = ngals
    parameters.snapshot = input_parameters['save_snapshot']
    parameters.outputdir = input_parameters["outputdir"]
    # only save test/tracer particles if there are any!
    if TRACERS:
        parameters.write_tracers = input_parameters["save_tracers"]
    else:
        parameters.write_tracers = 0

    cdef int nsnaps = 0
    if (parameters.tpast < 0.0):
        if (parameters.tfuture <= 0.0):
            nsnaps = int(parameters.tpast/(-1.0*parameters.dtout)+1)
        else:  # otherwise just save a snapshot at the end of the backward integration
            nsnaps = int(parameters.tpast/(-1.0*parameters.dtout) +
                         parameters.tfuture/parameters.dtout + 1)

        output_snapshots = <Snapshot **>malloc(sizeof(Snapshot *) * nsnaps)
        for i in range(nsnaps):
            output_snapshots[i] = <Snapshot *>malloc(sizeof(Snapshot) * ngals)


    if (parameters.tfuture > 0.0):
        if (nsnaps == 0):
            nsnaps = int(parameters.tfuture/parameters.dtout + 1)
            output_snapshots = <Snapshot **>malloc(sizeof(Snapshot *) * nsnaps)
            for i in range(nsnaps):
                output_snapshots[i] = <Snapshot *>malloc(sizeof(Snapshot) * ngals)


    cdef np.ndarray[double, ndim=1, mode="c"] output_pos = np.zeros(3*ngals)
    cdef np.ndarray[double, ndim=1, mode="c"] output_vel = np.zeros(3*ngals)

    err = orbit(ngals, parameters, gal, output_snapshots)

    if (err > 0):
        raise RuntimeError("Problem in integration. Error code: %d" % err)

    #err = orbit(mode, ngals, input_parameters, &output_pos[0], &output_vel[0])
    # Had a weird error with the string representation of the output
    try:
        _ = output_pos.__str__()
    except:
        pass

    output = []
    for i in xrange(nsnaps):
        output.append([])
        for j in xrange(ngals):
            s = output_snapshots[i][j]
            output[i].append({'name': s.name,
                              'stripped': s.stripped,
                              'pos': [p for p in s.pos],
                              'vel': [v for v in s.vel], 't': s.t})

    for i in range(nsnaps):
        free(output_snapshots[i])
    free(output_snapshots)
    return output


def likelihood(int ngals,
               np.ndarray[double, ndim=2, mode="c"] model_position,
               np.ndarray[double, ndim=2, mode="c"] model_velocity,
               np.ndarray[double, ndim=2, mode="c"] data_position,
               np.ndarray[double, ndim=2, mode="c"] data_velocity,
               double error_pos,
               double error_vel):

    cdef np.ndarray[double, ndim=1, mode="c"] dist2 = np.zeros(ngals)

    for i in xrange(ngals):
        # squared distance formula
        dist2[i] = np.sum((model_position[i, :]-data_position[i, :])**2)

    cdef np.ndarray[double, ndim=1, mode="c"] veldist2 = np.zeros(ngals)

    for i in xrange(ngals):
        # squared distance formula for velocities
        veldist2[i] = np.sum((model_velocity[i, :]-data_velocity[i, :])**2)


    cdef double ln_likelihood_gal_pos = (ngals*np.log(1.0/np.sqrt(2.0*np.pi*error_pos**2)) +
                                         np.sum(-0.5*dist2/error_pos**2))

    cdef double ln_likelihood_gal_vel = (ngals*np.log(1.0/np.sqrt(2.0*np.pi*error_vel**2)) +
                                         np.sum(-0.5*veldist2/error_vel**2))

    cdef double ln_likelihood = (ln_likelihood_gal_pos +
                                 ln_likelihood_gal_vel)

    return ln_likelihood


def likelihood2(np.ndarray[double, ndim=2, mode="c"] model_pos,
                np.ndarray[double, ndim=2, mode="c"] data_pos,
                double error_pos):
    cdef int i

    cdef double Small = 1e-5

    cdef int n_model = model_pos.shape[0]
    cdef int n_data = data_pos.shape[0]

    cdef double ln_likelihood_pos = 0.0
    #cdef double ln_likelihood_vel = 0.0

    for i in xrange(n_data):
        # Don't add in any NaN values
        if (any(model_pos[i, :] != model_pos[i, :])) or (any(data_pos[i, :] != data_pos[i, :])):
            continue
        ln_likelihood_pos += np.log(1.0/n_model*np.sum(np.exp(-0.5*((model_pos[:, 0]-data_pos[i, 0])**2 +
                                                                    (model_pos[:, 1]-data_pos[i, 1])**2 +
                                                                    (model_pos[:, 2]-data_pos[i, 2])**2)/
                                                                    error_pos**2)+Small))

         #   ln_likelihood_vel += np.log(1.0/n_model*np.sum(np.exp(-0.5*((model_vel[j, 0]-data_vel[i, 0])**2 +
         #                                                               (model_vel[j, 1]-data_vel[i, 1])**2 +
         #                                                               (model_vel[j, 2]-data_vel[i, 2])**2)/
         #                                                               error_vel**2)+Small))

    cdef double ln_likelihood = (ln_likelihood_pos) #+
                                # ln_likelihood_vel)
    return ln_likelihood


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
        results = run(input_parameters)
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

    ln_likelihood = likelihood2(model_pos,
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
        results = run(input_parameters)
    except RuntimeError, e:
        print({"Runtime Error: {:s}".format(e)})
        return None
    except KeyError, e:
        print("Problem in initial conditions. Missing parameter {:s}".format(e))
        return None

    cdef int ngals = len(input_parameters['galaxies'])
    cdef double ln_likelihood

    model_pos = np.array([results[-1][i]['pos'][j]
                          for i in xrange(ngals)
                          for j in xrange(3)])

    model_vel = np.array([results[-1][i]['vel'][j]
                          for i in xrange(ngals)
                          for j in xrange(3)])


    data_pos = np.empty((ngals, 3))
    data_vel = np.empty((ngals, 3))
    for i, (galname, gal) in enumerate(input_parameters['galaxies'].iteritems()):
        data_pos[i, :] = gal['pos']
        data_vel[i, :] = gal['vel']

    ln_likelihood = likelihood(ngals,
                               model_pos,
                               model_vel,
                               input_parameters['pos'],
                               input_parameters['vel'],
                               input_parameters['pos_err'],
                               input_parameters['vel_err'])
    return ln_likelihood


def orbit_statistics(dict input_parameters,
                     str gal_name,
                     str gal_name2):

    cdef list results

    try:
        results = run(input_parameters)
    except RuntimeError, e:
        print("Runtime Error: {:s}".format(e))
        return None
    except KeyError, e:
        print("Problem in initial conditions. Missing parameter {:s}".format(e))
        return None

    return test_orbit_statistics(results, len(input_parameters['galaxies']), gal_name, gal_name2)


def test_stream(dict input_parameters,
                np.ndarray[double, ndim=2, mode="c"] input_position):

    cdef list results

    try:
        results = run(input_parameters)
    except RuntimeError, e:
        print("Runtime Error: {:s}".format(e))
        return None
    except KeyError, e:
        print("Problem in initial conditions. Missing parameter {:s}".format(e))
        return None

    mindist, maxdist, apos, peris, stripped = test_orbit_statistics(results,
                                                                    len(input_parameters['galaxies']),
                                                                    "MW", "LMC")

    if (peris != 2):
        return -1e+5

    cdef double ln_likelihood = 0.0
    cdef int ngals = len(input_parameters['galaxies'])

    model_pos = np.array([results[-1][i]['pos'][j]
                          for i in xrange(ngals)
                          for j in xrange(3)])

    model_vel = np.array([results[-1][i]['vel'][j]
                          for i in xrange(ngals)
                          for j in xrange(3)])

    data_pos = np.empty((ngals, 3))
    data_vel = np.empty((ngals, 3))
    for i, (galname, gal) in enumerate(input_parameters['galaxies'].iteritems()):
        data_pos[i, :] = gal['pos']
        data_vel[i, :] = gal['vel']

    ln_likelihood += likelihood(ngals,
                                model_pos,
                                model_vel,
                                input_parameters['pos'],
                                input_parameters['vel'],
                                input_parameters['pos_err'],
                                input_parameters['vel_err'])

    #model_pos = np.zeros((1000, 3))
    #for i in xrange(1000):
    #    for j in xrange(3):
    #        model_pos[i, j] =
    # convert x, y, z to l,b

    #ln_likelihood += likelihood2(model_pos,
    #                             input_position,
    #                             input_parameters['pos_err'])

    return ln_likelihood


def test_orbit_statistics(list results,
                          int ngals,
                          str gal_name,
                          str gal_name2):

    cdef double ln_likelihood
    cdef int nsnaps, g, g2, s
    cdef double d

    nsnaps = len(results)
    g = np.where([results[0][i]['name'] == gal_name for i in xrange(ngals)])[0]
    g2 = np.where([results[0][i]['name'] == gal_name2 for i in xrange(ngals)])[0]

    dist = [np.sqrt((results[i][g]['pos'][0]-results[i][g2]['pos'][0])**2 +
                    (results[i][g]['pos'][1]-results[i][g2]['pos'][1])**2 +
                    (results[i][g]['pos'][2]-results[i][g2]['pos'][2])**2) for i in xrange(nsnaps)]

    stripped = [results[i][g2]['stripped'] for i in xrange(nsnaps)]
    old_dist = None
    direction = None

    cdef int apocenters = 0
    cdef int pericenters = 0

    for s, d in zip(stripped, dist):
        if s == 1:
            break  # if the galaxy becomes stirpped, stop
        if old_dist is not None:
            if direction is not None:
            # already have a direction set
                if ((direction == 1) and (d >= old_dist)):
                # moving in and got further away
                    direction = 0
                    pericenters += 1
                elif ((direction == 0) and (d <= old_dist)):
                # moving out and got closer
                    direction = 1
                    apocenters += 1
            else:
            # set direction
                if (d < old_dist):
                    direction = 1  # inward
                else:
                    direction = 0  # outward
        old_dist = d

    return np.min(dist), np.max(dist), apocenters, pericenters, any(stripped)

