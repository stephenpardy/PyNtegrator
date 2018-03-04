#encoding: utf-8
#cython: profile=True

import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc, free
import copy


cdef extern from 'orbit.c':
    int orbit(Params parameters,
              Gal *gal,
              Snapshot **output_snapshots)

cdef extern from *:

    struct Tracer:
        int nparticles
        double *pos
        double *vel

    struct Gal:
        double pos[3]
        double vel[3]
        double post[3]
        double velt[3]
        double mhalo
        double minit
        double r_halo
        double gamma
        int halo_type
        double a_disk
        double b_disk
        double M_disk
        double M_bulge
        double b_bulge
        int dyn_fric
        int mass_growth
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

    struct Params:
        double tpast
        double tfuture
        double dt0
        double dtout
        int ngals
        char *outputdir
        int snapshot
        int write_tracers
        int variabletimesteps

    struct Snapshot:
        char *name
        int stripped
        double pos[3]
        double vel[3]
        double t


def run(dict input_parameters,
        tracers=None):
    """
    Main run function
    args:
        dict input_parameters:
            A dictionary of input parameters. See below for structure.
        dict tracers:
            Either None or a dictionary containing the tracer particles for each galaxy. 
                Should have structure:
                {'<GALNAME>': (np.array(<SIZEOF N, 3>), np.array(<SIZEOF N, 3>))}
                The first element of the array is for the positions and the second element is for the velocities

    """
    cdef str param
    cdef dict galaxy
    cdef int i, n, j
    TRACERS = False
    cdef Params parameters
    cdef Snapshot **output_snapshots

    if 'galaxies' not in input_parameters.keys():
        raise ValueError("Must define galaxies to integrate.")

    #Initialize gal array
    cdef int ngals = len(input_parameters['galaxies'])
    cdef Gal *gal = <Gal *> malloc(ngals*sizeof(Gal))
    #Read galaxy parameters
    try:
        for n, (name, galaxy) in enumerate(input_parameters['galaxies'].items()):
            #temporary backwards compatibility check
            if 'halo_type' not in galaxy.keys():
                galaxy['halo_type'] = 0  # default is Hernquist
            name = name.encode('utf-8')
            gal[n].name = name
            gal[n].mhalo = galaxy['mass']
            gal[n].minit = galaxy['mass']
            gal[n].r_halo = galaxy['rad']
            gal[n].gamma = galaxy['gamma']
            gal[n].halo_type = galaxy['halo_type']
            gal[n].a_disk = galaxy['a2']
            gal[n].b_disk = galaxy['b2']
            gal[n].M_disk = galaxy['m2']
            gal[n].M_bulge = galaxy['m1']
            gal[n].b_bulge = galaxy['b1']
            gal[n].dyn_fric = galaxy['dynamical_friction']
            gal[n].mass_growth = galaxy['mass_growth']
            gal[n].inplace = galaxy['inplace']
            if galaxy['dynamical_friction'] == 1:
                gal[n].dyn_C_eq = galaxy['dyn_C_eq']
                gal[n].dyn_L_eq = galaxy['dyn_L_eq']
                gal[n].dyn_alpha_eq = galaxy['dyn_alpha_eq']
                gal[n].dyn_C_uneq = galaxy['dyn_C_uneq']
                gal[n].dyn_L_uneq = galaxy['dyn_L_uneq']
                gal[n].dyn_alpha_uneq = galaxy['dyn_alpha_uneq']

            gal[n].tidal_trunc = galaxy['tidal_truncation']
            gal[n].rt = galaxy['tidal_radius']

            gal[n].stripped = 0
            for i in range(3):
                gal[n].pos[i] = galaxy['pos'][i]
                gal[n].vel[i] = galaxy['vel'][i]
                gal[n].post[i] = galaxy['pos'][i]
                gal[n].velt[i] = galaxy['vel'][i]

            # TESTING tracers
            if galaxy['tracers'] == 1:
                TRACERS = True
                if tracers is not None:
                    try:
                        pos, vel = tracers[name]
                        ntracers = pos.shape[0]
                        gal[n].test_particles.nparticles = ntracers
                        gal[n].test_particles.pos = <double *>malloc(sizeof(double) * 3 * ntracers)
                        gal[n].test_particles.vel = <double *>malloc(sizeof(double) * 3 * ntracers)
                        for i in xrange(ntracers):
                            for j in xrange(3):
                                gal[n].test_particles.pos[i*3+j] = pos[i, j]
                                gal[n].test_particles.vel[i*3+j] = vel[i, j]

                    except KeyError, e:
                        free(gal)
                        print('No tracers found for galaxy {:s}'.format(name))
                        raise KeyError, e
                else:
                    free(gal)
                    print('No tracers found!')
                    raise RuntimeError
            else:
                gal[n].test_particles.nparticles = 0

    except KeyError, e:
        free(gal)
        print('Missing parameter from galaxy %s' % gal[n].name)
        raise KeyError, e
    # Read integration parameters
    parameters.tpast = input_parameters["tpast"]
    parameters.dtout = input_parameters["dtout"]
    parameters.tfuture = input_parameters["tfuture"]
    parameters.dt0 = input_parameters['dt0']
    parameters.ngals = ngals
    parameters.snapshot = input_parameters['save_snapshot']
    outputdir = input_parameters["outputdir"].encode('utf-8')
    parameters.outputdir = outputdir
    parameters.variabletimesteps = input_parameters["variable_timesteps"]
    # only save test/tracer particles if there are any!
    if TRACERS:
        parameters.write_tracers = input_parameters["save_tracers"]
        output_tracers = input_parameters["output_tracers"]
    else:
        parameters.write_tracers = 0
        output_tracers = 0

    #Initialize output arrays based on whether we want to integrate forward to backward
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

    elif (parameters.tfuture > 0.0):
        if (nsnaps == 0):
            nsnaps = int(parameters.tfuture/parameters.dtout + 1)
            output_snapshots = <Snapshot **>malloc(sizeof(Snapshot *) * nsnaps)
            for i in range(nsnaps):
                output_snapshots[i] = <Snapshot *>malloc(sizeof(Snapshot) * ngals)
    else:
        free(gal)
        raise RuntimeError("You must either set tpast < 0 or tfuture > 0")

    print('running, nsnaps: {:d}'.format(nsnaps))

    err = orbit(parameters, gal, output_snapshots)

    if (err > 0):

        raise RuntimeError("Problem in integration. Error code: %d" % err)

    #Generate an array of galaxy snapshots
    output = []
    for i in xrange(nsnaps):
        output.append([])
        for j in xrange(ngals):
            s = output_snapshots[i][j]
            output[i].append({'name': s.name,
                              'stripped': s.stripped,
                              'pos': [p for p in s.pos],
                              'vel': [v for v in s.vel], 't': s.t})

    # If we want to output final position of tracers then make that array
    tracers = None
    for i in xrange(ngals):
        npart = gal[i].test_particles.nparticles
        if (npart > 0):
            tracers_temp = np.empty((npart, 3))

            for j in xrange(npart):
                for k in  xrange(3):
                    tracers_temp[j, k] = gal[i].test_particles.pos[j*3+k]
            if (tracers is None):
                tracers = tracers_temp
            else:
                tracers = np.concatenate((tracers, tracers_temp), axis=0)
            free(gal[i].test_particles.pos)
            free(gal[i].test_particles.vel)


    #Free now that we are done with these
    for i in range(nsnaps):
        free(output_snapshots[i])
    free(output_snapshots)
    free(gal)

    if output_tracers:
        return output, tracers
    # Just output the galaxy positions
    else:
        return output


cpdef double likelihood(int ngals,
               np.ndarray[double, ndim=2, mode="c"] model_position,
               np.ndarray[double, ndim=2, mode="c"] model_velocity,
               np.ndarray[double, ndim=2, mode="c"] data_position,
               np.ndarray[double, ndim=2, mode="c"] data_velocity,
               double error_pos,
               double error_vel):

    cdef np.ndarray[double, ndim=1, mode="c"] dist2 = np.zeros(ngals)

    for i in xrange(ngals):
        # squared distance formula
        dist2[i] = np.sum((model_position[i, :] - data_position[i, :])**2)

    cdef np.ndarray[double, ndim=1, mode="c"] veldist2 = np.zeros(ngals)

    for i in xrange(ngals):
        # squared distance formula for velocities
        veldist2[i] = np.sum((model_velocity[i, :]-data_velocity[i, :])**2)


    cdef double ln_likelihood_gal_pos = (ngals*np.log(1.0/np.sqrt(2.0*np.pi*error_pos**2)) +
                                         np.sum(-0.5*dist2/error_pos**2))

    cdef double ln_likelihood_gal_vel = (ngals*np.log(1.0/np.sqrt(2.0*np.pi*error_vel**2)) +
                                         np.sum(-0.5*veldist2/error_vel**2))

    cdef double ln_likelihood = ln_likelihood_gal_pos + ln_likelihood_gal_vel

    return ln_likelihood


cpdef double likelihood2(np.ndarray[double, ndim=2, mode="c"] model_pos,
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

    cdef double ln_likelihood = ln_likelihood_pos #+ ln_likelihood_vel
    return ln_likelihood


def test_orbit_statistics(list results,
                          int ngals,
                          str gal_name,
                          str gal_name2,
                          float max_time=1000,
                          float min_time=-1000):

    cdef double ln_likelihood
    cdef int nsnaps, g, g2, s
    cdef double d

    nsnaps = len(results)
    g = np.where([results[0][i]['name'] == gal_name for i in xrange(ngals)])[0]
    g2 = np.where([results[0][i]['name'] == gal_name2 for i in xrange(ngals)])[0]
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

    idx = (times > min_time)*(times < max_time)

    for s, d in zip(stripped[idx], dist[idx]):
        if s == 1:
            break  # if the galaxy becomes stripped, stop
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

    # if we are currently moving in set that as a pericenter
    if d == 1:
        pericenters += 1
    # equivalently for moving out
    elif d == 0:
        apocenters += 1

    return np.min(dist), np.max(dist), apocenters, pericenters, any(stripped)

