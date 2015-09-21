import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc


cdef extern from "orbit.c":
    int orbit(int ngals,
              Params parameters,
              Gal *gal,
              Snapshot **output_snapshots)


cdef extern from *:
    struct Gal:
        double pos[3]
        double vel[3]
        double post[3]
        double velt[3]
        int ID
        double mhalo
        double minit
        double r_halo
        double gamma
        double a2_LMJ
        double b2_LMJ
        double M2_LMJ
        double M1_LMJ
        double b1_LMJ
        double c_halo
        int dyn_fric
        double dyn_C_eq
        double dyn_L_eq
        double dyn_alpha_eq
        double dyn_C_uneq
        double dyn_L_uneq
        double dyn_alpha_uneq
        int tidal_trunc
        double rt
        int halo_type
        int inplace
        char *name


    struct Params:
        double tpast
        double tfuture
        double dt0
        double dtout
        int ngals
        char *outputdir

    struct Snapshot:
        char *name
        double pos[3]
        double vel[3]
        double t


def run(dict input_parameters):
    cdef str param
    cdef dict galaxy

    if 'galaxies' not in input_parameters.keys():
        raise ValueError("Must define galaxies to integrate.")
        return None

    cdef Snapshot **output_snapshots

    cdef int ngals = len(input_parameters['galaxies'])
    cdef Gal *gal = <Gal *> malloc(ngals*sizeof(Gal))
    for n, (gal_name, galaxy) in enumerate(input_parameters['galaxies'].iteritems()):
        gal[n].name = gal_name
        gal[n].mhalo = galaxy['mass']
        gal[n].minit = galaxy['mass']
        gal[n].r_halo = galaxy['rad']
        gal[n].gamma = galaxy['gamma']
        gal[n].c_halo = galaxy['c']
        gal[n].a2_LMJ = galaxy['a2']
        gal[n].b2_LMJ = galaxy['b2']
        gal[n].M2_LMJ = galaxy['m2']
        gal[n].M1_LMJ = galaxy['m1']
        gal[n].b1_LMJ = galaxy['b1']
        gal[n].halo_type = galaxy['type']
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
        for i in range(3):
            gal[n].pos[i] = galaxy['pos'][i]
            gal[n].vel[i] = galaxy['vel'][i]


    cdef Params parameters

    # Read parameters
    parameters.tpast = input_parameters["tpast"]
    parameters.dtout = input_parameters["dtout"]
    parameters.tfuture = input_parameters["tfuture"]
    parameters.dt0 = input_parameters['dt0']
    parameters.ngals = ngals

    parameters.outputdir = input_parameters["outputdir"]
    cdef int nsnaps = 0
    if (parameters.tpast < 0.0):
        if (parameters.tfuture <= 0.0):
            nsnaps = int(parameters.tpast/(-1.0*parameters.dtout)+1)
        else:  # otherwise just save a snapshot at the end of the backward integration
            nsnaps = int(parameters.tfuture/parameters.dtout+1)

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
        return None

    #err = orbit(mode, ngals, input_parameters, &output_pos[0], &output_vel[0])
    # Had a weird error with the string representation of the output
    try:
        _ = output_pos.__str__()
    except:
        pass
    #free(gal)
    output = []
    for i in xrange(nsnaps):
        output.append([])
        for j in xrange(ngals):
            s = output_snapshots[i][j]
            output[i].append({'name': s.name, 'pos': [p for p in s.pos],
                              'vel': [v for v in s.vel], 't': s.t})

    return output


def likelihood(int ngals,
               np.ndarray[double, ndim=1, mode="c"] model_position,
               np.ndarray[double, ndim=1, mode="c"] model_velocity,
               np.ndarray[double, ndim=2, mode="c"] data_position,
               np.ndarray[double, ndim=2, mode="c"] data_velocity,
               double error_pos,
               double error_vel):

    cdef np.ndarray[double, ndim=1, mode="c"] dist2 = np.zeros(ngals)

    for i in xrange(ngals):
        # squared distance formula
        dist2[i] = np.sum((model_position[i*3:(i+1)*3]-data_position[i, :])**2)

    cdef np.ndarray[double, ndim=1, mode="c"] veldist2 = np.zeros(ngals)

    for i in xrange(ngals):
        # squared distance formula for velocities
        veldist2[i] = np.sum((model_velocity[i*3:(i+1)*3]-data_velocity[i, :])**2)


    cdef double ln_likelihood_gal_pos = (ngals*np.log(1.0/np.sqrt(2.0*np.pi*error_pos**2)) +
                                         np.sum(-0.5*dist2/error_pos**2))

    cdef double ln_likelihood_gal_vel = (ngals*np.log(1.0/np.sqrt(2.0*np.pi*error_vel**2)) +
                                         np.sum(-0.5*veldist2/error_vel**2))

    cdef double ln_likelihood = (ln_likelihood_gal_pos +
                                 ln_likelihood_gal_vel)

    return ln_likelihood


def test_location(dict input_parameters):
    """
    Test against the current location of galaxies.
    Args:
        dict input_parameters
        See above
    """
    cdef list results = run(input_parameters)
    cdef int ngals = len(input_parameters['galaxies'])
    cdef double ln_likelihood
    if results is not None:
        model_pos = np.array([results[-1][i]['pos'][j]
                              for i in xrange(ngals)
                              for j in xrange(3)])

        model_vel = np.array([results[-1][i]['vel'][j]
                              for i in xrange(ngals)
                              for j in xrange(3)])

        ln_likelihood = likelihood(ngals,
                                   model_pos,
                                   model_vel,
                                   input_parameters['pos'],
                                   input_parameters['vel'],
                                   input_parameters['pos_err'],
                                   input_parameters['vel_err'])
        return ln_likelihood
    else:
        raise RuntimeError('Incorrect parameters or problem with run')



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
    cdef list results = run(input_parameters)
    cdef double ln_likelihood
    cdef int ngals = len(input_parameters['galaxies'])
    cdef int nsnaps, g
    if results is not None:
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
    else:
        raise RuntimeError('Incorrect parameters or problem with run')






def likelihood2(np.ndarray[double, ndim=2, mode="c"] model_pos,
                np.ndarray[double, ndim=2, mode="c"] data_pos,
                double error_pos):
    cdef int i, j

    cdef double Small = 1e-5

    cdef int n_model = model_pos.shape[0]
    cdef int n_data = data_pos.shape[0]

    cdef double ln_likelihood_pos = 0.0
   # cdef double ln_likelihood_vel = 0.0

    for i in xrange(n_model):
        # Don't add in any NaN values
        if (any(model_pos[j, :] != model_pos[j, :])) or (any(data_pos[i, :] != data_pos[i, :])):
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


def orbit_statistics(dict input_parameters,
                     str gal_name,
                     str gal_name2,
                     str output_file):

    cdef list results = run(input_parameters)
    cdef double ln_likelihood
    cdef int ngals = len(input_parameters['galaxies'])
    cdef int nsnaps, g
    if results is not None:
        nsnaps = len(results)
        g = np.where([results[0][i]['name'] == gal_name for i in xrange(ngals)])[0]
        g2 = np.where([results[0][i]['name'] == gal_name2 for i in xrange(ngals)])[0]

        dist = [np.sqrt(np.sum((results[i][g]['pos']-results[i][g2]['pos'])**2)) for i in xrange(nsnaps)]

        old_dist = None
        direction = None

        apocenters = 0
        pericenters = 0

        for d in dist:
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

        return np.min(dist), np.max(dist), apocenters, pericenters
    else:
        raise RuntimeError('Incorrect parameters or problem with run')


# On hold until we add tracer particles
#    cdef int i
#    cdef double chi2 = 0.0
#
#    cdef double sigma_x2 = sigma_x**2
#    cdef double sigma_v2 = sigma_v**2
#    cdef double sigma_vx2 = sigma_vx**2
#    cdef double sigma_mu2 = sigma_mu**2
#    cdef int nr_OD = n_OD.shape[0]
#    cdef int nr_VR = n_VR.shape[0]
#    cdef double sigma_OD[nr_OD]
#    cdef double normvel[nr_VR]
#
#    cdef double dl, db, dl2, db2, dvx, dvx2
#    dl = np.sqrt(0.125*0.125+sigma_x2)
#    db = sqrt(0.125*0.125+sigma_x2)
#    dvx = sqrt(0.125*0.125+sigma_vx2)
#    dl2 = dl**2  # squared values
#    db2 = db**2
#    dvx2 = dvx**2
#
#    sigma_OD[:] =  sigma_x2+n_OD[:, 4]**2/(8.0*np.log(2.0))  # //squared values of sigma_x plus sigma_obs
#
#    for i in xrange(N):
#        normvel[:] += np.exp(-0.5*((star[i, 0]-n_VR[:, 4])**2/dvx2 +
#                                    (star[i, 1]-n_VR[:, 3])**2/dvx2 +
#                                    (star[i, 2]-n_VR[:, 0])**2/(n_VR[:, 1]**2+sigma_v2) +
#                                    (star[i, 3]-n_VR[:, 6])**2/(n_VR[:, 8]**2+sigma_mu2) +
#                                    (star[i, 4]-n_VR[:, 7])**2/(n_VR[:, 9]**2+sigma_mu2)))
#
#        #  overdensities
#        n_OD[:, 2] +=  np.exp(-0.5*((star[i, 0]-n_OD[:, 0])**2/sigma_OD[:] +
#                                     (star[i, 1]-n_OD[:, 1])**2/sigma_OD[:]))
#
#
#    # normalization velocities
#    for i in xrange(nr_VR):
#        normvel[5] *= (1.0/(1.0*N*dvx2*sqrt(n_VR[i, 1]**2+sigma_v2)
#                            *np.sqrt(n_VR[i, 8]**2+sigma_mu2)
#                            *np.sqrt(n_VR[i, 9]**2sigma_mu2)))
#
#
#    #normalization overdensities
#    n_OD[:, 2] *=  (1.0/(1.0*N*sigma_OD[:]))
#
#
#    #construction of final likelihood value
#
#    chi2 += np.sum(n_OD[:, 3]*np.log(n_OD[:, 2]+SMALL))
#    chi2 += np.sum(np.log(normvel[:]+SMALL))
#
#    return chi2
