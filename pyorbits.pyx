import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc, free


cdef extern from "orbit.c":
    int orbit(int int_mode,
              int ngals,
              dict input_parameters,
              double* output_pos,
              double* output_vel)

    int orbit_new(int int_mode,
              int ngals,
              Params parameters,
              Gal *gal,
              double* output_pos,
              double* output_vel)


cdef extern from *:
    struct Gal:
        double pos[3]
        double vel[3]
        double post[3]
        double velt[3]
        int ID
        double mhalo
        double r_halo
        double gamma
        double a2_LMJ
        double b2_LMJ
        double M2_LMJ
        double M1_LMJ
        double b1_LMJ
        double c_halo
        int halo_type


    struct Params:
        double b1_LMJ
        double M1_LMJ
        double a2_LMJ
        double b2_LMJ
        double M2_LMJ
        double Mhalo
        double q_halo
        double r_halo
        double gamma
        double c_halo
        double halo_type
        double tpast
        double tfuture
        double dt0
        double dtout
        int ngals
        char *outputdir


def run(int mode, dict input_parameters):
    cdef str param
    cdef list galaxy_params = ['rad_gal',
                               'mass_gal',
                               'gamma_gal',
                               'a2_gal',
                               'b2_gal',
                               'm2_gal',
                               'pos',
                               'vel']
    for param in galaxy_params:
        if param not in input_parameters.keys():
            print '{} not in input parameters'.format(param)
            return None
        elif input_parameters[param].shape[0] != input_parameters['rad_gal'].shape[0]:
            print '{} not the correct size'.format(param)
            return None


    cdef int ngals = input_parameters['rad_gal'].shape[0]
    cdef Gal *gal = <Gal *> malloc(ngals*sizeof(Gal))
    cdef Params parameters

    # Read parameters
    parameters.b1_LMJ = input_parameters["b1_LMJ"]
    parameters.M1_LMJ = input_parameters["M1_LMJ"]
    parameters.a2_LMJ = input_parameters["a2_LMJ"]
    parameters.b2_LMJ = input_parameters["b2_LMJ"]
    parameters.M2_LMJ = input_parameters["M2_LMJ"]
    parameters.halo_type = input_parameters["halo_type"]
    parameters.Mhalo = input_parameters["Mhalo"]
    parameters.q_halo = input_parameters["q_halo"]
    parameters.r_halo = input_parameters["r_halo"]
    parameters.tpast = input_parameters["tpast"]
    parameters.dtout = input_parameters["dtout"]
    parameters.tfuture = input_parameters["tfuture"]
    parameters.dt0 = input_parameters['dt0']
    parameters.ngals = ngals

    if (parameters.halo_type == 1): # NFW
        parameters.c_halo = input_parameters["c_halo"]
    else: # Dehnen
        parameters.gamma = input_parameters["gamma_halo"]

    parameters.outputdir = input_parameters["outputdir"]

    for n in range(ngals):
        gal[n].mhalo = input_parameters['mass_gal'][n]
        gal[n].r_halo = input_parameters['rad_gal'][n]
        gal[n].gamma = input_parameters['gamma_gal'][n]
        gal[n].c_halo = input_parameters['c_gal'][n]
        gal[n].a2_LMJ = input_parameters['a2_gal'][n]
        gal[n].b2_LMJ = input_parameters['b2_gal'][n]
        gal[n].M2_LMJ = input_parameters['m2_gal'][n]
        gal[n].M1_LMJ = input_parameters['m1_gal'][n]
        gal[n].b1_LMJ = input_parameters['b1_gal'][n]
        gal[n].halo_type = input_parameters['gal_types'][n]
        for i in range(3):
            gal[n].pos[i] = input_parameters['pos'][n, i]
            gal[n].vel[i] = input_parameters['vel'][n, i]

    cdef np.ndarray[double, ndim=1, mode="c"] output_pos = np.zeros(3*ngals)
    cdef np.ndarray[double, ndim=1, mode="c"] output_vel = np.zeros(3*ngals)

    err = orbit_new(mode, ngals, parameters, gal, &output_pos[0], &output_vel[0])
    #err = orbit(mode, ngals, input_parameters, &output_pos[0], &output_vel[0])
    try:
        _ = output_pos.__str__()
    except:
        pass
    return {'pos': output_pos, 'vel': output_vel}


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


def test_orbit(dict input_parameters):
    cdef dict results = run(1, input_parameters)
    cdef int ngals = input_parameters['rad_gal'].shape[0]
    cdef double ln_likelihood
    if results is not None:
        ln_likelihood = likelihood(ngals,
                                   results['pos'],
                                   results['vel'],
                                   input_parameters['pos'],
                                   input_parameters['vel'],
                                   input_parameters['pos_err'],
                                   input_parameters['vel_err'])
        return ln_likelihood
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
