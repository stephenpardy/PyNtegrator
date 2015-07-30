import numpy as np
cimport numpy as np

cdef extern from "orbit.c":
    int orbit(int int_mode,
              int ngals,
              dict input_parameters,
              double* output_pos,
              double* output_vel)

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
    cdef np.ndarray[double, ndim=1, mode="c"] output_pos = np.zeros(3*ngals)
    cdef np.ndarray[double, ndim=1, mode="c"] output_vel = np.zeros(3*ngals)
    err = orbit(mode, ngals, input_parameters, &output_pos[0], &output_vel[0])
    print err
    try:
        _ = output_pos.__str__()
    except:
        pass
    return {'pos': np.array(output_pos), 'vel': output_vel}
