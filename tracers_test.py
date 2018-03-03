import numpy as np
import pyorbits
import matplotlib.pyplot as plt


gal1 = {"rad": 29.8228,
          "mass": 120*0.96,
          "gamma": 1.0,
          "a2": 1.0,
          "m2": 0.0,
          "b2": 1.0,
          "b1": 2.82925,
          "m1": 120.0*0.04,
          "dynamical_friction": 0,
          "mass_growth": 0,
          "tidal_truncation": 1,
          "tidal_radius": np.nan,
          "inplace": 0,
          "tracers": 1,
          "pos": np.array([-47.1429, -16.6599, 0.0]),
          "vel": np.array([91.808, -66.916, 0])
        }

gal2 = {"rad": 29.8228,
          "mass": 120*0.96,
          "gamma": 1.0,
          "a2": 1.0,
          "m2": 0.0,
          "b2": 1.0,
          "b1": 2.82925,
          "m1": 120.0*0.04,
          "dynamical_friction": 0,
          "mass_growth": 0,
          "tidal_truncation": 0,
          "tidal_radius": np.nan,
          "inplace": 0,
          "tracers": 1,
          "pos": np.array([47.1429, 16.6599, -0]),
          "vel": np.array([-91.808, 66.916, -0])
        }

base_params = {"galaxies":{"GAL1": gal1,
                           "GAL2": gal2,
                            },
                  "pos_err": 1.0,
                  "vel_err": 1.0,
                  "tpast": 0.0,
                  "tfuture": 5.0,
                  "dt0": 1e-3,
                  "outputdir":"./",
                  "dtout":0.05,
                  "variable_timesteps": 0,
                  "save_snapshot": 0,
                  "save_tracers": 0,
                  "output_tracers": 1}


init_tracers_pos1 = np.empty((1000, 3), dtype='float_')
init_tracers_vel1 = np.empty((1000, 3), dtype='float_')

init_tracers_pos2 = np.empty((1000, 3), dtype='float_')
init_tracers_vel2 = np.empty((1000, 3), dtype='float_')

velx, vely, velz, posx, posy, posz = np.loadtxt('MW.csv',
                                                delimiter=',',
                                                unpack=True,
                                                skiprows=1)

init_tracers_pos1[:, 0] = posx + gal1['pos'][0]
init_tracers_pos1[:, 1] = posy + gal1['pos'][1]
init_tracers_pos1[:, 2] = posz

init_tracers_vel1[:, 0] = velx + gal1['vel'][0]
init_tracers_vel1[:, 1] = vely + gal1['vel'][1]
init_tracers_vel1[:, 2] = velz

init_tracers_pos2[:, 0] = posx + gal2['pos'][0]
init_tracers_pos2[:, 1] = posy + gal2['pos'][1]
init_tracers_pos2[:, 2] = posz

init_tracers_vel2[:, 0] = velx + gal2['vel'][0]
init_tracers_vel2[:, 1] = vely + gal2['vel'][1]
init_tracers_vel2[:, 2] = velz


output, tracers_final = pyorbits.run(base_params, {'GAL1': (init_tracers_pos1, init_tracers_vel1),
                                                   'GAL2': (init_tracers_pos2, init_tracers_vel2)})


#gal1:
x = np.array([o[0]['pos'][0] for o in output])
y = np.array([o[0]['pos'][1] for o in output])
plt.plot(x, y, color='black')

#gal2:
x = np.array([o[1]['pos'][0] for o in output])
y = np.array([o[1]['pos'][1] for o in output])
plt.plot(x, y, color='red')

plt.plot(tracers_final[:1000, 0], tracers_final[:1000, 1], 'k.')
plt.plot(tracers_final[1000:, 0], tracers_final[1000:, 1], 'r.')

plt.show()
