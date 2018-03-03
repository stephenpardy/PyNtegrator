import pyorbits
import numpy as np
#

mw_pos = np.array([0.0, 0.0, 0.0])
mw_vel = np.array([0.0, 0.0, 0.0])

lmc_pos = np.array([-1.0, -41.0, -28.0])
#lmc_vel = np.array([-57.0, -226.0, 221.0])  # New solar
lmc_vel = [-78, -236, 196]

smc_vel = np.array([19.0, -153.0, 153.0])  # New Solar
smc_pos = np.array([ 14.88893139, -38.08632213, -44.20286079])

#LMCMASS = 9.9
LMCMASS = 11.63

gadget_LMC = {"rad": 14.0013,#7.93195, #14.0013#21.4,
          "mass": LMCMASS,
          "gamma": 1.0,
          "a2": 1.43,
          "m2": 0.02*LMCMASS,
          "b2": 0.312,
          "b1": 1.0,
          "m1": 0.0,
          "dynamical_friction": 0,
          "tidal_truncation": 0,
          "inplace": 0,
          "tracers": 1,
          "pos": lmc_pos,
          "vel": lmc_vel
        }


gadget_MW = {"rad": 29.8,
          #"mass": 120.0,
          "mass": 115.6,
          "gamma": 1.0,
          "a2": 3.5,
          "m2": 4.8,
          "b1": 0.7,
          "b2": 0.53,
          "m1": 1.2,
          "dynamical_friction": 1,
          "dyn_C_uneq": 1.6*3.0/29.8, #1.22
          "dyn_L_uneq": 0.0,
          "dyn_alpha_uneq": 1.5,#uneq_alpha,
          "dyn_C_eq": 0.17,
          "dyn_L_eq": 0.02,
          "dyn_alpha_eq": 1.0,
          "tidal_truncation": 0,
          "inplace": 0,
          "tracers": 0,
          "pos": np.array([0.0, 0.0, 0.0]),
          "vel": np.array([0.0, 0.0, 0.0])
        }


gadget_SMC = {"rad": 2.99579,
              "mass": 0.3,
              "gamma": 1.0,
              "a2": 0.279,
              "m2": 0.015,
              "b2": 0.062,
              "b1": 1.0,
              "m1": 0.0,
              "dynamical_friction": 0,
              "tidal_truncation": 0,
              "inplace": 0,
              "tracers": 1,
              "pos": smc_pos,
              "vel": smc_vel
              }


params_stream = {"galaxies":{"MW": gadget_MW,
                             "LMC": gadget_LMC,
                             "SMC": gadget_SMC
                             },
                  "pos_err": 100.0,  # 0.1 kpc randomly choosen for now
                  "vel_err": 20.0,  # km/s rough avg. from 3rd epoch data
                  "tpast": -8.0,
                  "tfuture": 0.01,
                  "dt0": 1e-3,
                  "outputdir": "./output_stream/",
                  "dtout":0.05,
                  "save_snapshot": 1,
                  "save_tracers": 1,
                  "output_tracers": 1}

stream_data = np.loadtxt('stream.pos.csv', delimiter=',')
#print _orbit.orbit(1,0,2,params)
print pyorbits.test_stream(params_stream, stream_data, True)

