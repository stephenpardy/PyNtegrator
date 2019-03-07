import numpy as np
import pyorbits
import copy

gal1 = {"rad": 10,
          "mass": 100,
          "gamma": 1.0,
          "a2": 1.0,
          "m2": 0.0,
          "b2": 1.0,
          "b1": 1.0,
          "m1": 0.0,
          "dynamical_friction": 0,
          "mass_growth": 0,
          "tidal_truncation": 1,
          "tidal_radius": np.nan,
          "inplace": 0,
          "tracers": 0,
          "pos": np.array([10.0, 0.0, 0.0]),
          "vel": np.array([0.0, 10.0, 0.0])
        }

gal2 = {"rad": 10.0,
          "mass": 100.0,
          "gamma": 1.0,
          "a2": 1.0,
          "m2": 0.0,
          "b1": 1.0,
          "b2": 1.0,
          "m1": 0.0,
          "dynamical_friction": 0,
          "mass_growth": 0,
          "tidal_truncation": 0,
          "tidal_radius": np.nan,
          "inplace": 0,
          "tracers": 0,
          "pos": np.array([-10.0, 0.0, 0.0]),
          "vel": np.array([0.0, -10.0, 0.0])
        }

base_params = {"galaxies":{"GAL1": gal1,
                           "GAL2": gal2,
                            },
                  "pos_err": 1.0,
                  "vel_err": 1.0,
                  "tpast": -5.0,
                  "tfuture": 0.0,
                  "dt0": 1e-3,
                  "outputdir":"./",
                  "dtout":0.05,
                  "variable_timesteps": 0,
                  "save_snapshot": 0,
                  "save_tracers": 0,
                  "output_tracers": 0}


def test_integration():
    #params = copy.deepcopy(base_params)
    output = pyorbits.run(base_params)
    #check last snapshot for each galaxy
    #They should be equal and opposite
    assert output[-1][0]['pos'] == [35.88936905187977, -31.04336783591637, 0.0]
    assert output[-1][1]['pos'] == [-35.88936905187977, 31.04336783591637, 0.0]
