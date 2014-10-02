import numpy as np
import matplotlib as mpl
import _orbit
from mcorbit_plot_utils import *
from mcorbit_data import *


def draw_best_model(name,
                    streamname,
                    folder='./',
                    radec=True,
                    datatype=3,
                    errors=False):

    mpl.rcParams.update({'font.size': 12})

    print "File name: chain%s.dat" % name

    prob,\
        mass_cluster,\
        pm_mu_delta,\
        pm_mu_alphacosdelta,\
        mass_halo,\
        distance_cluster,\
        mass_loss_rate,\
        q_halo,\
        r_halo,\
        tpast,\
        sigma_x,\
        sigma_v,\
        sigma_vx,\
        sigma_mu,\
        rgalsun,\
        vLSR = np.loadtxt("chain%s.dat" % name, unpack=True)

    ind = np.argmax(prob)
    print "Best model in chain is model no %i with P = %f" % (ind, prob[ind])

    _orbit.orbit(mass_cluster[ind],
                 pm_mu_delta[ind],
                 pm_mu_alphacosdelta[ind],
                 mass_halo[ind],
                 distance_cluster[ind],
                 mass_loss_rate[ind],
                 q_halo[ind],
                 r_halo[ind],
                 tpast[ind],
                 rgalsun[ind],
                 vLSR[ind],
                 sigma_x[ind],
                 sigma_v[ind],
                 sigma_vx[ind],
                 sigma_mu[ind],
                 1)

    stream_data = np.loadtxt("temp.dat")

    ODdata = NGC5466_OD()
    VRdata = NGC5466_VR()

    if datatype == 1:
        VRdata = VRdata[VRdata[:, 7] == 1]
        ODdata = ODdata[ODdata[:, 4] == 1]
    elif datatype == 2:
        VRdata = VRdata[VRdata[:, 7] == 2]
        ODdata = ODdata[ODdata[:, 4] == 2]
    elif datatype == 3:  # All data
        VRdata = VRdata[(VRdata[:, 7] == 2) | (VRdata[:, 7] == 1)]
        ODdata = ODdata[(ODdata[:, 4] == 2) | (ODdata[:, 4] == 1)]

    clusterpos = ([42.15, 73.59], [211.36370833, 28.53444444])[radec]

    plot_model(stream_data,
               VRdata,
               ODdata,
               'Best',
               name=name,
               streamname=streamname,
               folder=folder,
               radec=radec,
               clusterpos=clusterpos,
               errors=errors)


def draw_model_no(name,
                  modelnumber,
                  streamname,
                  folder='./',
                  radec=True,
                  datatype=3,
                  errors=False):
    mpl.rcParams.update({'font.size': 12})

    print "File name: chain%s.dat" % name
    print "Model number: %i" % modelnumber

    prob,\
        mass_cluster,\
        pm_mu_delta,\
        pm_mu_alphacosdelta,\
        mass_halo,\
        distance_cluster,\
        mass_loss_rate,\
        q_halo,\
        r_halo,\
        tpast,\
        sigma_x,\
        sigma_v,\
        sigma_vx,\
        sigma_mu,\
        rgalsun,\
        vLSR = np.loadtxt("chain%s.dat" % name, unpack=True)

    ind = np.argmax(prob)
    print "Best model in chain is model no %i with P = %f" % (ind, prob[ind])

    _orbit.orbit(mass_cluster[modelnumber],
                 pm_mu_delta[modelnumber],
                 pm_mu_alphacosdelta[modelnumber],
                 mass_halo[modelnumber],
                 distance_cluster[modelnumber],
                 mass_loss_rate[modelnumber],
                 q_halo[modelnumber],
                 r_halo[modelnumber],
                 tpast[modelnumber],
                 rgalsun[modelnumber],
                 vLSR[modelnumber],
                 sigma_x[modelnumber],
                 sigma_v[modelnumber],
                 sigma_vx[modelnumber],
                 sigma_mu[modelnumber],
                 1)

    stream_data = np.loadtxt("temp.dat")
    print(stream_data.shape)

    ODdata = NGC5466_OD()
    VRdata = NGC5466_VR()

    if datatype == 1:
        VRdata = VRdata[VRdata[:, 7] == 1]
        ODdata = ODdata[ODdata[:, 4] == 1]
    elif datatype == 2:
        VRdata = VRdata[VRdata[:, 7] == 2]
        ODdata = ODdata[ODdata[:, 4] == 2]
    elif datatype == 3:
        VRdata = VRdata[(VRdata[:, 7] == 2) | (VRdata[:, 7] == 1)]
        ODdata = ODdata[(ODdata[:, 4] == 2) | (ODdata[:, 4] == 1)]

    clusterpos = ([42.15, 73.59], [211.36370833, 28.53444444])[radec]

    plot_model(stream_data,
               VRdata,
               ODdata,
               str(modelnumber),
               name=name,
               streamname=streamname,
               folder=folder,
               radec=radec,
               clusterpos=clusterpos,
               errors=errors)


def draw_median_model(name,
                      streamname,
                      folder='./',
                      radec=True,
                      datatype=3,
                      errors=False):

    mpl.rcParams.update({'font.size': 12})

    print "File name: chain%s.dat" % name

    prob,\
        mass_cluster,\
        pm_mu_delta,\
        pm_mu_alphacosdelta,\
        mass_halo,\
        distance_cluster,\
        mass_loss_rate,\
        q_halo,\
        r_halo,\
        tpast,\
        sigma_x,\
        sigma_v,\
        sigma_vx,\
        sigma_mu,\
        rgalsun,\
        vLSR = np.loadtxt("chain%s.dat" % name, unpack=True)

    _orbit.orbit(np.median(mass_cluster),
                 np.median(pm_mu_delta),
                 np.median(pm_mu_alphacosdelta),
                 np.median(mass_halo),
                 np.median(distance_cluster),
                 np.median(mass_loss_rate),
                 np.median(q_halo),
                 np.median(r_halo),
                 np.median(tpast),
                 np.median(rgalsun),
                 np.median(vLSR),
                 np.median(sigma_x),
                 np.median(sigma_v),
                 np.median(sigma_vx),
                 np.median(sigma_mu),
                 1)

    stream_data = np.loadtxt("temp.dat")

    ODdata = NGC5466_OD()
    VRdata = NGC5466_VR()

    if datatype == 1:
        VRdata = VRdata[VRdata[:, 7] == 1]
        ODdata = ODdata[ODdata[:, 4] == 1]
    elif datatype == 2:
        VRdata = VRdata[VRdata[:, 7] == 2]
        ODdata = ODdata[ODdata[:, 4] == 2]
    elif datatype == 3:
        VRdata = VRdata[(VRdata[:, 7] == 2) | (VRdata[:, 7] == 1)]
        ODdata = ODdata[(ODdata[:, 4] == 2) | (ODdata[:, 4] == 1)]

    clusterpos = ([42.15, 73.59], [211.36370833, 28.53444444])[radec]

    plot_model(stream_data,
               VRdata,
               ODdata,
               'Median',
               name=name,
               streamname=streamname,
               folder=folder,
               radec=radec,
               clusterpos=clusterpos,
               errors=errors)


def draw_last_model(folder='./', radec=True, datatype=3, errors=False):
    mpl.rcParams.update({'font.size': 12})

    stream_data = np.loadtxt("temp.dat")

    ODdata = NGC5466_OD()
    VRdata = NGC5466_VR()

    if datatype == 1:
        VRdata = VRdata[VRdata[:, 7] == 1]
        ODdata = ODdata[ODdata[:, 4] == 1]
    elif datatype == 2:
        VRdata = VRdata[VRdata[:, 7] == 2]
        ODdata = ODdata[ODdata[:, 4] == 2]
    elif datatype == 3:
        VRdata = VRdata[(VRdata[:, 7] == 2) | (VRdata[:, 7] == 1)]
        ODdata = ODdata[(ODdata[:, 4] == 2) | (ODdata[:, 4] == 1)]

    clusterpos = ([42.15, 73.59], [211.36370833, 28.53444444])[radec]

    plot_model(stream_data,
               VRdata,
               ODdata,
               'Last',
               folder=folder,
               radec=radec,
               clusterpos=clusterpos,
               errors=errors)


def draw_model(mass_cluster,
               pm_mu_delta,
               pm_mu_alphacosdelta,
               mass_halo,
               distance_cluster,
               mass_loss_rate,
               q_halo,
               r_halo,
               tpast,
               rgalsun,
               vLSR,
               folder='./',
               radec=True,
               datatype=3,
               errors=False):

    mpl.rcParams.update({'font.size': 12})

    _orbit.orbit(mass_cluster,
                 pm_mu_delta,
                 pm_mu_alphacosdelta,
                 mass_halo,
                 distance_cluster,
                 mass_loss_rate,
                 q_halo,
                 r_halo,
                 tpast,
                 rgalsun,
                 vLSR,
                 0.0,
                 0.0,
                 0.0,
                 0.0,
                 1)

    print("Done with modeling")

    stream_data = np.loadtxt("temp.dat")

    ODdata = NGC5466_OD()
    VRdata = NGC5466_VR()

    if datatype == 1:
        VRdata = VRdata[VRdata[:, 7] == 1]
        ODdata = ODdata[ODdata[:, 4] == 1]
    elif datatype == 2:
        VRdata = VRdata[VRdata[:, 7] == 2]
        ODdata = ODdata[ODdata[:, 4] == 2]
    elif datatype == 3:
        VRdata = VRdata[(VRdata[:, 7] == 2) | (VRdata[:, 7] == 1)]
        ODdata = ODdata[(ODdata[:, 4] == 2) | (ODdata[:, 4] == 1)]

    clusterpos = ([42.15, 73.59], [211.36370833, 28.53444444])[radec]

    plot_model(stream_data,
               VRdata,
               ODdata,
               'Temp',
               name='NGC5466',
               streamname='NGC5466',
               folder=folder,
               radec=radec,
               clusterpos=clusterpos,
               errors=errors)


def draw_temp_model(folder='./', radec=True, datatype=3, errors=False):
    mpl.rcParams.update({'font.size': 12})

    stream_data = np.loadtxt("temp.dat")

    ODdata = NGC5466_OD()
    VRdata = NGC5466_VR()

    if datatype == 1:
        VRdata = VRdata[VRdata[:, 7] == 1]
        ODdata = ODdata[ODdata[:, 4] == 1]
    elif datatype == 2:
        VRdata = VRdata[VRdata[:, 7] == 2]
        ODdata = ODdata[ODdata[:, 4] == 2]
    elif datatype == 3:
        VRdata = VRdata[(VRdata[:, 7] == 2) | (VRdata[:, 7] == 1)]
        ODdata = ODdata[(ODdata[:, 4] == 2) | (ODdata[:, 4] == 1)]

    clusterpos = ([42.15, 73.59], [211.36370833, 28.53444444])[radec]

    plot_model(stream_data,
               VRdata,
               ODdata,
               'Temp',
               name='NGC5466',
               streamname='NGC5466',
               folder=folder,
               radec=radec,
               clusterpos=clusterpos,
               errors=errors)


def plot_all(name, folder='./'):

    streamname = "NGC5466"
    gal_latitude = 42.15
    gal_longitude = 73.59
    plot_convergence(name, streamname, folder=folder)
    param_histo_stream(name, streamname, folder=folder)
    param_histo_Halo(name, folder=folder)
    param_histo_acc(name,
                    streamname,
                    gal_latitude,
                    gal_longitude,
                    folder=folder)
    plot_mstream(name, streamname, folder=folder)
    plot_mudelta(name, folder=folder)
    plot_mualpha(name, folder=folder)
    plot_distance(name, folder=folder)
    plot_mhalo(name, folder=folder)
    plot_rhalo(name, folder=folder)
    plot_qhalo(name, folder=folder)
    plot_dmdt(name, folder=folder)
    plot_vLSR(name, folder=folder)
    plot_rgalsun(name, folder=folder)
    plot_triangle_cluster(name, streamname, folder=folder)
    plot_triangle_halo(name, folder=folder)
    plot_triangle_orbital(name, streamname, folder=folder)
    plot_triangle_derived(name,
                          streamname,
                          gal_latitude,
                          gal_longitude,
                          folder=folder)
