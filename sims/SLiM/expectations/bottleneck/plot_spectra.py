#!/usr/bin/env python3

import dadi
import numpy as np
import matplotlib.pyplot as plt

def bottleneck(params, ns, pts):
    """
    Model 2
    -------
    Bottleneck followed by growth.
    Parameters
    ----------
    nu0: Relative size of pop after bottleneck.
    T:   Time of bottleneck.
    F:   Inbreeding coefficient (0 < F < 1).
    """
    nu0, T, F = params

    # check that F is in range
    if F <= 0.0:
        F = 1e-4
    elif F >= 1.0:
        F = 1.0 - 1e-4
    else:
        pass

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)

    nu_func = lambda t: nu0 * (1.0 / nu0) ** (t/T)
    phi = dadi.Integration.one_pop(phi, xx, T, nu_func)

    fs = dadi.Spectrum.from_phi_inbreeding(phi, ns, (xx,), (F,), (2,))
    return fs

if __name__ == "__main__":
    theta = 10000
    pts_l = [70,80,90]
    func_ex = dadi.Numerics.make_extrap_log_func(bottleneck)
    
    # Start reading in spectra
    # F = 0.1
    #dadi_F10 = theta / (1+0.1) * func_ex([0.25,0.2,0.1],[50], pts_l)
    dadi_F10 = theta / (1+0.1) * func_ex([0.25,0.2 * (1+0.1),0.1],[50], pts_l)
    slim_F10 = dadi.Spectrum.from_file("SLiM_F0.1_mean_bottleneck.fs")
    print("F=0.1 RMSD:  {}".format(np.sqrt(np.mean((dadi_F10 - slim_F10)**2))))
    print("F=0.1 RRMSD: {}\n".format(np.sqrt(np.mean((dadi_F10 - slim_F10)**2))/np.sum(slim_F10)))
    plt.plot(dadi_F10, '-ob')
    plt.plot(slim_F10, '-og')
    plt.show()
    
    # F = 0.25
    #dadi_F25 = theta / (1+0.25) * func_ex([0.25,0.2,0.25],[50], pts_l)
    dadi_F25 = theta / (1+0.25) * func_ex([0.25,0.2*(1+0.25),0.25],[50], pts_l)
    slim_F25 = dadi.Spectrum.from_file("SLiM_F0.25_mean_bottleneck.fs")
    print("F=0.25 RMSD:  {}".format(np.sqrt(np.mean((dadi_F25 - slim_F25)**2))))
    print("F=0.25 RRMSD: {}\n".format(np.sqrt(np.mean((dadi_F25 - slim_F25)**2))/np.sum(slim_F25)))
    plt.plot(dadi_F25, '-ob')
    plt.plot(slim_F25, '-og')
    plt.show()
    
    # F = 0.5
    #dadi_F50 = theta / (1+0.5) * func_ex([0.25,0.2,0.5],[50], pts_l)
    dadi_F50 = theta / (1+0.5) * func_ex([0.25,0.2*(1+0.5),0.5],[50], pts_l)
    slim_F50 = dadi.Spectrum.from_file("SLiM_F0.5_mean_bottleneck.fs")
    print("F=0.5 RMSD:  {}".format(np.sqrt(np.mean((dadi_F50 - slim_F50)**2))))
    print("F=0.5 RRMSD: {}\n".format(np.sqrt(np.mean((dadi_F50 - slim_F50)**2))/np.sum(slim_F50)))
    plt.plot(dadi_F50, '-ob')
    plt.plot(slim_F50, '-og')
    plt.show()
    
    # F = 0.75
    #dadi_F75 = theta / (1+0.75) * func_ex([0.25,0.2,0.75],[50], pts_l)
    dadi_F75 = theta / (1+0.75) * func_ex([0.25,0.2*(1+0.75),0.75],[50], pts_l)
    slim_F75 = dadi.Spectrum.from_file("SLiM_F0.75_mean_bottleneck.fs")
    print("F=0.75 RMSD:  {}".format(np.sqrt(np.mean((dadi_F75 - slim_F75)**2))))
    print("F=0.75 RRMSD: {}\n".format(np.sqrt(np.mean((dadi_F75 - slim_F75)**2))/np.sum(slim_F75)))
    plt.plot(dadi_F75, '-ob')
    plt.plot(slim_F75, '-og')
    plt.show()
    
    # F = 0.9
    #dadi_F90 = theta / (1+0.9) * func_ex([0.25,0.2,0.9],[50], pts_l)
    dadi_F90 = theta / (1+0.9) * func_ex([0.25,0.2*(1+0.9),0.9],[50], pts_l)
    slim_F90 = dadi.Spectrum.from_file("SLiM_F0.9_mean_bottleneck.fs")
    print("F=0.9 RMSD:  {}".format(np.sqrt(np.mean((dadi_F90 - slim_F90)**2))))
    print("F=0.9 RRMSD: {}\n".format(np.sqrt(np.mean((dadi_F90 - slim_F90)**2))/np.sum(slim_F90)))
    plt.plot(dadi_F90, '-ob')
    plt.plot(slim_F90, '-og')
    plt.show()
