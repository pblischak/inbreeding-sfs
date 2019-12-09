#!/usr/bin/env python3

import dadi
import dadi.Godambe
import numpy as np

def three_epoch_noF(params, ns, pts):
    """
    params = (nuB,nuF,TB,TF)
    ns = (n1,)

    nuB: Ratio of bottleneck population size to ancient pop size
    nuF: Ratio of contemporary to ancient pop size
    TB: Length of bottleneck (in units of 2*Na generations)
    TF: Time since bottleneck recovery (in units of 2*Na generations)

    n1: Number of samples in resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nuB,nuF,TB,TF = params

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)

    phi = dadi.Integration.one_pop(phi, xx, TB, nuB)
    phi = dadi.Integration.one_pop(phi, xx, TF, nuF)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
    return fs

if __name__ == '__main__':
    # Read in the bootstrapped spectra and then fold
    all_boot = [dadi.Spectrum.from_file("boot/cabbage_test_boot{}.fs".format(i)) for i in range(100)]
    all_boot = [sfs.fold() for sfs in all_boot]

    # Now read in the original data set
    data = dadi.Spectrum.from_file("cabbage.fs")
    data = data.fold()
    pts_l = [100,110,120]
    func = three_epoch_noF
    popt = [6.4524615672350958,0.0309347139612217,0.153264805381591,0.00100000000000000]

    func_ex = dadi.Numerics.make_extrap_log_func(func)

    # Calculate parameter uncertainties
    uncerts = dadi.Godambe.GIM_uncert(func_ex, pts_l, all_boot, popt, data, log=True, multinom=True)
    print("Estimated parameter standard deviations from GIM: {}".format(uncerts))

    # Set conversion parameters
    theta = 472576.2690794509 # estimated from model
    L = 411560319
    mu = 1.5e-8
    g  = 1 # generation time
    Nref = theta / L / mu / g / 4
    print("\n\nConverting parameters to actual units...\n")
    print("Using the following values for conversion:")
    print("  Sequence length = {}".format(L))
    print("  Mutation rate =   {}".format(mu))
    print("  Generation time = {}".format(g))

    print("\nNref = {:0.2f} ({:0.2f}--{:0.2f})".format(Nref, np.exp(np.log(theta)-1.96*uncerts[-1]) / L / mu / g / 4, np.exp(np.log(theta)+1.96*uncerts[-1]) / L / mu / g / 4))
    print("N1   = {:0.2f} ({:0.2f}--{:0.2f})".format(popt[0]*Nref, np.exp(np.log(popt[0])-1.96*uncerts[0])*Nref, np.exp(np.log(popt[0])+1.96*uncerts[0])*Nref))
    print("N2   = {:0.2f} ({:0.2f}--{:0.2f})".format(popt[1]*Nref, np.exp(np.log(popt[1])-1.96*uncerts[1])*Nref, np.exp(np.log(popt[1])+1.96*uncerts[1])*Nref))
    print("T1   = {:0.2f} ({:0.2f}--{:0.2f})".format(popt[2]*2*Nref, np.exp(np.log(popt[2])-1.96*uncerts[2])*2*Nref, np.exp(np.log(popt[2])+1.96*uncerts[2])*2*Nref))
    print("T2   = {:0.2f} ({:0.2f}--{:0.2f})".format(popt[3]*2*Nref, np.exp(np.log(popt[3])-1.96*uncerts[3])*2*Nref, np.exp(np.log(popt[3])+1.96*uncerts[3])*2*Nref))
