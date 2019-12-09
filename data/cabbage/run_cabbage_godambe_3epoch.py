#!/usr/bin/env python3

import dadi
import dadi.Godambe
from scipy import stats
import numpy as np

def three_epoch(params, ns, pts):
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
    nuB,nuF,TB,TF,F = params

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)

    phi = dadi.Integration.one_pop(phi, xx, TB, nuB)
    phi = dadi.Integration.one_pop(phi, xx, TF, nuF)

    fs = dadi.Spectrum.from_phi_inbreeding(phi, ns, (xx,), (F,), (2,))
    return fs

if __name__ == '__main__':
    # Read in the bootstrapped spectra and then fold
    all_boot = [dadi.Spectrum.from_file("boot/cabbage_test_boot{}.fs".format(i)) for i in range(100)]
    all_boot = [sfs.fold() for sfs in all_boot]

    # Now read in the original data set
    data = dadi.Spectrum.from_file("cabbage.fs")
    data = data.fold()
    pts_l = [100,110,120]
    func = three_epoch
    popt = [1.810449088130342,12.2790194725467110,0.47393521119534737,0.00921096365957015,0.577870722504928]

    func_ex = dadi.Numerics.make_extrap_log_func(func)

    # Calculate parameter uncertainties
    uncerts = dadi.Godambe.GIM_uncert(func_ex, pts_l, all_boot, popt, data, log=True, multinom=True)
    print('Estimated parameter standard deviations from GIM: {0}'.format(uncerts))

    llik_F   = -4281.14456624734
    llik_noF = -24330.4027006874
    adj      = dadi.Godambe.LRT_adjust(func_ex, pts_l, all_boot, popt, data, nested_indices=[4], multinom=True)
    D_adj    = 2 * adj * (llik_F - llik_noF)
    print("\n\nLikelihood ratio test:")
    print("  No inbreeding likelihood = {}".format(llik_noF))
    print("  Inbreeding likelihood    = {}".format(llik_F))
    print("  LRT adjustment           = {}".format(adj))
    print("  LRT statistic            = {}".format(D_adj))
    print("  LRT p-value (X^2, 1 df)  = {}".format(dadi.Godambe.sum_chi2_ppf(D_adj, (0.5,0.5))))


    # Set conversion parameters
    theta = 431524.124873428 # estimated from model
    L = 411560319
    mu = 1.5e-8
    g  = 1 # generation time
    Nref = theta / L / mu / g / 4
    print("\n\nConverting parameters to actual units...\n")
    print("Using the following values for conversion:")
    print("  Sequence length = {}".format(L))
    print("  Mutation rate   = {}".format(mu))
    print("  Generation time = {}".format(g))

    print("\nNref = {} ({}--{})".format(Nref, np.exp(np.log(theta)-1.96*uncerts[-1]) / L / mu / g / 4, np.exp(np.log(theta)+1.96*uncerts[-1]) / L / mu / g / 4))
    print("N1   = {} ({}--{})".format(popt[0]*Nref, np.exp(np.log(popt[0])-1.96*uncerts[0])*Nref, np.exp(np.log(popt[0])+1.96*uncerts[0])*Nref))
    print("N2   = {} ({}--{})".format(popt[1]*Nref, np.exp(np.log(popt[1])-1.96*uncerts[1])*Nref, np.exp(np.log(popt[1])+1.96*uncerts[1])*Nref))
    print("T1   = {} ({}--{})".format(popt[2]*2*Nref, np.exp(np.log(popt[2])-1.96*uncerts[2])*2*Nref, np.exp(np.log(popt[2])+1.96*uncerts[2])*2*Nref))
    print("T2   = {} ({}--{})".format(popt[3]*2*Nref, np.exp(np.log(popt[3])-1.96*uncerts[3])*2*Nref, np.exp(np.log(popt[3])+1.96*uncerts[3])*2*Nref))
    print("Fis  = {} ({}--{})".format(popt[4], np.exp(np.log(popt[4])-1.96*uncerts[4]), np.exp(np.log(popt[4])+1.96*uncerts[4])))
