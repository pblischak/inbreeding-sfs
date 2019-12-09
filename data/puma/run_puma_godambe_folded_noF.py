#!/usr/bin/env python3

import dadi
import dadi.Godambe
import numpy as np

def divergence_noF(params, ns, pts):
    """
    Divergence only model.
    """
    nuTX,nuFL,T1,T2 = params

    xx = yy = dadi.Numerics.default_grid(pts)
    n1,n2 = ns[0],ns[1]

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop(phi,xx,T1,nu=nuTX)
    phi = dadi.PhiManip.phi_1D_to_2D(xx,phi)

    phi = dadi.Integration.two_pops(phi,xx,T2,nu1=nuTX,nu2=nuFL)
    sfs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,yy))
    return sfs

if __name__ == '__main__':
    # Read in the bootstrapped spectra and then fold
    all_boot = [dadi.Spectrum.from_file("puma_test_boot{}.fs".format(i)) for i in range(100)]
    all_boot = [sfs.fold() for sfs in all_boot]

    # Now read in the original data set
    data = dadi.Spectrum.from_file("puma.fs")
    data = data.fold()
    pts_l = [40,50,60]
    func = divergence_noF
    popt = [0.197150919878390,0.0100564473768224,0.0371335637849756,0.0113941805000981]

    func_ex = dadi.Numerics.make_extrap_log_func(func)

    uncerts = dadi.Godambe.GIM_uncert(func_ex, pts_l, all_boot, popt, data, log=True, multinom=True)
    print('Estimated parameter standard deviations from GIM: {0}'.format(uncerts))
    theta = 2718260.477394347 # estimated from model
    L = 2564692624
    mu = 2.2e-9
    g  = 3 # generation time
    Nref = theta / L / mu / g / 4
    print("\n\nConverting parameters to actual units...\n")
    print("Using the following values for conversion:")
    print("  Sequence length = {}".format(L))
    print("  Mutation rate   = {}".format(mu))
    print("  Generation time = {}".format(g))

    print("\nNref = {} ({}--{})".format(Nref, np.exp(np.log(theta)-1.96*uncerts[-1]) / L / mu / g / 4, np.exp(np.log(theta)+1.96*uncerts[-1]) / L / mu / g / 4))
    print("N_TX = {} ({}--{})".format(popt[0]*Nref, np.exp(np.log(popt[0])-1.96*uncerts[0])*Nref, np.exp(np.log(popt[0])+1.96*uncerts[0])*Nref))
    print("N_FL = {} ({}--{})".format(popt[1]*Nref, np.exp(np.log(popt[1])-1.96*uncerts[1])*Nref, np.exp(np.log(popt[1])+1.96*uncerts[1])*Nref))
    print("T1   = {} ({}--{})".format(popt[2]*2*Nref, np.exp(np.log(popt[2])-1.96*uncerts[2])*2*Nref, np.exp(np.log(popt[2])+1.96*uncerts[2])*2*Nref))
    print("T2   = {} ({}--{})".format(popt[3]*2*Nref, np.exp(np.log(popt[3])-1.96*uncerts[3])*2*Nref, np.exp(np.log(popt[3])+1.96*uncerts[3])*2*Nref))
