#!/usr/bin/env python3

import dadi
import dadi.Godambe
from scipy import stats
import numpy as np

def divergence(params, ns, pts):
    """
    Divergence only model.
    """
    nuTX,nuFL,T,F1,F2 = params

    xx = yy = dadi.Numerics.default_grid(pts)
    n1,n2 = ns[0],ns[1]

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx,phi)

    phi = dadi.Integration.two_pops(phi,xx,T,nu1=nuTX,nu2=nuFL)
    sfs = dadi.Spectrum.from_phi_inbreeding(phi, (n1,n2), (xx,yy), (F1,F2), (2,2))
    return sfs


if __name__ == '__main__':
    # Read in the bootstrapped spectra and then fold
    all_boot = [dadi.Spectrum.from_file("boot/puma_test_boot{}.fs".format(i)) for i in range(100)]
    all_boot = [sfs.fold() for sfs in all_boot]

    # Now read in the original data set
    data = dadi.Spectrum.from_file("puma.fs")
    data = data.fold()
    pts_l = [40,50,60]
    func = divergence
    popt = [0.45043163,0.01252053,0.01002769,0.4814076365,0.6278800995]

    func_ex = dadi.Numerics.make_extrap_log_func(func)

    uncerts = dadi.Godambe.GIM_uncert(func_ex, pts_l, all_boot, popt, data, log=True, multinom=True)
    print('Estimated parameter standard deviations from GIM: {0}'.format(uncerts))

    llik_F   = -343218.201692238
    llik_noF = -499397.564085949
    adj      = dadi.Godambe.LRT_adjust(func_ex, pts_l, all_boot, popt, data, nested_indices=[3,4], multinom=True)
    D_adj    = 2 * adj * (llik_F - llik_noF)
    print("\n\nLikelihood ratio test:")
    print("  No inbreeding likelihood = {}".format(llik_noF))
    print("  Inbreeding likelihood    = {}".format(llik_F))
    print("  LRT adjustment           = {}".format(adj))
    print("  LRT statistic            = {}".format(D_adj))
    print("  LRT p-value              = {}".format(dadi.Godambe.sum_chi2_ppf(D_adj, (1.0/3.0,1.0/3.0,1.0/3.0))))

    # Set conversion parameters
    theta = 2307096.160 # estimated from model
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
    print("T    = {} ({}--{})".format(popt[2]*2*Nref, np.exp(np.log(popt[2])-1.96*uncerts[2])*2*Nref, np.exp(np.log(popt[2])+1.96*uncerts[2])*2*Nref))
    print("F_TX = {} ({}--{})".format(popt[3], np.exp(np.log(popt[3])-1.96*uncerts[3]), np.exp(np.log(popt[3])+1.96*uncerts[3])))
    print("F_FL = {} ({}--{})".format(popt[4], np.exp(np.log(popt[4])-1.96*uncerts[4]), np.exp(np.log(popt[4])+1.96*uncerts[4])))
