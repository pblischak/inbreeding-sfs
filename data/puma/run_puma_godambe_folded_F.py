#!/usr/bin/env python3

import dadi
import dadi.Godambe
from scipy import stats
import numpy as np

def divergence(params, ns, pts):
    """
    Divergence only model.
    """
    nuTX,nuFL,T1,T2,F1,F2 = params

    xx = yy = dadi.Numerics.default_grid(pts)
    n1,n2 = ns[0],ns[1]

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop(phi,xx,T1,nu=nuTX)
    phi = dadi.PhiManip.phi_1D_to_2D(xx,phi)

    phi = dadi.Integration.two_pops(phi,xx,T2,nu1=nuTX,nu2=nuFL)
    sfs = dadi.Spectrum.from_phi_inbreeding(phi, (n1,n2), (xx,yy), (F1,F2), (2,2))
    return sfs

if __name__ == '__main__':
    # Read in the bootstrapped spectra and then fold
    all_boot = [dadi.Spectrum.from_file("puma_test_boot{}.fs".format(i)) for i in range(100)]
    all_boot = [sfs.fold() for sfs in all_boot]

    # Now read in the original data set
    data = dadi.Spectrum.from_file("puma.fs")
    data = data.fold()
    pts_l = [40,50,60]
    func = divergence
    popt = [0.5439309410920344,0.0120104498902649,0.3154710478935099,0.0100122945339300,0.43951783360569285,0.606535498834911]

    func_ex = dadi.Numerics.make_extrap_log_func(func)

    uncerts = dadi.Godambe.GIM_uncert(func_ex, pts_l, all_boot, popt, data, log=True, multinom=True)
    print('Estimated parameter standard deviations from GIM: {0}'.format(uncerts))

    llik_F   = -318058.079230279
    llik_noF = -453003.047885215
    adj      = dadi.Godambe.LRT_adjust(func_ex, pts_l, all_boot, popt, data, nested_indices=[4,5], multinom=True)
    D_adj    = 2 * adj * (llik_F - llik_noF)
    print("\n\nLikelihood ratio test:")
    print("  No inbreeding likelihood = {:0.3f}".format(llik_noF))
    print("  Inbreeding likelihood    = {:0.3f}".format(llik_F))
    print("  LRT adjustment           = {:0.3f}".format(adj))
    print("  LRT statistic            = {:0.3f}".format(D_adj))
    print("  LRT p-value              = {:0.3f}".format(dadi.Godambe.sum_chi2_ppf(D_adj, (1.0/3.0,1.0/3.0,1.0/3.0))))

    # Set conversion parameters
    theta = 2939439.410209113 # estimated from model
    L = 2564692624
    mu = 2.2e-9
    g  = 3 # generation time
    Nref = theta / L / mu / g / 4
    print("\n\nConverting parameters to actual units...\n")
    print("Using the following values for conversion:")
    print("  Sequence length = {}".format(L))
    print("  Mutation rate   = {}".format(mu))
    print("  Generation time = {}".format(g))

    print("\nNref = {:0.2f} ({:0.2f}--{:0.2f})".format(Nref, np.exp(np.log(theta)-1.96*uncerts[-1]) / L / mu / g / 4, np.exp(np.log(theta)+1.96*uncerts[-1]) / L / mu / g / 4))
    print("N_TX = {:0.2f} ({:0.2f}--{:0.2f})".format(popt[0]*Nref, np.exp(np.log(popt[0])-1.96*uncerts[0])*Nref, np.exp(np.log(popt[0])+1.96*uncerts[0])*Nref))
    print("N_FL = {:0.2f} ({:0.2f}--{:0.2f})".format(popt[1]*Nref, np.exp(np.log(popt[1])-1.96*uncerts[1])*Nref, np.exp(np.log(popt[1])+1.96*uncerts[1])*Nref))
    print("T1   = {:0.2f} ({:0.2f}--{:0.2f})".format(popt[2]*2*Nref, np.exp(np.log(popt[2])-1.96*uncerts[2])*2*Nref, np.exp(np.log(popt[2])+1.96*uncerts[2])*2*Nref))
    print("T2   = {:0.2f} ({:0.2f}--{:0.2f})".format(popt[3]*2*Nref, np.exp(np.log(popt[3])-1.96*uncerts[3])*2*Nref, np.exp(np.log(popt[3])+1.96*uncerts[3])*2*Nref))
    print("F_TX = {:0.3f} ({:0.3f}--{:0.3f})".format(popt[4], np.exp(np.log(popt[4])-1.96*uncerts[4]), np.exp(np.log(popt[4])+1.96*uncerts[4])))
    print("F_FL = {:0.3f} ({:0.3f}--{:0.3f})".format(popt[5], np.exp(np.log(popt[5])-1.96*uncerts[5]), np.exp(np.log(popt[5])+1.96*uncerts[5])))
