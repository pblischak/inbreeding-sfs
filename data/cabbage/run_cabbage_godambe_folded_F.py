#!/usr/bin/env python3

import dadi
import dadi.Godambe

def bottlegrowth(params, ns, pts):
    nuBot,nuCur,T,F = params

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)

    nu_func = lambda t: nuBot * (nuCur/nuBot) ** (t/T)
    phi = dadi.Integration.one_pop(phi,xx,T,nu=nu_func)
    sfs = dadi.Spectrum.from_phi_inbreeding(phi, ns, (xx,), (F,), (2,))
    return sfs

def expand_bottlegrowth(params,ns,pts):
    """
    This is the model that we used.
    """
    nuExp,nuBot,nuCur,T1,T2,F = params

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)

    phi = dadi.Integration.one_pop(phi,xx,T1,nu=nuExp)
    nu_func = lambda t: nuBot * (nuCur/nuBot) ** (t/T2)
    phi = dadi.Integration.one_pop(phi,xx,T2,nu=nu_func)
    sfs = dadi.Spectrum.from_phi_inbreeding(phi, ns, (xx,), (F,), (2,))
    return sfs

if __name__ == '__main__':
    # Read in the bootstrapped spectra and then fold
    all_boot = [dadi.Spectrum.from_file("boot/cabbage_boot{}.fs".format(i)) for i in range(100)]
    all_boot = [sfs.fold() for sfs in all_boot]

    # Now read in the original data set
    data = dadi.Spectrum.from_file("cabbage.fs")
    data = data.fold()
    pts_l = [100,110,120]
    func = expand_bottlegrowth
    popt = [1.8034444276865402,1.5406082627009385,7.128657409544192,0.4769670348935206,0.0143092989018097,0.577605971746217]

    func_ex = dadi.Numerics.make_extrap_log_func(func)

    # If we've already calculated uncertainties, there's no reason to do it again
    # Just assign to the uncerts variable from a previous run
    calc_uncerts = False
    if calc_uncerts:
        uncerts = dadi.Godambe.GIM_uncert(func_ex, pts_l, all_boot, popt, data,
                                          multinom=True)
    else:
        uncerts = [5.51189193e-03,1.07895751e+00,2.38999449e+00,7.50757659e-03,5.70797205e-03,9.23088358e-04,1.16563746e+03]
    print('Estimated parameter standard deviations from GIM: {0}'.format(uncerts))

    # Set conversion parameters
    theta = 430717.205081515 # estimated from model
    L = 554977060 # Sequence length
    mu = 1.5e-8
    g  = 1 # generation time
    Nref = theta / L / mu / g / 4
    print("\n\nConverting parameters to actual units...\n")
    print("Using the following values for conversion:")
    print("  Sequence length = {}".format(L))
    print("  Mutation rate =   {}".format(mu))
    print("  Generation time = {}".format(g))

    print("\nNref = {} (CI: {} -- {})".format(Nref, (theta-uncerts[-1]) / L / mu / g / 4, (theta+uncerts[-1]) / L / mu / g / 4))
    print("NExp = {} (CI: {} -- {})".format(popt[0]*Nref, (popt[0]-uncerts[0])*Nref, (popt[0]+uncerts[0])*Nref))
    print("NBot = {} (CI: {} -- {})".format(popt[1]*Nref, (popt[1]-uncerts[1])*Nref, (popt[1]+uncerts[1])*Nref))
    print("NCur = {} (CI: {} -- {})".format(popt[2]*Nref, (popt[2]-uncerts[2])*Nref, (popt[2]+uncerts[2])*Nref))
    print("T1    = {} (CI: {} -- {})".format(popt[3]*2*Nref, (popt[3]-uncerts[3])*2*Nref, (popt[3]+uncerts[3])*2*Nref))
    print("T2    = {} (CI: {} -- {})".format(popt[4]*2*Nref, (popt[4]-uncerts[4])*2*Nref, (popt[4]+uncerts[4])*2*Nref))
    print("Fis  = {} (CI: {} -- {})".format(popt[5], popt[5]-uncerts[5], popt[5]+uncerts[5]))
