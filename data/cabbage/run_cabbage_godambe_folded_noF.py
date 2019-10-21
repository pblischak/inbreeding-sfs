#!/usr/bin/env python3

import dadi
import dadi.Godambe

def bottlegrowth_noF(params, ns, pts):
    nuBot,nuCur,T = params

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)

    nu_func = lambda t: nuBot * (nuCur/nuBot) ** (t/T)
    phi = dadi.Integration.one_pop(phi,xx,T,nu=nu_func)
    sfs = dadi.Spectrum.from_phi(phi, ns, (xx,))
    return sfs

def expand_bottlegrowth_noF(params,ns,pts):
    """
    This is the model that we used.
    """
    nuExp,nuBot,nuCur,T1,T2 = params

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)

    phi = dadi.Integration.one_pop(phi,xx,T1,nu=nuExp)
    nu_func = lambda t: nuBot * (nuCur/nuBot) ** (t/T2)
    phi = dadi.Integration.one_pop(phi,xx,T2,nu=nu_func)
    sfs = dadi.Spectrum.from_phi(phi, ns, (xx,))
    return sfs

if __name__ == '__main__':
    # Read in the bootstrapped spectra and then fold
    all_boot = [dadi.Spectrum.from_file("boot/cabbage_boot{}.fs".format(i)) for i in range(100)]
    all_boot = [sfs.fold() for sfs in all_boot]

    # Now read in the original data set
    data = dadi.Spectrum.from_file("cabbage.fs")
    data = data.fold()
    pts_l = [100,110,120]
    func = expand_bottlegrowth_noF
    popt = [10.0000000000000000,1.9296649206007304,0.1031051605174129,0.1937327487823060,0.0100000000000000]

    func_ex = dadi.Numerics.make_extrap_log_func(func)

    # If we've already calculated uncertainties, there's no reason to do it again
    # Just assign to the uncerts variable from a previous run
    calc_uncerts = False
    if calc_uncerts:
        uncerts = dadi.Godambe.GIM_uncert(func_ex, pts_l, all_boot, popt, data,
                                          multinom=True)
    else:
        uncerts = [3.95715610e+00,5.61924059e-01,2.27036898e-03,2.20838193e-02,1.99685862e-03,6.02129250e+03]
    print("Estimated parameter standard deviations from GIM: {}".format(uncerts))

    # Set conversion parameters
    theta = 459187.876565839 # estimated from model
    #L = 554977060 # Sequence length
    L = 411560319 # Sequence length
    mu = 1.5e-8
    g  = 1 # generation time
    Nref = theta / L / mu / g / 4
    print("\n\nConverting parameters to actual units...\n")
    print("Using the following values for conversion:")
    print("  Sequence length = {}".format(L))
    print("  Mutation rate =   {}".format(mu))
    print("  Generation time = {}".format(g))

    print("\nNref = {} ({} -- {})".format(Nref, (theta-uncerts[-1]) / L / mu / g / 4, (theta+uncerts[-1]) / L / mu / g / 4))
    print("NExp = {} ({} -- {})".format(popt[0]*Nref, (popt[0]-uncerts[0])*Nref, (popt[0]+uncerts[0])*Nref))
    print("NBot = {} ({} -- {})".format(popt[1]*Nref, (popt[1]-uncerts[1])*Nref, (popt[1]+uncerts[1])*Nref))
    print("NCur = {} ({} -- {})".format(popt[2]*Nref, (popt[2]-uncerts[2])*Nref, (popt[2]+uncerts[2])*Nref))
    print("T1    = {} ({} -- {})".format(popt[3]*2*Nref, (popt[3]-uncerts[3])*2*Nref, (popt[3]+uncerts[3])*2*Nref))
    print("T2    = {} ({} -- {})".format(popt[4]*2*Nref, (popt[4]-uncerts[4])*2*Nref, (popt[4]+uncerts[4])*2*Nref))
