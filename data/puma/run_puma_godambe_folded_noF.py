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
    all_boot = [dadi.Spectrum.from_file("boot/puma_test2_boot{}.fs".format(i)) for i in range(100)]
    all_boot = [sfs.fold() for sfs in all_boot]

    # Now read in the original data set
    data = dadi.Spectrum.from_file("puma.fs")
    data = data.fold()
    pts_l = [40,50,60]
    func = divergence_noF
    popt = [0.197150919878390,0.0100564473768224,0.0371335637849756,0.0113941805000981]

    func_ex = dadi.Numerics.make_extrap_log_func(func)

    # Calculate parameter uncertainties
    """
    We include this section to check uncertainties across different step sizes (eps).
    This can be turned off by setting check_grid_size=False.
    """
    check_grid_sizes = False
    if check_grid_sizes:
        print("\nChecking uncertainties with different grid sizes:")
        for e in [1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7]:
    	    u = dadi.Godambe.GIM_uncert(func_ex, pts_l, all_boot, popt, data, log=True, multinom=True, eps=e)
    	    print("{} = {}".format(e,u))
        exit(0)
    
    eps=1e-2 # the default is 1e-2
    uncerts,GIM = dadi.Godambe.GIM_uncert(func_ex, pts_l, all_boot, popt, data, log=True, multinom=True, return_GIM=True,eps=eps)
    vcov = np.linalg.inv(GIM)
    
    # Set conversion parameters
    theta = 2718260.477394347 # estimated from model
    L = 2564692624
    mu = 2.2e-9
    g  = 3 # generation time
    scalar = L*mu*4
    Nref = theta / L / mu / 4
    scalar = L*mu*4
    print("\n\nConverting parameters to actual units (eps={})...\n".format(eps))
    print("Using the following values for conversion:")
    print("  Sequence length = {}".format(L))
    print("  Mutation rate   = {}".format(mu))
    print("  Generation time = {}".format(g))
    print("")
    
    uncerts2 = [
        np.sqrt(vcov[-1,-1] + vcov[0,0] + 2*vcov[0,-1]),
        np.sqrt(vcov[-1,-1] + vcov[1,1] + 2*vcov[1,-1]),
        np.sqrt(vcov[-1,-1] + vcov[2,2] + 2*vcov[2,-1]),
        np.sqrt(vcov[-1,-1] + vcov[3,3] + 2*vcov[3,-1]),
        uncerts[-1]
    ]
    
    log_params = [
        np.log(theta) + np.log(popt[0]) + np.log(1/scalar),   # N_TX
        np.log(theta) + np.log(popt[1]) + np.log(1/scalar),   # N_FL
        np.log(theta) + np.log(popt[2]) + np.log(2*g/scalar), # T1
        np.log(theta) + np.log(popt[3]) + np.log(2*g/scalar), # T2
    ]
    
    print("Estimated parameter standard deviations from GIM:\n{}\n".format(uncerts))
    print("Estimated parameter standard deviations from error propagation:\n{}\n".format(uncerts2))
    print("Variance-Covariance Matrix:\n{}\n".format(vcov))

    """
    With propogation of uncertainty
    """
    print("\nParameter estimates and 95% confidence intervals:")
    print("Nref = {} ({}--{})".format(Nref, np.exp(np.log(theta)-1.96*uncerts2[-1])/scalar, np.exp(np.log(theta)+1.96*uncerts2[-1])/scalar))
    print("N_TX = {} ({}--{})".format(popt[0]*Nref, np.exp(log_params[0]-1.96*uncerts2[0]), np.exp(log_params[0]+1.96*uncerts2[0])))
    print("N_FL = {} ({}--{})".format(popt[1]*Nref, np.exp(log_params[1]-1.96*uncerts2[1]), np.exp(log_params[1]+1.96*uncerts2[1])))
    print("T1   = {} ({}--{})".format(popt[2]*2*g*Nref, np.exp(log_params[2]-1.96*uncerts2[2]), np.exp(log_params[2]+1.96*uncerts2[2])))
    print("T2   = {} ({}--{})".format(popt[3]*2*g*Nref, np.exp(log_params[3]-1.96*uncerts2[3]), np.exp(log_params[3]+1.96*uncerts2[3])))
