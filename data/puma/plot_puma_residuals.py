#!/usr/bin/env python3

import dadi
import matplotlib.pyplot as plt

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

def divergence_noF(params, ns, pts):
    """
    Divergence only model.
    """
    nuTX,nuFL,T = params

    xx = yy = dadi.Numerics.default_grid(pts)
    n1,n2 = ns[0],ns[1]

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx,phi)

    phi = dadi.Integration.two_pops(phi,xx,T,nu1=nuTX,nu2=nuFL)
    sfs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,yy))
    return sfs

if __name__ == "__main__":
    data = dadi.Spectrum.from_file("puma.fs")
    data = data.fold()
    pts_l = [40,50,60]
    func1 = divergence
    func1_ex = dadi.Numerics.make_extrap_log_func(func1)
    func2 = divergence_noF
    func2_ex = dadi.Numerics.make_extrap_log_func(func2)

    popt = [0.45043163,0.01252053,0.01002769,0.4814076365,0.6278800995]
    popt_noF = [0.11219582,0.01005386,0.01114349]

    model = func1_ex(popt, data.sample_sizes, pts_l)
    model = model.fold()
    model_noF = func2_ex(popt_noF, data.sample_sizes, pts_l)
    model_noF = model_noF.fold()


    dadi.Plotting.plot_2d_comp_multinom(model,data)
    #plt.savefig("puma_fit.pdf")
    #plt.close()

    dadi.Plotting.plot_2d_comp_multinom(model_noF, data)
    #plt.savefig("puma_fit_noF.pdf")
    #plt.close()
