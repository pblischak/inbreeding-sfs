#!/usr/bin/env python3

import dadi
import matplotlib.pyplot as plt

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

if __name__ == "__main__":
    data = dadi.Spectrum.from_file("cabbage.fs")
    data = data.fold()
    pts_l = [100,110,120]
    func1 = expand_bottlegrowth
    func1_ex = dadi.Numerics.make_extrap_log_func(func1)
    func2 = expand_bottlegrowth_noF
    func2_ex = dadi.Numerics.make_extrap_log_func(func2)

    popt = [1.8034444276865402,1.5406082627009385,7.128657409544192,0.4769670348935206,0.0143092989018097,0.577605971746217]
    popt_noF = [10.0000000000000000,1.9296649206007304,0.1031051605174129,0.1937327487823060,0.0100000000000000]

    model = func1_ex(popt, data.sample_sizes, pts_l)
    model = model.fold()
    model_noF = func2_ex(popt_noF, data.sample_sizes, pts_l)
    model_noF = model_noF.fold()

    dadi.Plotting.plot_1d_comp_multinom(model,data)
    #plt.savefig("puma_fit.pdf")
    #plt.close()

    dadi.Plotting.plot_1d_comp_multinom(model_noF, data)
    #plt.savefig("puma_fit_noF.pdf")
    #plt.close()
