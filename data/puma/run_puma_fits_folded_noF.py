#!/usr/bin/env python3

import dadi

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


if __name__ == '__main__':
    # Read in data
    data = dadi.Spectrum.from_file("puma.fs")
    data = data.fold()
    ns = data.sample_sizes
    pts_l = [40,50,60]
    func = divergence_noF
    upper_bound = [10.0,10.0,10.0]
    lower_bound = [1e-2,1e-2,1e-2]
    p0 = [1.0,1.0,1.0]

    func_ex = dadi.Numerics.make_extrap_log_func(func)

    p0 = dadi.Misc.perturb_params(p0,fold=1,upper_bound=upper_bound,lower_bound=lower_bound)

    popt = dadi.Inference.optimize_log(p0, data, func_ex, pts_l,
                                   lower_bound=lower_bound,
                                   upper_bound=upper_bound,
                                   verbose=len(p0))
    model = func_ex(popt, ns, pts_l)
    ll_model = dadi.Inference.ll_multinom(model, data)
    theta = dadi.Inference.optimal_sfs_scaling(model, data)

    with open("puma_fits_folded_noF.csv", 'a') as f:
        for p in popt:
            print("{},".format(p), sep='',end='', file=f)
        print(ll_model,theta, sep=",", file=f)