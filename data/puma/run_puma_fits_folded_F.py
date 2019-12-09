#!/usr/bin/env python3

import dadi

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
    # Read in
    dd = dadi.Misc.make_data_dict("puma.dadi")
    data = dadi.Spectrum.from_data_dict(dd, ['Texas','Florida'], [10,4])
    data = data.fold()
    ns = data.sample_sizes
    pts_l = [40,50,60]
    func = divergence
    upper_bound = [10.0,10.0,10.0,0.9999,0.9999]
    lower_bound = [1e-2,1e-2,1e-2,0.0001,0.0001]
    p0 = [1.0,1.0,1.0,0.1,0.1]

    func_ex = dadi.Numerics.make_extrap_log_func(func)

    p0 = dadi.Misc.perturb_params(p0,fold=1,upper_bound=upper_bound,lower_bound=lower_bound)

    popt = dadi.Inference.optimize_log(p0, data, func_ex, pts_l,
                                   lower_bound=lower_bound,
                                   upper_bound=upper_bound,
                                   verbose=len(p0))
    model = func_ex(popt, ns, pts_l)
    ll_model = dadi.Inference.ll_multinom(model, data)
    theta = dadi.Inference.optimal_sfs_scaling(model, data)

    with open("puma_fits_folded.csv", 'a') as f:
        for p in popt:
            print("{},".format(p), sep='',end='', file=f)
        print(ll_model,theta, sep=",", file=f)
