#!/usr/bin/env python3

import dadi
import dadi.NLopt_mod
import nlopt

def three_epoch_noF(params, ns, pts):
    """
    params = (nuB,nuF,TB,TF)
    ns = (n1,)

    nuB: Ratio of bottleneck population size to ancient pop size
    nuF: Ratio of contemporary to ancient pop size
    TB: Length of bottleneck (in units of 2*Na generations)
    TF: Time since bottleneck recovery (in units of 2*Na generations)

    n1: Number of samples in resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nuB,nuF,TB,TF = params

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)

    phi = dadi.Integration.one_pop(phi, xx, TB, nuB)
    phi = dadi.Integration.one_pop(phi, xx, TF, nuF)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
    return fs

if __name__ == '__main__':
    # Read in
    data = dadi.Spectrum.from_file("cabbage.fs")
    data_f = data.fold()
    ns = data.sample_sizes
    pts_l = [100,110,120]
    func = three_epoch_noF
    upper_bound = [50.0,50.0,50.0,50.0]
    lower_bound = [1e-3,1e-3,1e-3,1e-3]
    p0 = [1.0,1.0,0.5,0.5]

    func_ex = dadi.Numerics.make_extrap_log_func(func)

    p0 = dadi.Misc.perturb_params(p0,fold=1,upper_bound=upper_bound,lower_bound=lower_bound)

    print("\n\n**** Starting NLopt optimization ****\n\n")
    popt,LLopt,result = dadi.NLopt_mod.opt(p0, data_f, func_ex, pts_l,
                                            lower_bound=lower_bound, upper_bound=upper_bound,
                                            algorithm=nlopt.LN_BOBYQA)
    model = func_ex(popt, ns, pts_l)

    model_f  = model.fold()

    ll_optF  = dadi.Inference.ll_multinom(model_f, data_f)

    thetaF  = dadi.Inference.optimal_sfs_scaling(model_f, data_f)

    with open("cabbage_fits_3epoch_noF.csv", 'a') as f1_out:
        for p in range(len(popt)):
            print("{}".format(popt[p]), ",", end='', sep='', file=f1_out)
        print(ll_optF, ",", thetaF, sep='', file=f1_out)
