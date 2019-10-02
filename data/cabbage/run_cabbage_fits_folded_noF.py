#!/usr/bin/env python3

import dadi
import dadi.NLopt_mod
import nlopt

def bottlegrowth(params, ns, pts):
    nuBot,nuCur,T,F = params

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
    # Read in
    data = dadi.Spectrum.from_file("cabbage.fs")
    data_f = data.fold()
    ns = data.sample_sizes
    pts_l = [100,110,120]
    func = expand_bottlegrowth_noF
    upper_bound = [10.0,10.0,10.0,10.0,10.0]
    lower_bound = [1e-2,1e-2,1e-2,1e-2,1e-2]
    p0 = [2.0,0.5,1.0,1.0,1.0]

    func_ex = dadi.Numerics.make_extrap_log_func(func)

    p0 = dadi.Misc.perturb_params(p0,fold=1,upper_bound=upper_bound,lower_bound=lower_bound)

    popt = dadi.Inference.optimize_log(p0, data_f, func_ex, pts_l,
                                   lower_bound=lower_bound,
                                   upper_bound=upper_bound,
                                   verbose=len(p0))

    print("\n\n**** Starting NLopt optimization ****\n\n")
    popt2,LLopt,result = dadi.NLopt_mod.opt(p0, data_f, func_ex, pts_l,
                                            lower_bound=lower_bound, upper_bound=upper_bound,
                                            algorithm=nlopt.LN_BOBYQA)
    model = func_ex(popt, ns, pts_l)
    model2 = func_ex(popt2, ns, pts_l)

    model_f  = model.fold()
    model2_f = model2.fold()

    ll_optF  = dadi.Inference.ll_multinom(model_f, data_f)
    ll2_optF = dadi.Inference.ll_multinom(model2_f, data_f)

    thetaF  = dadi.Inference.optimal_sfs_scaling(model_f, data_f)
    theta2F = dadi.Inference.optimal_sfs_scaling(model2_f, data_f)

    print("\n\n  Is there a likelihood difference: {}\n\n".format(ll2_optF - LLopt))
    with open("cabbage_fits_noF_oldOpt.csv", 'a') as f1_out, open("cabbage_fits_noF_newOpt.csv", 'a') as f2_out:
        for p in range(len(popt)):
            print("{}".format(popt[p]), ",", end='', sep='', file=f1_out)
            print("{}".format(popt2[p]), ",", end='', sep='', file=f2_out)
        print(ll_optF, ",", thetaF, sep='', file=f1_out)
        print(LLopt, ",", theta2F, sep='', file=f2_out)
