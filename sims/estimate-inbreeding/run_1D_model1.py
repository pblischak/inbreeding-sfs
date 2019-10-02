#!/usr/bin/env python3

import dadi
import numpy as np
from sys import argv
from models_1D import model1

#F_array = [0.1, 0.25, 0.5, 0.75, 0.9]
F = float(argv[1])

# Simulation rep settings
nsims = 20
reps_per_sim = 50
f_out = open("1D_model1_{0}.out".format(argv[1]), 'w')
print("Sim", "F_true", "F_est", "Theta", "Loglik", sep=",", file=f_out)
f_out.flush()

# Set up dadi optimization
N = 50
L = 1e4
pts_l = [70,80,90]
upper_bound = [0.9999]
lower_bound = [1.0-0.9999]
func = model1
func_ex = dadi.Numerics.make_extrap_log_func(func)

for s in range(nsims):
    res = {}
    for r in range(reps_per_sim):
        print(F, s+1, r+1, sep='\t')
        tmp = func_ex([F], [N], pts_l)
        data = tmp * L
        data = data.sample()
        p0 = dadi.Misc.perturb_params([F], fold=1, upper_bound=upper_bound, lower_bound=lower_bound)
        popt = dadi.Inference.optimize_log(p0, data, func_ex, pts_l,
                                           lower_bound=lower_bound,
                                           upper_bound=upper_bound,
                                           verbose=len(p0), maxiter=3)
        model = func_ex(popt, [N], pts_l)
        ll_model = dadi.Inference.ll_multinom(model, data)
        theta = dadi.Inference.optimal_sfs_scaling(model, data)
        res[ll_model] = [popt[0],theta]
    ll_max = np.max(list(res))
    print(s+1, F, res[ll_max][0], res[ll_max][1], ll_max, sep=",", file=f_out)
    f_out.flush()
