#!/usr/bin/env python3

import dadi
import numpy as np
from sys import argv
from models_1D import model2, model2_noF

#F_array   = [0.05, 0.25, 0.50, 0.75, 0.95]
F = float(argv[1])
nu0_array = [0.1, 0.25, 0.50]
T_array   = [0.1, 0.2, 0.3]

# Simulation rep settings
nsims = 20
reps_per_sim = 50
f_out = open("1D_model2_noF_{0}.out".format(argv[1]), 'w')
print("Sim", "F_true", "Nu0_true", "Nu0_est", "T_true", "T_est", "Theta", "Loglik",  sep=",", file=f_out)
f_out.flush()

# Set up dadi optimization
N = 50
L = 1e4
pts_l = [70,80,90]
upper_bound = [10.0, 10.0]
lower_bound = [1.0e-3, 1.0e-3]
func = model2_noF
func_ex = dadi.Numerics.make_extrap_log_func(func)

func2 = model2
func2_ex = dadi.Numerics.make_extrap_log_func(func2)

for nu in nu0_array:
    for tt in T_array:
        for s in range(nsims):
            res = {}
            for r in range(reps_per_sim):
                print(F, nu, tt, s+1, r+1, sep='\t')
                tmp = func2_ex([nu,tt,F], [N], pts_l)
                data = tmp * L
                data = data.sample()
                p0 = dadi.Misc.perturb_params([nu,tt], fold=1, upper_bound=upper_bound, lower_bound=lower_bound)
                popt = dadi.Inference.optimize_log(p0, data, func_ex, pts_l,
                                                   lower_bound=lower_bound,
                                                   upper_bound=upper_bound,
                                                   verbose=len(p0), maxiter=3)
                model = func_ex(popt, [N], pts_l)
                ll_model = dadi.Inference.ll_multinom(model, data)
                theta = dadi.Inference.optimal_sfs_scaling(model, data)
                res[ll_model] = [popt[0],popt[1],theta]
            ll_max = np.max(list(res))
            print(s+1, F, nu, res[ll_max][0], tt, res[ll_max][1], res[ll_max][2], ll_max, sep=",", file=f_out)
            f_out.flush()
