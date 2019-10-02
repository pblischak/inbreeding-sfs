#!/usr/bin/env python3

import dadi
import numpy as np
from sys import argv
from models_2D import model1, model1_noF

#F_array   = [0.05, 0.25, 0.50, 0.75, 0.95]
F = float(argv[1])
nu2_array = [0.1, 0.25, 0.5]
T_array   = [0.1, 0.2, 0.3]
m21_array = [0.5, 1.0, 1.5]

# Simulation rep settings
nsims = 20
reps_per_sim = 50
f_out = open("2D_model1_{0}.out".format(argv[1]), 'w')
print("Sim", "F_true", "Nu2_true", "Nu2_est", "T_true", "T_est", "M21_true", "M21_est", "Theta", "Loglik",  sep=",", file=f_out)
f_out.flush()

# Set up dadi optimization
N = 50
L = 1e4
pts_l = [70,80,90]
upper_bound = [10.0, 10.0, 10.0,0.9999]
lower_bound = [1.0e-3, 1.0e-3, 1.0e-3,1.0-0.9999]
func = model1
func_ex = dadi.Numerics.make_extrap_log_func(func)
func2 = model1_noF
func2_ex = dadi.Numerics.make_extrap_log_func(func2)

for nu in nu2_array:
    for tt in T_array:
        for mm in m21_array:
            for s in range(nsims):
                res = {}
                for r in range(reps_per_sim):
                    print(F,nu,tt,mm,s+1,r+1,sep='\t')
                    tmp = func2_ex([nu,tt,mm], [N,N], pts_l)
                    data = tmp * L
                    data = data.sample()
                    p0 = dadi.Misc.perturb_params([nu,tt,mm,F], fold=1, upper_bound=upper_bound, lower_bound=lower_bound)
                    popt = dadi.Inference.optimize_log(p0, data, func_ex, pts_l,
                                                       lower_bound=lower_bound,
                                                       upper_bound=upper_bound,
                                                       verbose=len(p0), maxiter=3)
                    model = func_ex(popt, data.sample_sizes, pts_l)
                    ll_model = dadi.Inference.ll_multinom(model, data)
                    theta = dadi.Inference.optimal_sfs_scaling(model, data)
                    res[ll_model] = [popt[0],popt[1],popt[2],popt[3],theta]
                ll_max = np.max(list(res))
                print(s+1, F, res[ll_max][3], nu, res[ll_max][0], tt, res[ll_max][1], mm, res[ll_max][2], res[ll_max][4], ll_max, sep=",", file=f_out)
                f_out.flush()
