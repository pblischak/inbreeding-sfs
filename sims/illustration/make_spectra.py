#!/usr/bin/env python3

import dadi

def snm_inbreeding(params, ns, pts):
    """
    Model 1
    -------

    Standard neutral model with inbreeding.

    Parameters
    ----------

    F: inbreeding coefficient (0 < F < 1).
    """
    F = params[0]
    # double-check that F is in range
    if F <= 0.0:
        F = 1e-4
    elif F >= 1.0:
        F = 1.0 - 1e-4
    else:
        pass

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)

    fs = dadi.Spectrum.from_phi_inbreeding(phi, ns, (xx,), (F,), (2,))
    return fs

if __name__ == "__main__":

    pts_l = [70,80,90]
    func = snm_inbreeding
    func_ex = dadi.Numerics.make_extrap_log_func(func)
    fs_F00 = func_ex([1e-6], [50], pts_l)
    fs_F10 = func_ex([0.1], [50], pts_l)
    fs_F25 = func_ex([0.25], [50], pts_l)
    fs_F50 = func_ex([0.5], [50], pts_l)
    fs_F75 = func_ex([0.75], [50], pts_l)
    fs_F90 = func_ex([0.9], [50], pts_l)

    with open("spectra.csv", "w") as f_out:  
        print("Fis","X","Y", sep=",", file=f_out)
        for i in range(1,50):
            print("0.0",i,fs_F00[i], sep=",", file=f_out)
            print("0.1",i,fs_F10[i], sep=",", file=f_out)
            print("0.25",i,fs_F25[i], sep=",", file=f_out)
            print("0.5",i,fs_F50[i], sep=",", file=f_out)
            print("0.75",i,fs_F75[i], sep=",", file=f_out)
            print("0.9",i,fs_F90[i], sep=",", file=f_out)
