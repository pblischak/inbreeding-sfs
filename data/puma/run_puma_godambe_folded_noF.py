#!/usr/bin/env python3

import dadi
import dadi.Godambe

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
    # Make SFS files so we don't have to read in
    # big data files every time.
    make_spectra = False # Swith this to False after SFS files have been generated
    for i in range(100):
        print(i)
        if make_spectra:
            boot_dd = dadi.Misc.make_data_dict("puma_boot{}.dadi".format(i))
            boot_data = dadi.Spectrum.from_data_dict(boot_dd, ['Texas','Florida'], [10,4])
            boot_data.to_file("puma_boot{}.fs".format(i))

    # Read in the bootstrapped spectra and then fold
    all_boot = [dadi.Spectrum.from_file("puma_boot{}.fs".format(i)) for i in range(100)]
    all_boot = [sfs.fold() for sfs in all_boot]

    # Now read in the original data set
    data = dadi.Spectrum.from_file("puma.fs")
    data = data.fold()
    pts_l = [40,50,60]
    func = divergence_noF
    popt = [0.11219582,0.01005386,0.01114349]

    func_ex = dadi.Numerics.make_extrap_log_func(func)

    uncerts = dadi.Godambe.GIM_uncert(func_ex, pts_l, all_boot, popt, data,
                                  multinom=True)
    print('Estimated parameter standard deviations from GIM: {0}'.format(uncerts))
