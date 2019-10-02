"""
2D models for testing inference with inbreeding.
"""

import dadi
import numpy as np

def model1(params, ns, pts):
    """
    Model 1
    -------

    A population splits from the ancestral population, pop1 stays the same size
    and pop2 is nu2 percent of the original size. The split happens at time T,
    and there is migration from population 1 into population 2 at a rate of m21.
    Pop1 is also inbred, with coefficient equal to F.
    """
    nu2, T, m21, F = params
    if F <= 0.0:
        F = 1e-4
    elif F >= 1.0:
        F = 1.0 - 1e-4
    else:
        pass

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    phi = dadi.Integration.two_pops(phi, xx, T, 1.0, nu2, m12=0.0, m21=m21)

    fs = dadi.Spectrum.from_phi_inbreeding(phi, ns, (xx,xx), (0.0001,F), (2,2))
    return fs

def model1_noF(params, ns, pts):
    """
    Model 1
    -------

    A population splits from the ancestral population, pop1 stays the same size
    and pop2 is nu2 percent of the original size. The split happens at time T,
    and there is migration from population 1 into population 2 at a rate of m21.
    """
    nu2, T, m21 = params

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    phi = dadi.Integration.two_pops(phi, xx, T, 1.0, nu2, m12=0.0, m21=m21)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs
