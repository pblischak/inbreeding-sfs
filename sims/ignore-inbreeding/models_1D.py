"""
1D models for testing inference with inbreeding.
"""

import dadi
import numpy as np

def model1(F, ns, pts):
    """
    Model 1
    -------

    Standard neutral model with inbreeding.

    Parameters
    ----------

    F: inbreeding coefficient (0 < F < 1).
    """
    # check that F is in range
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

def model2(params, ns, pts):
    """
    Model 2
    -------

    Bottleneck followed by growth.

    Parameters
    ----------

    nu0: Relative size of pop after bottleneck.
    T:   Time of bottleneck.
    F:   Inbreeding coefficient (0 < F < 1).
    """
    nu0, T, F = params

    # check that F is in range
    if F <= 0.0:
        F = 1e-4
    elif F >= 1.0:
        F = 1.0 - 1e-4
    else:
        pass

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)

    nu_func = lambda t: nu0 * (1.0 / nu0) ** (t/T)
    phi = dadi.Integration.one_pop(phi, xx, T, nu_func)

    fs = dadi.Spectrum.from_phi_inbreeding(phi, ns, (xx,), (F,), (2,))
    return fs

def model2_noF(params, ns, pts):
    """
    Model 2
    -------

    Bottleneck followed by growth.

    Parameters
    ----------

    nu0: Relative size of pop after bottleneck.
    T:   Time of bottleneck.
    """
    nu0, T = params

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)

    nu_func = lambda t: nu0 * (1.0 / nu0) ** (t/T)
    phi = dadi.Integration.one_pop(phi, xx, T, nu_func)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
    return fs
