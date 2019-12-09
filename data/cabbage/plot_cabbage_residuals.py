#!/usr/bin/env python3

import dadi
import matplotlib
import pylab
import numpy
from dadi import Numerics, Inference

def plot_1d_comp_multinom(model, data, fig_num=None, residual='Anscombe',
                          plot_masked=False):
    """
    Mulitnomial comparison between 1d model and data.


    model: 1-dimensional model SFS
    data: 1-dimensional data SFS
    fig_num: Clear and use figure fig_num for display. If None, an new figure
             window is created.
    residual: 'Anscombe' for Anscombe residuals, which are more normally
              distributed for Poisson sampling. 'linear' for the linear
              residuals, which can be less biased.
    plot_masked: Additionally plots (in open circles) results for points in the
                 model or data that were masked.

    This comparison is multinomial in that it rescales the model to optimally
    fit the data.
    """
    model = Inference.optimally_scaled_sfs(model, data)

    plot_1d_comp_Poisson(model, data, fig_num, residual,
                         plot_masked)

def plot_1d_comp_Poisson(model, data, fig_num=None, residual='Anscombe',
                         plot_masked=False, show=True):
    """
    Poisson comparison between 1d model and data.


    model: 1-dimensional model SFS
    data: 1-dimensional data SFS
    fig_num: Clear and use figure fig_num for display. If None, an new figure
             window is created.
    residual: 'Anscombe' for Anscombe residuals, which are more normally
              distributed for Poisson sampling. 'linear' for the linear
              residuals, which can be less biased.
    plot_masked: Additionally plots (in open circles) results for points in the
                 model or data that were masked.
    show: If True, execute pylab.show command to make sure plot displays.
    """
    if fig_num is None:
        f = pylab.gcf()
    else:
        f = pylab.figure(fig_num, figsize=(10,8))
    pylab.clf()

    if data.folded and not model.folded:
        model = model.fold()

    masked_model, masked_data = Numerics.intersect_masks(model, data)

    ax = pylab.subplot(2,1,1)
    pylab.semilogy(masked_data, '-ob')
    pylab.semilogy(masked_model, '-or')

    if plot_masked:
        pylab.semilogy(masked_data.data, '--ob', mfc='w', zorder=-100)
        pylab.semilogy(masked_model.data, '--or', mfc='w', zorder=-100)

    pylab.subplot(2,1,2, sharex = ax)
    if residual == 'Anscombe':
        resid = Inference.Anscombe_Poisson_residual(masked_model, masked_data)
    elif residual == 'linear':
        resid = Inference.linear_Poisson_residual(masked_model, masked_data)
    else:
        raise ValueError("Unknown class of residual '%s'." % residual)
    pylab.plot(resid, '-og')
    pylab.ylim(-160,120)
    if plot_masked:
        pylab.plot(resid.data, '--og', mfc='w', zorder=-100)

    ax.set_xlim(0, data.shape[0]-1)
    if show:
        pylab.show()

def three_epoch(params, ns, pts):
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
    nuB,nuF,TB,TF,F = params

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)

    phi = dadi.Integration.one_pop(phi, xx, TB, nuB)
    phi = dadi.Integration.one_pop(phi, xx, TF, nuF)

    fs = dadi.Spectrum.from_phi_inbreeding(phi, ns, (xx,), (F,), (2,))
    return fs

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

if __name__ == "__main__":
    data = dadi.Spectrum.from_file("cabbage.fs")
    data = data.fold()
    pts_l = [100,110,120]
    func1 = three_epoch
    func1_ex = dadi.Numerics.make_extrap_log_func(func1)
    func2 = three_epoch_noF
    func2_ex = dadi.Numerics.make_extrap_log_func(func2)

    popt = [1.810449088130342,12.2790194725467110,0.47393521119534737,0.00921096365957015,0.577870722504928]
    popt_noF = [6.4524615672350958,0.0309347139612217,0.153264805381591,0.00100000000000000]

    model = func1_ex(popt, data.sample_sizes, pts_l)
    model = model.fold()
    model_noF = func2_ex(popt_noF, data.sample_sizes, pts_l)
    model_noF = model_noF.fold()

    plot_1d_comp_multinom(model,data, fig_num=1)
    #plt.savefig("puma_fit.pdf")
    #plt.close()

    plot_1d_comp_multinom(model_noF, data, fig_num=2)
    #plt.savefig("puma_fit_noF.pdf")
    #plt.close()
