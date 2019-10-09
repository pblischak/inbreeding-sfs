#!/usr/bin/env python3

import dadi
import matplotlib.pyplot as plt

import matplotlib
import pylab
import numpy

#: Custom ticks that label only the lowest and highest bins in an FS plot.
class _sfsTickLocator(matplotlib.ticker.Locator):
    def __call__(self):
        'Return the locations of the ticks'

        try:
            vmin, vmax = self.axis.get_view_interval()
            dmin, dmax = self.axis.get_data_interval()
        except AttributeError:
            self.verify_intervals()
            vmin, vmax = self.viewInterval.get_bounds()
            dmin, dmax = self.dataInterval.get_bounds()

        tmin = max(vmin, dmin)
        tmax = min(vmax, dmax)

        return numpy.array([round(tmin)+0.5, round(tmax)-0.5])
#: Custom tick formatter
_ctf = matplotlib.ticker.FuncFormatter(lambda x,pos: '%i' % (x-0.4))


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
        f = pylab.figure(fig_num, figsize=(7,7))
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
    pylab.ylim(-175,115)
    if plot_masked:
        pylab.plot(resid.data, '--og', mfc='w', zorder=-100)

    ax.set_xlim(0, data.shape[0]-1)
    if show:
        pylab.show()

def bottlegrowth(params, ns, pts):
    nuBot,nuCur,T,F = params

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)

    nu_func = lambda t: nuBot * (nuCur/nuBot) ** (t/T)
    phi = dadi.Integration.one_pop(phi,xx,T,nu=nu_func)
    sfs = dadi.Spectrum.from_phi_inbreeding(phi, ns, (xx,), (F,), (2,))
    return sfs

def expand_bottlegrowth(params,ns,pts):
    """
    This is the model that we used.
    """
    nuExp,nuBot,nuCur,T1,T2,F = params

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)

    phi = dadi.Integration.one_pop(phi,xx,T1,nu=nuExp)
    nu_func = lambda t: nuBot * (nuCur/nuBot) ** (t/T2)
    phi = dadi.Integration.one_pop(phi,xx,T2,nu=nu_func)
    sfs = dadi.Spectrum.from_phi_inbreeding(phi, ns, (xx,), (F,), (2,))
    return sfs

def bottlegrowth_noF(params, ns, pts):
    nuBot,nuCur,T = params

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

if __name__ == "__main__":
    data = dadi.Spectrum.from_file("cabbage.fs")
    data = data.fold()
    pts_l = [100,110,120]
    func1 = expand_bottlegrowth
    func1_ex = dadi.Numerics.make_extrap_log_func(func1)
    func2 = expand_bottlegrowth_noF
    func2_ex = dadi.Numerics.make_extrap_log_func(func2)

    popt = [1.8034444276865402,1.5406082627009385,7.128657409544192,0.4769670348935206,0.0143092989018097,0.577605971746217]
    popt_noF = [10.0000000000000000,1.9296649206007304,0.1031051605174129,0.1937327487823060,0.0100000000000000]

    model = func1_ex(popt, data.sample_sizes, pts_l)
    model = model.fold()
    model_noF = func2_ex(popt_noF, data.sample_sizes, pts_l)
    model_noF = model_noF.fold()

    #dadi.Plotting.plot_1d_comp_multinom(model,data)
    plot_1d_comp_multinom(model,data)
    #plt.savefig("puma_fit.pdf")
    #plt.close()

    #dadi.Plotting.plot_1d_comp_multinom(model_noF, data)
    plot_1d_comp_multinom(model_noF, data)
    #plt.savefig("puma_fit_noF.pdf")
    #plt.close()
