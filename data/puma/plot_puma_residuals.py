#!/usr/bin/env python3

import matplotlib
import pylab
import numpy
import dadi
from dadi import Numerics, Inference

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

def plot_single_2d_sfs(sfs, vmin=None, vmax=None, ax=None,
                       pop_ids=None, extend='neither', colorbar=True,
                       cmap=pylab.cm.viridis_r):
    """
    Heatmap of single 2d SFS.

    If vmax is greater than a factor of 10, plot on log scale.

    Returns colorbar that is created.

    sfs: SFS to plot
    vmin: Values in sfs below vmin are masked in plot.
    vmax: Values in sfs above vmax saturate the color spectrum.
    ax: Axes object to plot into. If None, the result of pylab.gca() is used.
    pop_ids: If not None, override pop_ids stored in Spectrum.
    extend: Whether the colorbar should have 'extension' arrows. See
            help(pylab.colorbar) for more details.
    colorbar: Should we plot a colorbar?
    cmap: Pylab colormap to use for plotting.
    """
    if ax is None:
        ax = pylab.gca()

    if vmin is None:
        vmin = sfs.min()
    if vmax is None:
        vmax = sfs.max()

    pylab.cm.hsv.set_under('w')
    if vmax / vmin > 10:
        # Under matplotlib 1.0.1, default LogFormatter omits some tick lines.
        # This works more consistently.
        norm = matplotlib.colors.LogNorm(vmin=vmin*(1-1e-3), vmax=vmax*(1+1e-3))
        format = matplotlib.ticker.LogFormatterMathtext()
    else:
        norm = matplotlib.colors.Normalize(vmin=vmin*(1-1e-3),
                                           vmax=vmax*(1+1e-3))
        format = None
    mappable=ax.pcolor(numpy.ma.masked_where(sfs<vmin, sfs),
                       cmap=cmap, edgecolors='none',
                       norm=norm)
    cb = ax.figure.colorbar(mappable, extend=extend, format=format)
    if not colorbar:
        ax.figure.delaxes(ax.figure.axes[-1])
    else:
        # A hack so we can manually work around weird ticks in some colorbars
        try:
            ax.figure.dadi_colorbars.append(cb)
        except AttributeError:
            ax.figure.dadi_colorbars = [cb]

    ax.plot([0,sfs.shape[1]],[0, sfs.shape[0]], '-k', lw=0.2)

    if pop_ids is None:
        if sfs.pop_ids is not None:
            pop_ids = sfs.pop_ids
        else:
            pop_ids = ['pop0','pop1']
    ax.set_ylabel(pop_ids[0], verticalalignment='top')
    ax.set_xlabel(pop_ids[1], verticalalignment='bottom')

    ax.xaxis.set_major_formatter(_ctf)
    ax.xaxis.set_major_locator(_sfsTickLocator())
    ax.yaxis.set_major_formatter(_ctf)
    ax.yaxis.set_major_locator(_sfsTickLocator())
    for tick in ax.xaxis.get_ticklines() + ax.yaxis.get_ticklines():
        tick.set_visible(False)

    ax.set_xlim(0, sfs.shape[1])
    ax.set_ylim(0, sfs.shape[0])

    return cb

def plot_2d_resid(resid, resid_range=None, ax=None, pop_ids=None,
                  extend='neither', colorbar=True,cmap=pylab.cm.RdBu_r):
    """
    Linear heatmap of 2d residual array.

    sfs: Residual array to plot.
    resid_range: Values > resid range or < resid_range saturate the color
                 spectrum.
    ax: Axes object to plot into. If None, the result of pylab.gca() is used.
    pop_ids: If not None, override pop_ids stored in Spectrum.
    extend: Whether the colorbar should have 'extension' arrows. See
            help(pylab.colorbar) for more details.
    colorbar: Should we plot a colorbar?
    """
    if ax is None:
        ax = pylab.gca()

    if resid_range is None:
        resid_range = abs(resid).max()

    mappable=ax.pcolor(resid, cmap=cmap, vmin=-resid_range,
                       vmax=resid_range, edgecolors='none')

    cbticks = [-resid_range, 0, resid_range]
    format = matplotlib.ticker.FormatStrFormatter('%.2g')
    cb = ax.figure.colorbar(mappable, ticks=cbticks, format=format,
                            extend=extend)
    if not colorbar:
        ax.figure.delaxes(ax.figure.axes[-1])
    else:
        try:
            ax.figure.dadi_colorbars.append(cb)
        except AttributeError:
            ax.figure.dadi_colorbars = [cb]

    ax.plot([0,resid.shape[1]],[0, resid.shape[0]], '-k', lw=0.2)

    if pop_ids is None:
        if resid.pop_ids is not None:
            pop_ids = resid.pop_ids
        else:
            pop_ids = ['pop0','pop1']
    ax.set_ylabel(pop_ids[0], verticalalignment='top')
    ax.set_xlabel(pop_ids[1], verticalalignment='bottom')

    ax.xaxis.set_major_formatter(_ctf)
    ax.xaxis.set_major_locator(_sfsTickLocator())
    ax.yaxis.set_major_formatter(_ctf)
    ax.yaxis.set_major_locator(_sfsTickLocator())
    for tick in ax.xaxis.get_ticklines() + ax.yaxis.get_ticklines():
        tick.set_visible(False)

    ax.set_xlim(0, resid.shape[1])
    ax.set_ylim(0, resid.shape[0])

    return cb

# Used to determine whether colorbars should have 'extended' arrows
_extend_mapping = {(True, True): 'neither',
                   (False, True): 'min',
                   (True, False): 'max',
                   (False, False): 'both'}

def plot_2d_comp_multinom(model, data, vmin=None, vmax=None,
                          resid_range=None, fig_num=None,
                          pop_ids=None, residual='Anscombe',
                          adjust=True,show=True):
    """
    Mulitnomial comparison between 2d model and data.


    model: 2-dimensional model SFS
    data: 2-dimensional data SFS
    vmin, vmax: Minimum and maximum values plotted for sfs are vmin and
                vmax respectively.
    resid_range: Residual plot saturates at +- resid_range.
    fig_num: Clear and use figure fig_num for display. If None, an new figure
             window is created.
    pop_ids: If not None, override pop_ids stored in Spectrum.
    residual: 'Anscombe' for Anscombe residuals, which are more normally
              distributed for Poisson sampling. 'linear' for the linear
              residuals, which can be less biased.
    adjust: Should method use automatic 'subplots_adjust'? For advanced
            manipulation of plots, it may be useful to make this False.
    show: Display the figure? False is useful for saving many comparisons
            in a loop.

    This comparison is multinomial in that it rescales the model to optimally
    fit the data.
    """
    model = Inference.optimally_scaled_sfs(model, data)

    plot_2d_comp_Poisson(model, data, vmin=vmin, vmax=vmax,
                         resid_range=resid_range, fig_num=fig_num,
                         pop_ids=pop_ids, residual=residual,
                         adjust=adjust,show=show)

def plot_2d_comp_Poisson(model, data, vmin=None, vmax=None,
                         resid_range=None, fig_num=None,
                         pop_ids=None, residual='Anscombe',
                         adjust=True, show=True):
    """
    Poisson comparison between 2d model and data.


    model: 2-dimensional model SFS
    data: 2-dimensional data SFS
    vmin, vmax: Minimum and maximum values plotted for sfs are vmin and
                vmax respectively.
    resid_range: Residual plot saturates at +- resid_range.
    fig_num: Clear and use figure fig_num for display. If None, an new figure
             window is created.
    pop_ids: If not None, override pop_ids stored in Spectrum.
    residual: 'Anscombe' for Anscombe residuals, which are more normally
              distributed for Poisson sampling. 'linear' for the linear
              residuals, which can be less biased.
    adjust: Should method use automatic 'subplots_adjust'? For advanced
            manipulation of plots, it may be useful to make this False.
    """
    if data.folded and not model.folded:
        model = model.fold()

    masked_model, masked_data = Numerics.intersect_masks(model, data)

    if fig_num is None:
        f = pylab.gcf()
    else:
        f = pylab.figure(fig_num, figsize=(10,8))

    pylab.clf()
    if adjust:
        pylab.subplots_adjust(bottom=0.07, left=0.07, top=0.94, right=0.95,
                              hspace=0.26, wspace=0.26)

    max_toplot = max(masked_model.max(), masked_data.max())
    min_toplot = min(masked_model.min(), masked_data.min())
    if vmax is None:
        vmax = max_toplot
    if vmin is None:
        vmin = min_toplot
    extend = _extend_mapping[vmin <= min_toplot, vmax >= max_toplot]

    if pop_ids is not None:
        data_pop_ids = model_pop_ids = resid_pop_ids = pop_ids
        if len(pop_ids) != 2:
            raise ValueError('pop_ids must be of length 2.')
    else:
        data_pop_ids = masked_data.pop_ids
        model_pop_ids = masked_model.pop_ids
        if masked_model.pop_ids is None:
            model_pop_ids = data_pop_ids

        if model_pop_ids == data_pop_ids:
           resid_pop_ids = model_pop_ids
        else:
            resid_pop_ids = None

    ax = pylab.subplot(2,2,1)
    plot_single_2d_sfs(masked_data, vmin=vmin, vmax=vmax,
                       pop_ids=data_pop_ids, colorbar=False)
    ax.set_title('data')

    ax2 = pylab.subplot(2,2,2, sharex=ax, sharey=ax)
    plot_single_2d_sfs(masked_model, vmin=vmin, vmax=vmax,
                       pop_ids=model_pop_ids, extend=extend)
    ax2.set_title('model')

    if residual == 'Anscombe':
        resid = Inference.Anscombe_Poisson_residual(masked_model, masked_data,
                                              mask=vmin)
    elif residual == 'linear':
        resid = Inference.linear_Poisson_residual(masked_model, masked_data,
                                            mask=vmin)
    else:
        raise ValueError("Unknown class of residual '%s'." % residual)

    if resid_range is None:
        resid_range = max((abs(resid.max()), abs(resid.min())))
    resid_extend = _extend_mapping[-resid_range <= resid.min(),
                                   resid_range >= resid.max()]

    ax3 = pylab.subplot(2,2,3, sharex=ax, sharey=ax)
    plot_2d_resid(resid, resid_range, pop_ids=resid_pop_ids,
                  extend=resid_extend)
    ax3.set_title('residuals')

    ax = pylab.subplot(2,2,4)
    flatresid = numpy.compress(numpy.logical_not(resid.mask.ravel()),
                               resid.ravel())
    ax.hist(flatresid, bins=30, density=True)
    ax.set_title('residuals')
    ax.set_yticks([])
    if show:
        pylab.show()

def divergence(params, ns, pts):
    """
    Divergence only model.
    """
    nuTX,nuFL,T,F1,F2 = params

    xx = yy = dadi.Numerics.default_grid(pts)
    n1,n2 = ns[0],ns[1]

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx,phi)

    phi = dadi.Integration.two_pops(phi,xx,T,nu1=nuTX,nu2=nuFL)
    sfs = dadi.Spectrum.from_phi_inbreeding(phi, (n1,n2), (xx,yy), (F1,F2), (2,2))
    return sfs

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

if __name__ == "__main__":
    data = dadi.Spectrum.from_file("puma.fs")
    data = data.fold()
    pts_l = [40,50,60]
    func1 = divergence
    func1_ex = dadi.Numerics.make_extrap_log_func(func1)
    func2 = divergence_noF
    func2_ex = dadi.Numerics.make_extrap_log_func(func2)

    popt = [0.45043163,0.01252053,0.01002769,0.4814076365,0.6278800995]
    popt_noF = [0.11219582,0.01005386,0.01114349]

    model = func1_ex(popt, data.sample_sizes, pts_l)
    model = model.fold()
    model_noF = func2_ex(popt_noF, data.sample_sizes, pts_l)
    model_noF = model_noF.fold()


    plot_2d_comp_multinom(model,data, fig_num=1, resid_range=475)
    #plt.savefig("puma_fit.pdf")
    #plt.close()

    plot_2d_comp_multinom(model_noF, data, fig_num=2, resid_range=475)
    #plt.savefig("puma_fit_noF.pdf")
    #plt.close()
