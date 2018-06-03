"""
"""

import numpy as np
# import matplotlib as mpl
import matplotlib.pyplot as plt

from .const import PARTS, CONV_ILL_TO_SOL

ZSOL = 0.012

part_inds = [PARTS.GAS, PARTS.DM, PARTS.STAR, PARTS.BH]
part_names = ["Gas", "DM", "Star", "BH"]
part_cuts = [20, 0, 10, 0]


part_cut_str = "".join("{}{}".format(part_names[ii][0].lower(), part_cuts[ii])
                       for ii in range(len(part_names)))

part_cut_lab = ", ".join("{}:{}".format(part_names[ii], part_cuts[ii])
                         for ii in range(len(part_names)))


def num_particles(subhalos, nbins=20):
    nparts = subhalos['SubhaloLenType']
    nparts = np.array([nparts[:, ii] for ii in part_inds]).T
    names = part_names

    # idx = (nparts_type > 0)
    # nparts_type[idx] = np.log10(nparts_type[idx])

    npar = nparts.shape[-1]

    fig, axes = plt.subplots(figsize=[16, 6], ncols=npar)
    plt.subplots_adjust(left=0.04, right=0.96, bottom=0.1, top=0.95, wspace=0.30)

    ymax = 0
    twes = []
    for ii, ax in enumerate(axes):
        log_flag = (ii < npar - 1)
        xscale = 'log' if log_flag else 'linear'
        ax.set(xscale=xscale, yscale='log', xlabel=names[ii])
        ax.tick_params(axis='y', which='both', colors='blue')

        vals = nparts[:, ii]
        vals = vals[vals > 0]
        extr = [np.min(vals), np.max(vals)]
        edges = np.logspace(*np.log10(extr), nbins) if log_flag else np.linspace(*extr, nbins)
        ax.hist(vals, bins=edges, rwidth=0.8, alpha=0.8, color='dodgerblue', edgecolor='blue')

        tw = ax.twinx()
        twes.append(tw)
        tw.tick_params(axis='y', which='both', colors='red')
        edges_fine = np.logspace(*np.log10(extr), nbins*5) if log_flag else np.linspace(*extr, nbins)
        hist, _ = np.histogram(vals, bins=edges_fine)
        csum = np.cumsum(hist)
        tw.plot(edges_fine[1:], csum, color='0.5', lw=4.0, alpha=0.6)
        tw.plot(edges_fine[1:], csum, color='red', lw=2.0, alpha=0.6)

        ymax = np.maximum(ymax, tw.get_ylim()[-1])
        if part_cuts[ii] > 0:
            ax.axvline(part_cuts[ii], lw=3.0, alpha=0.5, color='0.3', ls=(0, (6, 3)))

    for ax, tw in zip(axes, twes):
        ax.set_ylim([0.5, ymax])
        tw.set_ylim([0.5, ymax])

    return fig


def stars_galaxies(subhalos, npart_select=False):

    masses_type = subhalos['SubhaloMassType'] * CONV_ILL_TO_SOL.MASS
    mass_tot = masses_type.sum(axis=-1)
    mass_gs = masses_type[:, PARTS.GAS]
    mass_dm = masses_type[:, PARTS.DM]
    mass_st = masses_type[:, PARTS.STAR]
    nparts_type = subhalos['SubhaloLenType']
    sfr = subhalos['SubhaloSFR']
    gas_met = subhalos['SubhaloGasMetallicity']/ZSOL

    nparts = np.array([nparts_type[:, ii] for ii in part_inds]).T
    cuts = np.array(part_cuts)
    select = (nparts >= cuts[np.newaxis, :])
    select = np.all(select, axis=-1)
    if npart_select:
        print("Selecting by num particles: {}".format(part_cuts))
        mass_tot = mass_tot[select]
        mass_gs = mass_gs[select]
        mass_dm = mass_dm[select]
        mass_st = mass_st[select]
        nparts_type  = nparts_type[select, :]
        sfr = sfr[select]
        gas_met = gas_met[select]

        num_subh = len(mass_gs)
        print("Num subhalos: {}".format(num_subh))

    # Plot
    # --------------

    xlabels = ["$M_\mathrm{DM}$", "$M_\mathrm{star}$", "$M$", "$M$"]
    ylabels = ["$M_\mathrm{star}/M_\mathrm{DM}$", "$M_\mathrm{gas}/M$",
               "$\mathrm{SFR}/M_\mathrm{star}$", "$Z_\mathrm{gas}/Z_\odot$"]
    fig, axes = plt.subplots(figsize=[16, 6], ncols=4)
    plt.subplots_adjust(left=0.05, right=0.98, bottom=0.1, top=0.93)
    for ii, ax in enumerate(axes):
        xlab = xlabels[ii]
        ylab = ylabels[ii]
        ax.set(xscale='log', xlabel=xlab, yscale='log', ylabel=ylab)

    KW = {'facecolor': 'dodgerblue', 'edgecolor': 'blue', 's': 40, 'alpha': 0.6}

    if npart_select:
        label_part_cut(fig, ax)

    ax = axes[0]
    # ax.scatter(mass_dm, mass_st/mass_dm, **KW)
    scatter_with_limits(ax, mass_dm, mass_st/mass_dm, **KW)

    ax = axes[1]
    # ax.scatter(mass_st, mass_gs/mass_tot, **KW)
    scatter_with_limits(ax, mass_st, mass_gs/mass_tot, **KW)

    ax = axes[2]
    # ax.scatter(mass_tot, sfr, **KW)
    scatter_with_limits(ax, mass_tot, sfr, **KW)

    ax = axes[3]
    # ax.scatter(mass_tot, gas_met, **KW)
    scatter_with_limits(ax, mass_tot, gas_met, **KW)

    return fig


def blackhole_galaxy(subhalos, nbins=10, npart_select=False):

    mass_bh = subhalos['SubhaloBHMass'] * CONV_ILL_TO_SOL.MASS
    # masses_type = subhalos['SubhaloMassInRadType'] * CONV_ILL_TO_SOL.MASS
    masses_type = subhalos['SubhaloMassType'] * CONV_ILL_TO_SOL.MASS
    mass_st = masses_type[:, PARTS.STAR]
    nparts_type = subhalos['SubhaloLenType']
    vel_disp = subhalos['SubhaloVelDisp']

    nparts = np.array([nparts_type[:, ii] for ii in part_inds]).T
    cuts = np.array(part_cuts)
    select = (nparts >= cuts[np.newaxis, :])
    select = np.all(select, axis=-1)
    if npart_select:
        print("Selecting by num particles: {}".format(part_cuts))
        mass_st = mass_st[select]
        mass_bh = mass_bh[select]
        vel_disp = vel_disp[select]
        nparts_type  = nparts_type[select, :]

        num_subh = len(mass_bh)
        print("Num subhalos: {}".format(num_subh))

    for jj, nn in zip(part_inds, part_names):
        print("{:d} : {:4s} {:6d}".format(jj, nn, nparts_type[:, jj].sum()))

    # print("Subhalos:")
    # for ii in range(num_subh):
    #     masses = masses_type[ii]
    #     nums = nparts_type[ii]
    #     print("\n{:02d} : {:6d}, {:.2e}".format(ii, nums.sum(), masses.sum()))
    #     for jj, nn in zip(part_inds, part_names):
    #         print("\t{:2d} - {:4s} : {:6d}, {:.2e}".format(jj, nn, nums[jj], masses[jj]))
    #
    # print("")

    fig, axes = plt.subplots(figsize=[16, 6], ncols=4)
    plt.subplots_adjust(left=0.05, right=0.98, bottom=0.1, top=0.93)
    xlabels = ["MBH Mass $[M_\odot]$", "Stellar Mass $[M_\odot]$",
               "Stellar Mass $[M_\odot]$", "Vel Disp"]
    ylabels = ["Number", "Number",
               "MBH Mass $[M_\odot]$", "MBH Mass $[M_\odot]$"]
    for ii, ax in enumerate(axes):
        ax.set(xscale='log', xlabel=xlabels[ii],
               yscale='log', ylabel=ylabels[ii])

    if npart_select:
        label_part_cut(fig, ax)

    # BH Mass Hist
    ax = axes[0]
    bh_idx = (mass_bh > 0.0)
    extr = [np.min(mass_bh[bh_idx]), np.max(mass_bh[bh_idx])]
    edges = np.logspace(*np.log10(extr), nbins)
    ax.hist(mass_bh, bins=edges, color='dodgerblue', edgecolor='blue', alpha=0.6)

    # Stellar Mass Hist
    ax = axes[1]
    st_idx = (mass_st > 0.0)
    extr = [np.min(mass_st[st_idx]), np.max(mass_st[st_idx])]
    edges = np.logspace(*np.log10(extr), nbins)
    ax.hist(mass_st, bins=edges, color='dodgerblue', edgecolor='blue', alpha=0.6)

    # Stellar-Mass vs. MBH-Mass
    bin_scatter(axes[2], mass_st, mass_bh, nbins)
    scatter_with_limits(axes[2], mass_st, mass_bh)
    scatter_with_limits(axes[2], mass_st, mass_bh)

    # Velocity-Dispersion vs. MBH-Mass
    bin_scatter(axes[3], vel_disp, mass_bh, nbins)
    scatter_with_limits(axes[3], vel_disp, mass_bh)
    scatter_with_limits(axes[3], vel_disp, mass_bh)

    return fig


def scatter_with_limits(ax, xx, yy, **kwargs):
    fc_lim = '0.5'

    kwargs.setdefault('alpha', 0.5)
    kwargs.setdefault('s', 40)
    kwargs.setdefault('lw', 1.0)
    cc = kwargs.pop('color', None)
    if cc is not None:
        print("disregarding `color` parameter!")

    fc = kwargs.pop('facecolor', None)
    ec = kwargs.pop('edgecolor', None)
    if fc is None:
        fc = 'dodgerblue'

    idx = (xx > 0.0) & (yy > 0.0)
    ax.scatter(xx[idx], yy[idx], facecolor=fc, edgecolor=ec, **kwargs)

    if np.any(yy > 0):
        y0 = np.min(yy[yy > 0.0]) * 0.5
        idx = (xx > 0.0) & np.isclose(yy, 0.0)
        num = np.sum(idx)
        ax.scatter(xx[idx], np.ones(num) * y0, marker='v',
                   facecolor=fc_lim, edgecolor=ec, **kwargs)

    if np.any(xx > 0):
        x0 = np.min(xx[xx > 0.0]) * 0.5
        idx = (yy > 0.0) & np.isclose(xx, 0.0)
        num = np.sum(idx)
        ax.scatter(np.ones(num) * x0, yy[idx], marker='<',
                   facecolor=fc_lim, edgecolor=ec, **kwargs)
    return


def bin_scatter(ax, xx_in, yy_in, nbins):
    """
    """
    idx = (xx_in > 0.0) & (yy_in > 0.0)
    xx = xx_in[idx]
    yy = yy_in[idx]

    PERCS = 100 * np.array([0.159, 0.5, 0.841])

    xedges = np.logspace(*np.log10([np.min(xx), np.max(xx)]), nbins+1)
    rr = np.zeros(nbins)
    zz = np.zeros((nbins, PERCS.size))
    ave = np.zeros(nbins)
    for ii in range(nbins):
        x0 = xedges[ii]
        x1 = xedges[ii+1]
        rr[ii] = np.power(10, 0.5*(np.log10(x0) + np.log10(x1)))
        if ii == 0:
            idx = ((x0 <= xx) | np.isclose(x0, xx)) & (xx < x1)
        else:
            idx = (x0 < xx) & ((xx <= x1) | np.isclose(xx, x1))

        if idx.sum() == 0:
            continue

        ave[ii] = np.average(yy[idx])
        zz[ii, :] = np.percentile(yy[idx], PERCS)

    idx = (ave > 0)
    ax.plot(rr[idx], ave[idx], 'r-', lw=1.0, alpha=0.8)
    ax.plot(rr[idx], zz[idx, 1], 'r--', lw=1.0, alpha=0.8)
    ax.fill_between(rr[idx], zz[idx, 0], zz[idx, -1], color='red', alpha=0.5)
    return


def label_part_cut(fig, ax):
    ax.text(0.5, 0.99, part_cut_lab,
            horizontalalignment='center', verticalalignment='top',
            transform=fig.transFigure)
    return
