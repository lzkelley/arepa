"""
"""

import os
import glob
import subprocess
import logging

import numpy as np
import matplotlib as mpl
import matplotlib.patches   # noqa
import matplotlib.pyplot as plt
import scipy as sp
import scipy.stats  # noqa

from numba import njit
import tqdm.notebook as tqdm

import zcode.plot as zplot
import zcode.math as zmath

from .const import PARTS, CONV_ILL_TO_SOL

# import arepo_analysis as lysis
from . import readio
from . import plot

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


# ===========================  NEW STUFF ========================================


def rad_proj_from_2d(path, snap_num, scatter_param='Density', hist_param='Masses',
                     hist_dens=True, extr=None, nbins=20):
    bin_sigma = 3

    _params = ['Coordinates', 'Masses']
    params = [pp for pp in _params]
    if (scatter_param is not None) and (scatter_param not in _params):
        params.append(scatter_param)
    if (hist_param is not None) and (hist_param not in _params):
        params.append(hist_param)

    snap = readio.snapshot(path, snap_num=snap_num, pars=params)
    pos = snap['Coordinates']
    mass = snap['Masses']

    scat = None
    if scatter_param is not None:
        scat = snap[scatter_param]
    hist = None
    if hist_param is not None:
        hist = snap[hist_param]

    com = np.sum(pos*mass[:, np.newaxis], axis=0) / np.sum(mass)
    rads = np.linalg.norm(pos - com, axis=-1)
    if extr is None:
        extr = zmath.quantiles(rads, sigmas=[-bin_sigma, bin_sigma])

    edges = None
    if hist is not None:
        stat = 'sum' if hist_dens else 'mean'

        edges = zmath.spacing(extr, 'log', nbins+1)
        hist, *_ = sp.stats.binned_statistic(rads, hist, statistic=stat, bins=edges)
        if hist_dens:
            da = np.diff(np.pi*edges**2)
            hist = hist / da

    return rads, scat, hist, edges


def plot_proj_1d(path, scatter_param='Density', hist_param='Masses', snaps=None, hist_dens=True,
                 extr=None, nbins=20, xlim=None, ylim=None, movie=True, framerate=None):
    num_snaps = len(readio.snap_files(path))
    fname = "proj1d_{}-{}.png".format(scatter_param.lower(), hist_param.lower())

    init = None
    mean = None
    if snaps is None:
        snaps = range(num_snaps)

    yextr = None
    for snap in tqdm.tqdm(snaps):
        rads, scat, hist, edges = plot.rad_proj_from_2d(
            path, snap, scatter_param=scatter_param, hist_param=hist_param,
            hist_dens=hist_dens, extr=extr, nbins=nbins)
        if extr is None:
            extr = zmath.minmax(edges)
        if init is None:
            init = np.copy(hist)

        idx = (hist > 0.0)
        if mean is None:
            mean = np.zeros_like(hist)
            mean[idx] = hist[idx]
        else:
            mean[idx] = mean[idx] + hist[idx]

        yextr = zmath.minmax(scat, prev=yextr, log_stretch=0.1, filter='>')

    mean = mean / num_snaps
    if ylim is None:
        ylim = yextr

    output_fnames = []
    for snap_num in tqdm.tqdm(snaps):
        fig = plot.plot_proj_1d_snap(
            path, snap_num, scatter_param=scatter_param, hist_param=hist_param,
            extr=extr, hist_dens=hist_dens, mean=mean, init=init)

        fig.axes[0].set(xlim=xlim, ylim=ylim)
        _fname = save_fig(fig, fname, path, snap_num=snap_num, verbose=False)
        if _fname is None:
            print("BAD!")
            print(_fname, fname)
            raise RuntimeError()

        output_fnames.append(_fname)
        plt.close('all')

    if framerate is None:
        framerate = len(snaps) / 20
        framerate = int(np.clip(framerate, 2, 15))

    print("Saved to (e.g.) '{}'".format(output_fnames[0]))
    if movie:
        movie_fname = make_movie(path, output_fnames, fname.replace('.png', '.mp4'),
                                 framerate=framerate)
        print("Saved movie to '{}'".format(movie_fname))

    return


def plot_proj_1d_snap(path, snap_num, scatter_param='Density', hist_param='Masses',
                      mean=None, init=None, hist_dens=True, extr=None, ave=None, nbins=20):

    fig, ax = zplot.figax(left=0.05, right=0.98, bottom=0.1, top=0.98)
    keys, head, pars = readio.snapshot(path, snap_num=snap_num, header=True)
    plot.draw_snap_time(ax, head, num=snap_num)

    rads, scat, hist, edges = rad_proj_from_2d(
        path, snap_num, scatter_param=scatter_param, hist_param=hist_param,
        hist_dens=hist_dens, extr=extr, nbins=nbins)

    if hist is not None:
        zplot.plot_hist_line(ax, edges, hist, color='r', zorder=10, alpha=0.6, lw=2.0)

    if mean is not None:
        zplot.plot_hist_line(ax, edges, mean, color='r', ls='--', zorder=10, alpha=0.6, lw=2.0)

    if init is not None:
        zplot.plot_hist_line(ax, edges, init, color='k', ls='--', zorder=10, alpha=0.3, lw=2.0)

    if scat is not None:
        ax.scatter(rads, scat, alpha=0.4, s=10, facecolor='b')

    return fig


def plot_mesh_2d(path, snap_num=0, param='Density', parse_param=None,
                 region=None, periodic=None, smap=None):

    try:
        path = path.output
    except:
        path = str(path)

    mesh = readio.voronoi_mesh(path, snap_num=snap_num)
    # snap = readio.snapshot(path_input, snap_num=snap_num, pars=['Coordinates', param])
    # pos = snap['Coordinates']
    snap = readio.snapshot(path, snap_num=snap_num, pars=[param])
    # pos = snap['Coordinates']
    vals = snap[param]
    if parse_param is not None:
        vals = parse_param(vals)

    fig, ax = plt.subplots(figsize=[12, 10])
    plt.subplots_adjust(left=0.05, right=0.95, bottom=0.0, top=0.98)
    smap, extr = draw_mesh(
        ax, mesh, vals=vals, fix_poly=True, lines_flag=True,
        lw=0.25, color='k', alpha=0.25, region=region, periodic=periodic, smap=smap)

    # ax.scatter(*pos[:, :2].T, color='r', s=2, alpha=0.25)
    if region is None and periodic is not None:
        region = periodic
    if region is not None:
        ax.set(xlim=region[0], ylim=region[1])

    plt.colorbar(smap, ax=ax, orientation='horizontal', pad=0.05)
    return fig


def draw_mesh(ax, mesh, vals=None, fix_poly=False,
              smap=None, region=None, periodic=None, lines_flag=True, **kwargs):

    if (not lines_flag) and (vals is None):
        raise ValueError("Nothing is being plotted!")

    ndt = mesh['ndt_tot']
    xyz = mesh['xyz_edges']
    edge_list = mesh['edge_list']
    nedge_offset = mesh['nedge_offset']
    nedges = mesh['nedges']
    NDIM = 2
    NSIDE_MAX = 20

    xyz = xyz.reshape(ndt, NDIM)

    poly = np.zeros((NSIDE_MAX, NDIM))
    if lines_flag:
        lines = np.full((2*len(edge_list), NDIM), np.nan)
    tot_num = len(nedge_offset)

    mult = 1 if (periodic is None) else 2

    if periodic is not None:
        periodic = [np.array(pp) if pp is not None else pp for pp in periodic]
        # np.atleast_2d(periodic)
    if region is not None:
        # region = [np.array(pp) if pp is not None else pp for pp in region]
        region = np.atleast_2d(region)

    if vals is not None:
        patches = np.empty(mult*tot_num, dtype=object)
        colors = np.zeros(mult*tot_num)

    cnt = 0
    valid = np.zeros(mult*tot_num, dtype=bool)

    def add_cell(ee, ne, poly, cnt, end=False):
        if end:
            ff = mult*tot_num - 1 - ee
        else:
            ff = ee

        if lines_flag:
            lines[cnt:cnt+ne, :] = poly[:ne, :]

        if vals is not None:
            inc = 0
            if fix_poly and np.allclose(poly[0, :], poly[ne-1, :]):
                ne = ne - 1
                inc = 1

            pat = mpl.patches.Polygon(poly[:ne])
            patches[ff] = pat
            colors[ff] = vals[ee]
            ne = ne + inc

        valid[ff] = True
        cnt = cnt + ne + 1
        return cnt

    pers = 0
    for ee in tqdm.tqdm(range(tot_num), total=tot_num, leave=False):
        oo = nedge_offset[ee]
        ne = nedges[ee]
        ll = edge_list[oo]
        lo = xyz[ll]
        poly[0] = lo
        if ne >= NSIDE_MAX:
            err = "Number of edges for element {} = {}, exceeds max {}!".format(ee, ne, NSIDE_MAX)
            raise ValueError(err)

        for ff in range(1, ne):
            hh = edge_list[oo+ff]
            hi = xyz[hh]
            poly[ff] = hi

        if (region is not None) and (not any_within(poly[:ne], region)):
            continue

        cnt = add_cell(ee, ne, poly, cnt, end=False)

        if periodic is None:
            continue

        for dd in range(NDIM):
            if periodic[dd] is None:
                continue

            if np.any((poly[:ne, dd] < periodic[dd][0])):
                dup = np.copy(poly[:ne, :])
                dup[:, dd] += (periodic[dd][1] - periodic[dd][0])
                cnt = add_cell(ee, ne, dup, cnt, end=True)
                pers += 1
            elif np.any(poly[:ne, dd] > periodic[dd][1]):
                dup = np.copy(poly[:ne, :])
                dup[:, dd] -= (periodic[dd][1] - periodic[dd][0])
                cnt = add_cell(ee, ne, dup, cnt, end=True)
                pers += 1

        # if cnt > 1000:
        #     break

    extr = zmath.minmax(colors[valid])
    if vals is not None:
        if smap is None:
            smap = zplot.smap(extr, cmap='viridis')

        p = mpl.collections.PatchCollection(patches[valid], cmap=smap.cmap, norm=smap.norm)
        p.set_array(colors[valid])
        ax.add_collection(p)

    if lines_flag:
        lines = lines[:cnt, :].T
        ax.plot(*lines, **kwargs)

    return smap, extr


def make_movie(path, fname, output, framerate=10):
    '''
    try:
        path_out = lysis.path_analysis_from_output(path)
    except ValueError:
        path_out = path
    '''
    try:
        path = path.output
    except:
        path = str(path)

    # If `fname` is not a list of complete-filenames then try to glob them by appending 'snap-*'
    if not isinstance(fname, list):
        fname = fname_append_snap(fname, '*')
        file_list = sorted(glob.glob(fname))
        if len(file_list) == 0:
            fname = fname_append_snap(os.path.join(path, fname), '*')
            file_list = sorted(glob.glob(fname))
    else:
        file_list = fname

    if len(file_list) < 2:
        raise RuntimeError("Forming file_list failed!")

    file_list_fname = "_temp_file_list.txt"
    with open(file_list_fname, 'w') as out:
        for fn in file_list:
            out.write("file " + fn + "\n")

    otpath = os.path.join(path, output)
    cmd = [
        'ffmpeg',
        '-f', 'concat',
        '-safe', '0',
        '-r', str(framerate),
        '-i', file_list_fname,  # inpath,
        '-r', str(framerate),
        '-pix_fmt', 'yuv420p', otpath, '-y'
    ]
    try:
        print(" ".join(cmd))
        output = subprocess.check_output(cmd)  # , stderr=logging.error)
    except Exception as err:
        logging.error("`plot.py:movie()` failed!  command: '{}'".format(" ".join(cmd)))
        logging.error(str(output))
        logging.error(str(err))
        return None
    finally:
        os.remove(file_list_fname)
        pass

    return otpath


def draw_snap_time(ax, head, num=None, format='5.1f', **kwargs):
    time = float(head['Time'])
    label = "{{:{:}}}".format(format)
    label = label.format(time)
    if num is not None:
        label = "{:03d}: {}".format(num, label)

    kwargs.setdefault('loc', 'ur')
    kwargs.setdefault('fontsize', 12)
    kwargs.setdefault('alpha', 0.5)
    return zplot.text(ax, label, **kwargs)


@njit()
def any_within(poly, bounds):
    pnts, dims = np.shape(poly)
    for ii in range(pnts):
        test = 0
        for jj in range(dims):
            if ((poly[ii, jj] > bounds[jj, 0]) and (poly[ii, jj] < bounds[jj, 1])):
                test += 1
        if test == dims:
            return True

    return False


@njit()
def any_without(poly, bounds):
    pnts, dims = np.shape(poly)
    for ii in range(pnts):
        for jj in range(dims):
            if (poly[ii, jj] < bounds[jj, 0]):
                return True
            if (bounds[jj, 1] < poly[ii, jj]):
                return True

    return False
