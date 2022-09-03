"""
"""

import os
import glob

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import h5py
import tqdm

import zcode.plot as zplot
import zcode.math as zmath
from zcode.constants import GYR, PC, MSOL, YR

import illpy

from arepa import utils, PTYPES, PNAMES, PCOLORS


'''
def get_num_snaps(sim_path):
    num_snaps = sorted(glob.glob(os.path.join(sim_path, 'snap_*.hdf5')))
    num_snaps = len(num_snaps)
    return num_snaps
'''


def draw_scatter_proj(ax, fname, ptype=1, kk=2, skip=100, param=None, bounds=None, **kwargs):
    col = PCOLORS[PTYPES.index(ptype)]
    ii, jj = utils.cross_inds(kk)

    kwargs.setdefault('s', 4)
    kwargs.setdefault('color', col)
    kwargs.setdefault('alpha', 0.05)
    kwargs.setdefault('ec', 'none')
    kwargs.setdefault('lw', 0.0)
        
    cut = slice(None) if (skip is None or ptype == illpy.PARTICLE.BH) else slice(None, None, skip)
    labels = 'xyz'
    ax.set(xlabel=labels[ii], ylabel=labels[jj])
    
    with h5py.File(fname, 'r') as h5:
        head = h5['Header']
        key = f"PartType{ptype:d}"
        try:
            part = h5[key]
        except KeyError:
            return None
            
        pos = part['Coordinates']
        pos = pos[cut]
        if bounds is not None:
            sel = True
            for ll, bb in enumerate(bounds):
                sel = sel & (bb[0] <= pos[:, ll]) & (pos[:, ll] <= bb[1])
            pos = pos[sel]
        
        aa = pos[:, ii]
        bb = pos[:, jj]            
        hh = ax.scatter(aa, bb, **kwargs)
        
    return hh


def add_time_axis(ax, times, nticks=8, format=".2f"):
    num_snaps = len(times)
        
    tw = ax.twiny()
    tick_skip = int(np.ceil(num_snaps / nticks))
    xticks = [ii for ii in range(0, num_snaps, tick_skip)]
    xticklabels = [f"{times[xt]:{format}}" for xt in xticks]
    tw.set_xlim(ax.get_xlim())
    tw.set_xticks(xticks)
    tw.set_xticklabels(xticklabels)
    return tw


class SnapImage:
    
    def __init__(self, sim_path, snap_num):
        self._sim_path = sim_path
        self._snap_num = snap_num

    def fname_snapshot(self):
        # return os.path.join(self._sim_path, f'snap_{self._snap_num:03d}.hdf5')
        return utils.fname_snapshot(self._sim_path, self._snap_num)
    
    def plot(self, bounds=None, zoom=None, track=None, **kwargs):
        nrows = 1
        if zoom is not None:
            nrows += 1
        if track:
            nrows += 1     
            
        #fig, axes = plt.subplots(figsize=[16, 5*nrows], ncols=3, nrows=nrows, squeeze=False)
        #plt.subplots_adjust(left=0.06, bottom=0.06, right=0.98, top=0.98)
        fig = plt.figure(figsize=[16, 5*nrows])  # , constrained_layout=True)
        # gs = mpl.gridspec.GridSpec(ncols=3, nrows=nrows, figure=fig, wspace=0.03, hspace=0.03)
        gs = fig.add_gridspec(ncols=3, nrows=nrows, wspace=0.15, hspace=0.15, left=0.09, bottom=0.04, right=0.88, top=0.98)

        if track is not None:
            ax_track = fig.add_subplot(gs[0, :])
            draw_bh_track(ax_track, track, highlight=self._snap_num)

        axrow_main = [fig.add_subplot(gs[-2, ii]) for ii in range(3)]
        self.plot_row(axrow_main, bounds=bounds, zoom=bounds, **kwargs)
        if zoom is not None:
            skip = kwargs.pop('skip', None)
            skip = None
            axrow_zoom = [fig.add_subplot(gs[-1, ii]) for ii in range(3)]
            self.plot_row(axrow_zoom, bounds=zoom, zoom=zoom, skip=skip, **kwargs)
            for jj in range(3):
                zplot.zoom_effect(axrow_main[jj], axrow_zoom[jj], lines=[['bl', 'tl'], ['br', 'tr']])
                axrow_main[jj].set_xlabel('')
                
            self.draw_time(axrow_zoom[1])
        else:
            self.draw_time(axrow_main[1])

        return fig
        
    def get_time(self):
        fname = self.fname_snapshot()
        with h5py.File(fname, 'r') as h5:
            head = h5['Header']
            unit_time = head.attrs['UnitLength_in_cm'] / head.attrs['UnitVelocity_in_cm_per_s']
            unit_mass = head.attrs['UnitMass_in_g']
            time = head.attrs['Time'] * unit_time / (GYR / 1000)
            
        return time
        
    def draw_time(self, ax):
        time = self.get_time()
        time = fr"$S = {self._snap_num:04d}, T = {time:7.2f} \; \mathrm{{Myr}}$"
        ax.set_title(time)
        return

    def plot_row(self, axes, bounds=None, zoom=None, **kwargs):
        fname = self.fname_snapshot()
        
        for kk, ax in enumerate(axes):
            hh = draw_scatter_proj(ax, fname, kk=kk, ptype=illpy.PARTICLE.GAS,
                                   bounds=bounds, param=None, alpha=0.1, **kwargs)

            hh = draw_scatter_proj(ax, fname, kk=kk, ptype=illpy.PARTICLE.STAR,
                                   bounds=bounds, param=None, s=10, alpha=0.25, **kwargs)

            hh = draw_scatter_proj(ax, fname, kk=kk, ptype=illpy.PARTICLE.BH,
                       bounds=bounds, param=None, s=40, marker='*', alpha=1.0, **kwargs)

            if (zoom is not None):
                ii, jj = utils.cross_inds(kk)
                ax.set(xlim=zoom[ii], ylim=zoom[jj])
                
        return
        


def load_bh_track(sim_path, skip=30, cent=[300, 300, 300]):
    snap_list = utils.get_snap_list(sim_path, skip=skip)
    num_snaps = len(snap_list)
    DIST = 0.1 * PC * 1000           # in [cm]
    fac10 = np.power(10.0, 1/3)
    
    time = np.zeros(num_snaps)
    rad = np.zeros_like(time)
    mdot = np.zeros_like(time)
    pot = np.zeros_like(time)
    num = np.zeros_like(time, dtype=int)
    num10 = np.zeros_like(time, dtype=int)   
    
    for ii, snap in enumerate(tqdm.tqdm(snap_list)):
        fname = utils.fname_snapshot(sim_path, snap)
        with h5py.File(fname, 'r') as h5:
            head = h5['Header']
            unit_dist = head.attrs['UnitLength_in_cm']
            unit_vel = head.attrs['UnitVelocity_in_cm_per_s']
            unit_time = unit_dist / unit_vel
            unit_mass = head.attrs['UnitMass_in_g']
            target_dist = DIST / unit_dist
            time[ii] = head.attrs['Time'] * unit_time / (GYR / 1000)

            part = h5['PartType5']
            pos_bh = part['Coordinates'][:]
            if np.shape(pos_bh)[0] != 1:
                raise ValueError(f"No unique BH found! {pos_bh.shape=}")

            mdot[ii] = part['BH_Mdot'][:] * (unit_mass / MSOL) / (unit_time / YR)
            pot[ii] = part['Potential'][:]
            rad[ii] = np.linalg.norm(pos_bh[0] - np.array(cent)) * unit_dist / (PC * 1000)

            pos = h5['PartType1']['Coordinates']
            sel = (
                (np.fabs(pos[:, 0] - pos_bh[0, 0]) <= target_dist*fac10) &
                (np.fabs(pos[:, 1] - pos_bh[0, 1]) <= target_dist*fac10) &
                (np.fabs(pos[:, 2] - pos_bh[0, 2]) <= target_dist*fac10)
            )

            num[ii] = np.count_nonzero(np.linalg.norm(pos[sel, :] - pos_bh, axis=1) < target_dist)
            num10[ii] = np.count_nonzero(np.linalg.norm(pos[sel, :] - pos_bh, axis=1) < target_dist*fac10)
        
    data = dict(time=time, rad=rad, mdot=mdot, pot=pot, num=num, num10=num10, snaps=snap_list)
    return data


def draw_bh_track(ax, data, highlight=None):
    time = data['time']
    rad = data['rad']
    mdot = data['mdot']
    pot = data['pot']

    kw_scatter = dict(alpha=0.5, ec='r', lw=1.0)
    if highlight is not None:
        try:
            highlight = data['snaps'].index(highlight)
        except ValueError as exception:
            warnings.warn("Snapshot for '{highlight=}' not in list {zmath.str_array(data['snaps'])}!")
            highlight = None

    hl = highlight
            
    # ax.plot(time, rad, 'k-')
    cc = 'k'
    zplot.plot_bg(ax, time, rad, color=cc)
    ax.set(xlabel='Time [Myr]', yscale='log', ylabel='Distance')
    if hl is not None:
        ax.scatter(time[hl], rad[hl], color=cc, **kw_scatter)
    
    cc = 'r'
    tw = zplot.twin_axis(ax, color=cc, pos=1.0)
    # tw.plot(time, mdot, color=cc)
    zplot.plot_bg(tw, time, mdot, color=cc)
    tw.set(yscale='log', ylabel='Mdot')
    if hl is not None:
        tw.scatter(time[hl], mdot[hl], color=cc, **kw_scatter)

    cc = 'b'
    tw = zplot.twin_axis(ax, color=cc, pos=1.05)
    # tw.plot(time, -pot, color=cc)
    zplot.plot_bg(tw, time, -pot, color=cc)
    tw.set(yscale='log', ylabel='potential')
    if hl is not None:
        tw.scatter(time[hl], -pot[hl], color=cc, **kw_scatter)

    cc = 'g'
    tw = zplot.twin_axis(ax, color=cc, pos=-0.05)
    #tw.plot(time, data['num'], color=cc, ls='--')
    #tw.plot(time, data['num10'], color=cc, ls=':')
    zplot.plot_bg(tw, time, data['num'], color=cc, ls='--')
    zplot.plot_bg(tw, time, data['num10'], color=cc, ls=':')
    tw.set(yscale='log', ylabel='neighbors')
    if hl is not None:
        tw.scatter(time[hl], data['num'][hl], color=cc, **kw_scatter)
        tw.scatter(time[hl], data['num10'][hl], color=cc, **kw_scatter)

    return 
    
