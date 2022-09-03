"""
"""

import os
import glob
import sys
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt

import zcode.inout as zio
import zcode.plot as zplot

import illpy

_CDIR = os.path.split(os.path.abspath(__file__))[0]

# DNAME_ANALYSIS_OUTPUT = "analysis"
PTYPES = [illpy.PARTICLE.DM, illpy.PARTICLE.GAS, illpy.PARTICLE.STAR, illpy.PARTICLE.BH]
PNAMES = ['dm', 'gs', 'st', 'bh']
# pcolors = ['0.5', 'r', 'blue', 'k']
PCOLORS = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3"]

from . import readio  # noqa
from . import plot  # noqa
from . import isolated # noqa
from .isolated import *





'''
import socket
_hostname = socket.gethostname().lower()
ILLUSTRIS_PYTHON_PATH = None
if 'daedalus' in _hostname:
    HOST = 'daedalus'
elif 'harvard' in _hostname:
    HOST = 'odyssey'
    ILLUSTRIS_PYTHON_PATHS = "/n/home00/lkelley/illustris/"
elif 'ufhpc' in _hostname:
    HOST = 'hipergator'
    ILLUSTRIS_PYTHON_PATH = "/home/lkelley/arepo_sims"

    import matplotlib
    matplotlib.use('Agg')

else:
    print("WARNING: unrecognized host!")

# Add Path to `illustris_python` module
# ILLUSTRIS_PYTHON_PATH = "/n/home00/lkelley/illustris/"
# if ILLUSTRIS_PYTHON_PATH not in sys.path:
#     sys.path.append(ILLUSTRIS_PYTHON_PATH)
#     print("Added '{}' `to sys.path`".format(ILLUSTRIS_PYTHON_PATH))

if (ILLUSTRIS_PYTHON_PATH is not None) and (ILLUSTRIS_PYTHON_PATH not in sys.path):
    if os.path.exists(ILLUSTRIS_PYTHON_PATH):
        sys.path.append(ILLUSTRIS_PYTHON_PATH)
        print("Added '{}' `to sys.path`".format(ILLUSTRIS_PYTHON_PATH))
    else:
        print("WARNING: illustris path '{}' does not exist!".format(ILLUSTRIS_PYTHON_PATH))

# Make sure required paths exist
for path in _check_paths:
    if not os.path.exists(path):
        os.makedirs(path)
        print("Created '{}'".format(path))
'''


class Arepo_Sim_Path:

    _DNAME_OUTPUT = "output"
    _DNAME_ANALYSIS = "analysis"

    def __init__(self, sim_path, create=True):
        # Make sure the path does NOT include a terminal path-separator
        sim_path = os.path.dirname(os.path.join(sim_path, ''))
        # If the 'output' folder is specified, move up a directory to the base-path
        if sim_path.endswith(self._DNAME_OUTPUT):
            sim_path = os.path.dirname(sim_path)

        self.sim = sim_path
        self.output = os.path.join(sim_path, self._DNAME_OUTPUT)
        # self.analysis = os.path.join(sim_path, self._DNAME_ANALYSIS)
        self.analysis = os.path.join(self.output, self._DNAME_ANALYSIS)

        _check_paths = [self.output, self.analysis]
        if create:
            for cp in _check_paths:
                zio.check_path(cp)

        return

    def __str__(self):
        return self.sim

    def set_snapshot_times(self, fname):
        scales = []
        with open(fname, 'r') as times:
            for line in times.readlines():
                line = line.strip()
                comps = line.split()
                scales.append(float(comps[0]))

        self.scales = np.array(scales)
        self.redz = 1.0/self.scales - 1.0
        return

    def snapshots(self, snap=None):
        if snap is None:
            fnames = snap_files(self.output)
        else:
            if np.isscalar(snap):
                squeeze = True
            snap = np.atleast_1d(snap)
            fnames = [os.path.join(self.output, "snap_{:03d}.hdf5".format(sn))
                      for sn in snap]
            if squeeze:
                fnames = fnames[0]

        return fnames

    def num_snaps(self):
        return len(self.snapshots())

    def fname(self, fname, path=None, subdir=None, snap_num=None, modify=False):
        if path is None:
            path = self.analysis

        if subdir is True:
            subdir = fname.split('.')[0]

        if subdir is not None:
            path = os.path.join(path, subdir)

        # Try to 'format' the snap_num in place (e.g. "voronoi_mesh_{:03d}" ==> "voronoi_mesh_023")
        #   otherwise, append snap (e.g. "voronoi_mesh" ==> "voronoi_mesh_snap-023")
        if snap_num is not None:
            import string
            if string.Formatter().parse(fname)[1] is not None:
                fname = fname.format(snap_num)
            else:
                fname = fname_append_snap(fname, snap_num)

        fname = os.path.join(path, fname)

        if modify:
            fname = zio.modify_exists(fname)

        return fname

    def save_fig(self, fig, fname, path=None, subdir=None, snap_num=None, modify=False,
                 verbose=False, **kwargs):
        fname = self.fname(
            fname, path=path, subdir=subdir, snap_num=snap_num, modify=modify)
        kwargs.setdefault('transparent', True)
        kwargs.setdefault('dpi', 100)
        zio.check_path(fname)
        fig.savefig(fname, **kwargs)
        if verbose:
            print("saved to '{}' size: {}".format(fname, zio.get_file_size(fname)))
        return fname

    '''
    def snap_file(self, fname, snap_num=None):
        if snap_num is not None:
            fname = fname.format(snap_num)

        if not os.path.exists(fname):
            raise ValueError("Path '{}' does not exist!".format(fname))

        if os.path.isdir(fname):
            if snap_num is None:
                raise ValueError("If path is provided, `snap_num` must be given!")
            fname = os.path.join(fname, pattern.format(snap_num))

        if not os.path.exists(fname):
            raise ValueError("File '{}' does not exist!".format(fname))

        return fname
    '''


def _snap_file_from_path(fname, pattern, snap_num=None):
    if not os.path.exists(fname):
        raise ValueError("Path '{}' does not exist!".format(fname))

    if os.path.isdir(fname):
        if snap_num is None:
            raise ValueError("If path is provided, `snap_num` must be given!")
        fname = os.path.join(fname, pattern.format(snap_num))

    if not os.path.exists(fname):
        raise ValueError("File '{}' does not exist!".format(fname))

    return fname


def snap_files(path):
    # patterns = ["snap_[0-9]{3}.hdf5", "snap_[0-9]{3}", "snapdir_[0-9]{3}"]
    patterns = ["snap_*.hdf5", "snap_" + "[0-9]"*3, "snapdir_" + "[0-9]"*3]
    for pat in patterns:
        fname = os.path.join(path, pat)
        snaps = sorted(glob.glob(fname))
        if len(snaps) > 0:
            break

        # print("No files found matching '{}'".format(fname))

    return snaps


'''
def path_analysis_from_output(path_output):
    if os.path.basename(path_output).lower() != 'output':
        raise ValueError("Path doesn't look right!  '{}'".format(path_output))

    path = os.path.join(path_output, os.path.pardir, DNAME_ANALYSIS_OUTPUT, '')
    return path
'''


def get_num_snaps(arepo_output_dir, pattern="snap*"):
    # Get number of snapshots
    pattern = os.path.join(arepo_output_dir, pattern)
    snap_fnames = sorted(glob.glob(pattern))
    if len(snap_fnames) < 1:
        print("Pattern: '{}'".format(pattern))
        print("Directory contents: ", os.listdir(arepo_output_dir))
        raise RuntimeError("No snapshot directories found!")

    num_snaps = len(snap_fnames)
    return num_snaps


def fname_append_snap(fname, snap):
    if isinstance(snap, str):
        app = "_snap-{}".format(snap)
    else:
        app = "_snap-{:03d}".format(snap)

    fname = zio.modify_filename(fname, append=app)
    return fname


'''
def save_fig(fig, fname, path=None, subdir=None, snap_num=None, modify=False,
             transparent=True, dpi=100, **kwargs):

    # if os.path.basename(path).lower() == 'output':
    #    path = os.path.join(path, os.path.pardir, DIR_ANALYSIS_OUTPUT, '')
    try:
        path = path_analysis_from_output(path)
    except:
        pass

    if (snap_num is not None) and (subdir is None):
        comps = fname.split('.')
        if len(comps) != 2:
            raise ValueError("Bad `fname` = '{}'  ({})!".format(fname, comps))
        subdir = comps[0]
        fname = subdir + "_snap-{:03d}".format(snap_num) + "." + comps[1]

    fname = zplot.save_fig(fig, fname, path=path, modify=modify,
                           subdir=subdir, transparent=transparent, dpi=dpi, **kwargs)
    return fname
'''


def main(arepo_sim_dir, snap=-1):
    arepo_output_dir = os.path.join(arepo_sim_dir, "output", "")
    sim_name = os.path.split(arepo_sim_dir)[-1]
    print("Sim name: '{}'".format(sim_name))
    print("Running on directory: '{}'".format(arepo_output_dir))
    num_snaps = get_num_snaps(arepo_output_dir)
    print("\t{} Snapshots".format(num_snaps))

    # print(os.listdir(arepo_output_dir))
    # image_slices(arepo_output_dir)
    if snap < 0:
        snap = num_snaps + snap

    subfind_at_snapshot(arepo_output_dir, snap, sim_name)

    return


def subfind_at_snapshot(arepo_output_dir, snapshot, sim_name):
    beg = datetime.now()
    subhalos = ill.groupcat.loadSubhalos(arepo_output_dir, snapshot)
    end = datetime.now()
    dur = (end - beg)
    print("Loaded {} subhalos after {}".format(subhalos['count'], str(dur)))
    print("\tkeys: '{}'\n".format(subhalos.keys()))

    fig = plot.stars_galaxies(subhalos, npart_select=False)
    save_fig(fig, 'subhalos_galaxies.pdf', sim_name)

    fig = plot.stars_galaxies(subhalos, npart_select=True)
    save_fig(fig, 'subhalos_galaxies_{}.pdf'.format(plot.part_cut_str), sim_name)

    fig = plot.num_particles(subhalos)
    save_fig(fig, 'subhalos_num-parts_{}.pdf'.format(plot.part_cut_str), sim_name)

    fig = plot.blackhole_galaxy(subhalos, npart_select=False)
    save_fig(fig, 'subhalos_blackholes-hosts.pdf', sim_name)

    fig = plot.blackhole_galaxy(subhalos, npart_select=True)
    save_fig(fig, 'subhalos_blackholes-hosts_{}.pdf'.format(plot.part_cut_str), sim_name)
    return


def image_slices(arepo_output_dir, dest=".", fname="density_field_",
                 log_flag=False, interp='nearest'):
    """Draw the image slices produced directly by Arepo (e.g. 'proj_density_field_000').

    These slices are enabled with the `VORONOI_IMAGES_FOREACHSNAPSHOT` option.

    """
    # fname = "proj_density_field_"
    fname = "density_slice_"    # for some reason, some of these values are negative?
    # fname = "metal_slice_"
    # fname = "chemelements_slice_"
    pattern = os.path.join(arepo_output_dir, fname + "*")
    files = sorted(glob.glob(pattern))
    print("Found '{}' files matching pattern '{}'".format(len(files), pattern))

    DOUBLE_TYPE = np.float32
    INT_TYPE = np.int32

    path_output = os.path.abspath(dest)
    print("Outputting to '{}'".format(path_output))

    # for fname in [files[-1]]:
    for fname in files:
        fnum = os.path.split(fname)[-1].split("_")[-1]

        with open(fname, 'rb') as data:
            # Read header
            nums = np.fromfile(data, INT_TYPE, 2)     # read two 32 bit integers
            print("num = ", nums)
            # Read density
            dens = np.fromfile(data, DOUBLE_TYPE, np.product(nums)).reshape(nums)
            print("density extr = ", np.min(dens), np.max(dens))

        fig, ax = plt.subplots(figsize=[10, 10])
        vals = np.log10(dens) if log_flag else dens
        ax.imshow(vals, interpolation=interp)
        # fname_out = '{}{}.pdf'.format(fname, fnum)
        fname_out = '{}.pdf'.format(os.path.basename(fname))
        # fname_out = os.path.join(PATH_OUTPUT, fname_out)
        fname_out = os.path.join(path_output, fname_out)
        fig.savefig(fname_out)
        print("Saved to '{}'".format(fname_out))
        plt.close('all')

    return


if __name__ == "__main__":
    print(sys.argv)
    if len(sys.argv) > 1:
        arepo_sim_dir = sys.argv[1]
    else:
        raise RuntimeError("No directory provided!")

    main(arepo_sim_dir)
