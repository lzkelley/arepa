"""
"""
import os
import glob
import sys
from datetime import datetime

import numpy as np
from matplotlib import pyplot as plt

import illustris_python as ill

from . import plot
from . import PATH_OUTPUT


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


def image_slices(arepo_output_dir, log_flag=False, interp='nearest'):
    """Draw the image slices produced directly by Arepo (e.g. 'proj_density_field_000').

    These slices are enabled with the `VORONOI_IMAGES_FOREACHSNAPSHOT` option.

    """
    # image_name = "proj_density_field_"
    # image_name = "density_slice_"    # for some reason, some of these values are negative?
    # image_name = "metal_slice_"
    image_name = "chemelements_slice_"
    files = sorted(glob.glob(os.path.join(arepo_output_dir, image_name + "*")))

    DOUBLE_TYPE = np.float32
    INT_TYPE = np.int32

    # for fname in [files[-1]]:
    for fname in files:
        fnum = os.path.split(fname)[-1].split("_")[-1]

        with open(fname, 'rb') as data:
            # Read header
            nums = np.fromfile(data, INT_TYPE, 2)     # read three 32 bit integers
            # print("num = ", nums)
            # Read density
            dens = np.fromfile(data, DOUBLE_TYPE, np.product(nums)).reshape(nums)
            print("density extr = ", np.min(dens), np.max(dens))

        fig, ax = plt.subplots(figsize=[10, 10])
        vals = np.log10(dens) if log_flag else dens
        ax.imshow(vals, interpolation=interp)
        fname_out = '{}{}.pdf'.format(image_name, fnum)
        fname_out = os.path.join(PATH_OUTPUT, fname_out)
        fig.savefig(fname_out)
        print("Saved to '{}'".format(fname_out))
        plt.close('all')

    return


def get_num_snaps(arepo_output_dir):
    # Get number of snapshots
    pattern = os.path.join(arepo_output_dir, "snapdir_*")
    snap_fnames = sorted(glob.glob(pattern))
    if len(snap_fnames) < 1:
        print("Pattern: '{}'".format(pattern))
        print("Directory contents: ", os.listdir(arepo_output_dir))
        raise RuntimeError("No snapshot directories found!")

    num_snaps = len(snap_fnames)
    return num_snaps


def save_fig(fig, fname, sim_name):
    fname = os.path.join(PATH_OUTPUT, sim_name, fname)
    path = os.path.split(fname)[0]
    if not os.path.exists(path):
        os.makedirs(path)

    fig.savefig(fname)
    print("wrote to '{}'".format(fname))
    return fname


if __name__ == "__main__":
    print(sys.argv)
    if len(sys.argv) > 1:
        arepo_sim_dir = sys.argv[1]
    else:
        raise RuntimeError("No directory provided!")

    main(arepo_sim_dir)
