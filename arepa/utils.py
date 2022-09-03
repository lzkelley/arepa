"""
"""

import os
import glob

import numpy as np
import h5py

# import illustris_python as ill
import illpy


def get_num_snaps(arepo_output_dir, pattern="snapdir_*"):
    # Get number of snapshots
    pattern = os.path.join(arepo_output_dir, pattern)
    snap_fnames = sorted(glob.glob(pattern))
    if len(snap_fnames) < 1:
        print("Pattern: '{}'".format(pattern))
        print("Directory contents: ", os.listdir(arepo_output_dir))
        raise RuntimeError("No snapshot directories found!")

    num_snaps = len(snap_fnames)
    return num_snaps


def fname_snapshot(sim_path, snap_num, format='snap_{:03d}.hdf5'):
    other_formats = ['snap_{:03d}.hdf5', 'snapshot_{:03d}.hdf5', 'snap_{:04d}.hdf5', 'snapshot_{:04d}.hdf5']
    all_formats = [format] + other_formats

    sim_path = sim_path.rstrip('/').lower()
    if not sim_path.endswith('output'):
        sim_path = os.path.join(sim_path, 'output', '')
    if not os.path.isdir(sim_path):
        raise FileNotFoundError("{sim_path=} does not exist!")
        
    for form in all_formats:
        fname = os.path.join(sim_path, form.format(snap_num))
        if os.path.isfile(fname):
            break
    else:
        raise FileNotFoundError("Could not find a snapshot '{snap_num}' in '{sim_path}'; tried '{all_formats}'!")

    return fname


def get_snap_list(sim_path, skip=1, pattern='snap_*'):
    num_snaps = get_num_snaps(sim_path, pattern=pattern)
    snap_list = list(range(0, num_snaps, skip))
    if num_snaps-1 not in snap_list:
        snap_list.append(num_snaps-1)
    return snap_list


def cross_inds(kk, reverse=True):
    if reverse:
        kk = 3-1-kk
    ii = (kk + 1) % 3
    jj = (kk + 2) % 3
    ii, jj = sorted([ii, jj])
    return ii, jj


def load_snapshot_redshifts(arepo_output_dir):
    fname = os.path.join(arepo_output_dir, os.path.pardir, "output_list.txt")

    num_snaps = get_num_snaps(arepo_output_dir)
    redshifts = []
    for snap in range(num_snaps):
        with h5py.File(illpy.groupcat.gcPath(arepo_output_dir, snap), 'r') as ff:
            header = dict(ff['Header'].attrs.items())
            redz = np.float(header['Redshift'])
            redshifts.append(redz)            

    return redshifts
