"""
"""

import os
import glob

import numpy as np
import h5py

import illustris_python as ill


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


def load_snapshot_redshifts(arepo_output_dir):
    fname = os.path.join(arepo_output_dir, os.path.pardir, "output_list.txt")

    num_snaps = get_num_snaps(arepo_output_dir)
    redshifts = []
    for snap in range(num_snaps):
        with h5py.File(ill.groupcat.gcPath(arepo_output_dir, snap), 'r') as ff:
            header = dict(ff['Header'].attrs.items())
            redz = np.float(header['Redshift'])
            redshifts.append(redz)            

    return redshifts
