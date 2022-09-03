"""Read Arepo IO Files
"""

import os
import glob

import numpy as np

import h5py


def image_slice(fname, int_type=np.int32, double_type=np.float32, verbose=True):
    """Read in the image-slices/projections files produced directly by arepo

    These are the "<parameter>_slice_###" and "proj_<parameter>_field_###" files.

    """

    with open(fname, 'rb') as data:
        # Read header
        nums = np.fromfile(data, int_type, 2)     # read two 32 bit integers
        if verbose:
            print("num = ", nums)
        # Read density
        vals = np.fromfile(data, double_type, np.product(nums)).reshape(nums)

    return vals


def voronoi_mesh(path, snap_num=None,
                 ndim=2, int_type=np.int32, double_type=np.float32, verbose=False):

    # fname = _snap_file_from_path(fname, "voronoi_mesh_{:03d}", snap_num=snap_num)
    fname = os.path.join(path, "voronoi_mesh_{:03d}".format(snap_num))

    mesh = dict()
    with open(fname, 'rb') as data:
        # Read header
        nums = np.fromfile(data, int_type, 3)
        ngas_tot, nel_tot, ndt_tot = nums
        mesh['ngas_tot'] = ngas_tot   # total number of gas cells
        mesh['nel_tot'] = nel_tot   # total number of edges
        mesh['ndt_tot'] = ndt_tot   # total number of delaunay tetrahedra

        if verbose:
            print("num = ", nums)

        mesh['nedges'] = np.fromfile(data, int_type, ngas_tot)
        mesh['nedge_offset'] = np.fromfile(data, int_type, ngas_tot)
        mesh['edge_list'] = np.fromfile(data, int_type, nel_tot)
        mesh['xyz_edges'] = np.fromfile(data, double_type, ndt_tot*ndim)

        done = np.fromfile(data, int_type, 1)
        if len(done) != 0:
            print("done = '{}'".format(done), done.size, len(done))
            raise RuntimeError("Did not reach expected end of file!  Wrong sizes for int/float?")

    return mesh


def snapshot(path, snap_num=None, part_type=0, pars=None, verbose=False, header=False):
    if header and snap_num is None:
        snap_num = 0
    elif snap_num is None:
        raise ValueError("`snap_num` must be provided unless `header` is True!")

    fname = os.path.join(path, "snap_{:03d}.hdf5".format(snap_num))
    # fname = _snap_file_from_path(fname, "snap_{:03d}.hdf5", snap_num=snap_num)

    single_flag = isinstance(pars, str)

    part = "PartType{:1d}".format(part_type)
    with h5py.File(fname) as h5in:
        keys = list(h5in[part].keys())
        if header:
            head = {kk: vv for kk, vv in h5in['Header'].attrs.items()}
            params = {kk: vv for kk, vv in h5in['Parameters'].attrs.items()}
            return keys, head, params

        if pars is None:
            pars = keys
        if verbose:
            top_keys = h5in.keys()
            print("File keys: ", list(top_keys))
            for kk in top_keys:
                try:
                    print("\t{} keys:".format(kk), list(h5in[kk].keys()))
                    print("\t{} attrs:".format(kk), list(h5in[kk].attrs.keys()))
                except:
                    print(kk, "failed")
                    continue
            print("Particle '{}' keys: ".format(part), keys)

        data = h5in[part][pars][:] if single_flag else {kk: h5in[part][kk][:] for kk in pars}

    return data


'''
def snap_files(path):
    pattern = "snap_*.hdf5"
    fname = os.path.join(path, pattern)
    snaps = sorted(glob.glob(fname))
    return snaps


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
'''
