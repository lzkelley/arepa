"""
"""

import os
import sys
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt
import h5py

import illustris_python as ill

import zcode.math as zmath

from . import utils
from . const import CONV_ILL_TO_SOL, PARTS


def subhalos_metallicity_redshift(arepo_sim_dir):
    beg_all = datetime.now()
    arepo_output_dir = os.path.join(arepo_sim_dir, "output", "")
    sim_name = os.path.split(arepo_sim_dir)[-1]
    print("subhalo_metals.main() : {} - {}".format(sim_name, arepo_output_dir))
    
    HALOS = True

    redshifts = utils.load_snapshot_redshifts(arepo_output_dir)
    # num_snaps = utils.get_num_snaps(arepo_output_dir)
    num_snaps = len(redshifts)
    print("Snaps: ", num_snaps)

    fields = ['SubhaloBHMass', 'SubhaloGasMetallicity', 'SubhaloGasMetallicityHalfRad',
              'SubhaloGasMetallicityMaxRad', 'SubhaloGasMetallicitySfr',
              'SubhaloGasMetallicitySfrWeighted', 'SubhaloLenType', 'SubhaloMass',
              'SubhaloMassInRadType', 'SubhaloSFR', 'SubhaloSFRinRad', 'SubhaloStarMetallicity',
              'SubhaloVelDisp']

    NBINS_MASSES = 60
    NBINS_METALS = 80
    edges_metals = np.logspace(-8, 0, NBINS_METALS+1)
    edges_masses = np.logspace(8, 14, NBINS_MASSES+1)
    
    # SIGMA = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
    SIGMA = [0.0, 1.0, 2.0, 3.0]
    PERCS = np.array(sorted(list(set(np.concatenate(zmath.sigma(SIGMA, boundaries=True))))))
    NUM_PERCS = len(PERCS)

    METALS_KEYS = ['SubhaloGasMetallicity', 'SubhaloGasMetallicityHalfRad',
                   'SubhaloGasMetallicitySfr', 'SubhaloGasMetallicitySfrWeighted']

    NUM_KEYS = len(METALS_KEYS)
    metals = np.zeros((NUM_KEYS, num_snaps, NBINS_METALS, NBINS_MASSES))

    PERCS_MASS_BIN_EDGES = np.logspace(8, 14, 7)
    NUM_PMB = PERCS_MASS_BIN_EDGES.size - 1
    metals_percs = np.zeros((NUM_KEYS, num_snaps, NUM_PMB, NUM_PERCS))
    metals_percs_num = np.zeros((NUM_KEYS, num_snaps, NUM_PMB), dtype=int)
    metals_percs_nonzero = np.zeros((NUM_KEYS, num_snaps, NUM_PMB, NUM_PERCS))
    metals_percs_nonzero_num = np.zeros((NUM_KEYS, num_snaps, NUM_PMB), dtype=int)
    # masses_percs = np.zeros((NUM_KEYS, num_snaps, NBINS_MASSES, NUM_PERCS))

    edges = [edges_metals, edges_masses]

    for snap in range(num_snaps):
        beg = datetime.now()
        print(snap, "  ------")

        subhalos = ill.groupcat.loadSubhalos(arepo_output_dir, snap, fields=fields)
        # print(subhalos.keys())

        masses_snap = subhalos['SubhaloMass'] * CONV_ILL_TO_SOL.MASS
        num_subh = len(masses_snap)
        print("\tLoaded {} subhalos after {}".format(num_subh, str(datetime.now()-beg)))

        for kk, key in enumerate(METALS_KEYS):
            # metals_snap = subhalos['SubhaloGasMetallicity']
            metals_snap = subhalos[key]
            print("\t", key)
            print("\t\tminmax = {:.2e}, {:.2e}".format(np.min(metals_snap), np.max(metals_snap)))
            metals[kk, snap, :, :], xe, ye = np.histogram2d(metals_snap, masses_snap, bins=edges)

            for jj in range(NUM_PMB):
                lo = PERCS_MASS_BIN_EDGES[jj]
                hi = PERCS_MASS_BIN_EDGES[jj+1]

                idx = (lo < masses_snap) & (masses_snap <= hi)
                if any(idx):
                    metals_percs[kk, snap, jj, :] = np.percentile(metals_snap[idx], 100*PERCS)
                    metals_percs_num[kk, snap, jj] = np.sum(idx)

                idx = (lo < masses_snap) & (masses_snap <= hi) & (metals_snap > 0.0)
                if any(idx):
                    metals_percs_nonzero[kk, snap, jj, :] = np.percentile(metals_snap[idx], 100*PERCS)
                    metals_percs_nonzero_num[kk, snap, jj] = np.sum(idx)

    fname = "{}_subhalos_metals".format(sim_name)
    shape_str = np.array("(NUM_KEYS, num_snaps, NBINS_METALS, NBINS)")
    keys_str = np.array("(" + ", ".join(METALS_KEYS) + ")")
    np.savez(fname, metals=metals, metals_percs=metals_percs, metals_percs_num=metals_percs_num,
             metals_percs_nonzero=metals_percs_nonzero, metals_percs_nonzero_num=metals_percs_nonzero_num,
             redshifts=redshifts,
             edges_metals=edges_metals, edges_masses=edges_masses, edges_percs_masses=PERCS_MASS_BIN_EDGES,
             sigma=SIGMA, percs=PERCS,
             shape=shape_str, metal_keys=keys_str)
    print("Saved data to '{}'".format(os.path.abspath(fname)))

    print("Done after {}".format(datetime.now()-beg_all))

    return


def halos_metallicity_redshift(arepo_sim_dir):
    beg_all = datetime.now()
    arepo_output_dir = os.path.join(arepo_sim_dir, "output", "")
    sim_name = os.path.split(arepo_sim_dir)[-1]
    print("subhalo_metals.halos_metallicity_redshift() : {} - {}".format(sim_name, arepo_output_dir))
    
    redshifts = utils.load_snapshot_redshifts(arepo_output_dir)
    num_snaps = len(redshifts)
    print("Snaps: ", num_snaps)

    fields = ['GroupGasMetallicity', 'GroupMass']

    NBINS_MASSES = 80
    NBINS_METALS = 90
    edges_metals = np.logspace(-8, 0, NBINS_METALS+1)
    edges_masses = np.logspace(8, 16, NBINS_MASSES+1)
    
    SIGMA = [0.0, 1.0, 2.0, 3.0]
    PERCS = np.array(sorted(list(set(np.concatenate(zmath.sigma(SIGMA, boundaries=True))))))
    NUM_PERCS = len(PERCS)

    metals = np.zeros((num_snaps, NBINS_METALS, NBINS_MASSES))

    PERCS_MASS_BIN_EDGES = np.logspace(8, 16, 9)
    NUM_PMB = PERCS_MASS_BIN_EDGES.size - 1
    metals_percs = np.zeros((num_snaps, NUM_PMB, NUM_PERCS))
    metals_percs_num = np.zeros((num_snaps, NUM_PMB), dtype=int)
    metals_percs_nonzero = np.zeros((num_snaps, NUM_PMB, NUM_PERCS))
    metals_percs_nonzero_num = np.zeros((num_snaps, NUM_PMB), dtype=int)

    edges = [edges_metals, edges_masses]

    for snap in range(num_snaps):
        beg = datetime.now()
        print(snap, "  ------")

        halos = ill.groupcat.loadHalos(arepo_output_dir, snap, fields=fields)
        # print(halos.keys())

        masses_snap = halos['GroupMass'] * CONV_ILL_TO_SOL.MASS
        num_halo = len(masses_snap)
        print("\tLoaded {} halos after {}".format(num_halo, str(datetime.now()-beg)))

        metals_snap = halos['GroupGasMetallicity']
        print("\t\tminmax = {:.2e}, {:.2e}".format(np.min(metals_snap), np.max(metals_snap)))
        metals[snap, :, :], xe, ye = np.histogram2d(metals_snap, masses_snap, bins=edges)

        for jj in range(NUM_PMB):
            lo = PERCS_MASS_BIN_EDGES[jj]
            hi = PERCS_MASS_BIN_EDGES[jj+1]

            idx = (lo < masses_snap) & (masses_snap <= hi)
            if any(idx):
                metals_percs[snap, jj, :] = np.percentile(metals_snap[idx], 100*PERCS)
                metals_percs_num[snap, jj] = np.sum(idx)

            idx = (lo < masses_snap) & (masses_snap <= hi) & (metals_snap > 0.0)
            if any(idx):
                metals_percs_nonzero[snap, jj, :] = np.percentile(metals_snap[idx], 100*PERCS)
                metals_percs_nonzero_num[snap, jj] = np.sum(idx)

    fname = "{}_halos_metals".format(sim_name)
    shape_str = np.array("(num_snaps, NBINS_METALS, NBINS)")
    np.savez(fname, metals=metals, metals_percs=metals_percs, metals_percs_num=metals_percs_num,
             metals_percs_nonzero=metals_percs_nonzero, metals_percs_nonzero_num=metals_percs_nonzero_num,
             redshifts=redshifts,
             edges_metals=edges_metals, edges_masses=edges_masses, edges_percs_masses=PERCS_MASS_BIN_EDGES,
             sigma=SIGMA, percs=PERCS, shape=shape_str)
    print("Saved data to '{}'".format(os.path.abspath(fname)))

    print("Done after {}".format(datetime.now()-beg_all))

    return


def halos_stars_metallicity_redshift(arepo_sim_dir):
    beg_all = datetime.now()
    arepo_output_dir = os.path.join(arepo_sim_dir, "output", "")
    sim_name = os.path.split(arepo_sim_dir)[-1]
    print("subhalo_metals.halos_stars_metallicity_redshift() : {} - {}".format(sim_name, arepo_output_dir))
    
    redshifts = utils.load_snapshot_redshifts(arepo_output_dir)
    num_snaps = len(redshifts)
    print("Snaps: ", num_snaps)

    fields = ['GroupGasMetallicity', 'GroupMassType', 'GroupLenType']

    NBINS_MASSES = 70
    NBINS_METALS = 60
    edges_metals = np.logspace(-8, 0, NBINS_METALS+1)
    edges_masses = np.logspace(6, 12, NBINS_MASSES+1)
    
    SIGMA = [0.0, 1.0, 2.0]
    PERCS = np.array(sorted(list(set(np.concatenate(zmath.sigma(SIGMA, boundaries=True))))))
    NUM_PERCS = len(PERCS)

    metals = np.zeros((num_snaps, NBINS_METALS, NBINS_MASSES))

    PERCS_MASS_BIN_EDGES = np.logspace(6, 12, 7)
    NUM_PMB = PERCS_MASS_BIN_EDGES.size - 1
    metals_percs = np.zeros((num_snaps, NUM_PMB, NUM_PERCS))
    metals_percs_num = np.zeros((num_snaps, NUM_PMB), dtype=int)
    metals_percs_nonzero = np.zeros((num_snaps, NUM_PMB, NUM_PERCS))
    metals_percs_nonzero_num = np.zeros((num_snaps, NUM_PMB), dtype=int)

    edges = [edges_metals, edges_masses]

    for snap in range(num_snaps):
        beg = datetime.now()
        print(snap, "  ------")

        halos = ill.groupcat.loadHalos(arepo_output_dir, snap, fields=fields)
        # print(halos.keys())

        # masses_snap = halos['GroupMass'] * CONV_ILL_TO_SOL.MASS
        masses_stars = halos['GroupMassType'][:, PARTS.STAR] * CONV_ILL_TO_SOL.MASS
        num_halo = len(masses_stars)
        idx = (masses_stars > 0)
        num_halo_stars = np.sum(idx)
        print("\tLoaded {} halos after {}".format(num_halo, str(datetime.now()-beg)))
        print("\tWith stars: {:.1e}/{:.1e} = {:.4f}".format(
            num_halo_stars, num_halo, num_halo_stars/num_halo))

        masses_snap = masses_stars[idx]
        metals_snap = halos['GroupGasMetallicity'][idx]
        print("\t\tminmax = {:.2e}, {:.2e}".format(np.min(metals_snap), np.max(metals_snap)))
        metals[snap, :, :], xe, ye = np.histogram2d(metals_snap, masses_snap, bins=edges)

        for jj in range(NUM_PMB):
            lo = PERCS_MASS_BIN_EDGES[jj]
            hi = PERCS_MASS_BIN_EDGES[jj+1]

            idx = (lo < masses_snap) & (masses_snap <= hi)
            if any(idx):
                metals_percs[snap, jj, :] = np.percentile(metals_snap[idx], 100*PERCS)
                metals_percs_num[snap, jj] = np.sum(idx)

            idx = (lo < masses_snap) & (masses_snap <= hi) & (metals_snap > 0.0)
            if any(idx):
                metals_percs_nonzero[snap, jj, :] = np.percentile(metals_snap[idx], 100*PERCS)
                metals_percs_nonzero_num[snap, jj] = np.sum(idx)

    fname = "{}_halos_stars_metals".format(sim_name)
    shape_str = np.array("(num_snaps, NBINS_METALS, NBINS)")
    np.savez(fname, metals=metals, metals_percs=metals_percs, metals_percs_num=metals_percs_num,
             metals_percs_nonzero=metals_percs_nonzero, metals_percs_nonzero_num=metals_percs_nonzero_num,
             redshifts=redshifts,
             edges_metals=edges_metals, edges_masses=edges_masses, edges_percs_masses=PERCS_MASS_BIN_EDGES,
             sigma=SIGMA, percs=PERCS, shape=shape_str)
    print("Saved data to '{}'".format(os.path.abspath(fname)))

    print("Done after {}".format(datetime.now()-beg_all))

    return


if __name__ == "__main__":
    print(sys.argv)
    if len(sys.argv) > 1:
        arepo_sim_dir = sys.argv[1]
    else:
        raise RuntimeError("No directory provided!")

    # main(arepo_sim_dir)
    # subhalos_metallicity_redshift(arepo_sim_dir)
    # halos_metallicity_redshift(arepo_sim_dir)
    halos_stars_metallicity_redshift(arepo_sim_dir)
