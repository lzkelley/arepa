"""Create a snapshot_times file for input to Arepo.
"""

import numpy as np
import scipy as sp
import astropy as ap

import zcode.math as zmath
import zcode.astro as zastro
from zcode.constants import *

cosmo = ap.cosmology.WMAP9

# FNAME = "snapshot_output_list.txt"
# FNAME = "snapshot_output_list_high-z.txt"
FNAME = "snapshot_output_list_hisnap.txt"
FLAG = 1
PAD = 1e-6


def main():
    # bounds = [20, 12, 8, 6, 4.0, 3.0]
    # steps = [-4, -2, -1, -0.5, -0.25, -0.1]
    bounds = [19, 13, 12, 9, 7.0]
    steps = [-1, -0.25, -0.2, -0.04]

    zero = (len(bounds) == len(steps))

    redz = []
    for ii in range(len(bounds)):
        last = (ii == len(bounds)-1)
        if ii == len(steps):
            break
        lo = bounds[ii]
        if zero:
            hi = -PAD if last else bounds[ii+1]
        else:
            hi = bounds[ii+1]

        zz = np.arange(lo, hi, steps[ii])
        print(ii, lo, hi)
        print("\t[{:5.2f}, {:5.2f}] - {:2d}".format(zz[0], zz[-1], zz.size))
        redz.append(zz)

    if zero:
        redz[-1] = 0.0
    else:
        redz.append([bounds[-1]])
        
    redz = np.concatenate(redz)

    print(len(redz), "total")
    print(redz)

    age = cosmo.lookback_time(redz).to('Gyr').value
    dt = -np.diff(age)

    with open(FNAME, 'w') as out:
        for ii in range(len(redz)):
            sc = 1/(redz[ii] + 1)
            _dt = dt[ii-1] if ii > 0 else 0.0
            vals = [sc, FLAG, ii, redz[ii], age[ii], _dt]
            vals = '{:11.8f}  {:1d}     # {:4d} {:8.5f} {:8.5f} {:8.5f}'.format(*vals)
            out.write(vals + "\n")
            print("{:3d}".format(ii), vals)

    return


if __name__ == "__main__":
    main()
