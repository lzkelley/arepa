#!/usr/bin/env python
""" routine for making hybrid nuclear disk ICs    """

__author__ = "Paul Torrey"
__copyright__ = "Copyright 2015, The Authors"
__credits__ = ["Paul Torrey"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Paul Torrey"
__email__ = "ptorrey@mit.harvard.edu"
__status__ = "Public Release.  v1.0."


# generic module imports
import numpy as np
import sys

# my module imports
import simread.readsnapHDF5 as ws


ic1 = str(sys.argv[1])	# it's assumed these are hdf5 files
ic2 = str(sys.argv[2])

h1 = ws.snapshot_header(ic1)    # load header 1
npart1 =  h1.nall
h2 = ws.snapshot_header(ic2)    # load header 2
npart2 = h2.nall

nall = h1.nall + h2.nall
all_blocks1 = ws.return_tags(ic1)

f = ws.openfile('new_ics.hdf5')
my_header = h1
last_ID = 1
for parttype in np.arange(6):
    this_part_blocks = all_blocks1[parttype]

    for block in this_part_blocks: 
	print parttype, block
        data = None
        if npart1[parttype] > 0:  data1 = ws.read_block(ic1, block, parttype=parttype)
        else:  data1 = np.zeros( 1 )

        if npart2[parttype] > 0:  data2 = ws.read_block(ic1, block, parttype=parttype)
        else:  data2 = np.zeros( 1 ) 

        if npart1[parttype] > 0 and npart2[parttype] > 0:
	    if parttype==5:
	        data = data1
		my_header.nall[5] = 1
	    else:
                data = np.append(data1, data2, axis=0)
        elif npart1[parttype] > 0:  data = data1
        elif npart2[parttype] > 0:  data = data2
        else: print "[error]:  You shouldnt be here. expect an error." 

	if block=="ID  ":
	    data = np.arange( data.shape[0] ) + last_ID
	    last_ID += data.shape[0]


	ws.write_block(f, block, parttype, data)


my_header.NumPart_ThisFile = nall
my_header.NumPart_Total    = nall
my_header.NumPart_Total_HighWord = nall * 0
my_header.Time = 0.0

ws.writeheader(f, my_header)
ws.closefile(f)



