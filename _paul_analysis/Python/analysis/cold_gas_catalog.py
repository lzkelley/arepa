import illustris
import readsubfHDF5
import readhaloHDF5
import numpy as np
import matplotlib.pyplot as plt
import tables 
from mpi4py import MPI
import time


# some mpi init commands
#=========================================================#
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
n_tasks = size
this_task = rank
#=========================================================#



# ============================================== #
target_redshifts = np.array([4, 2, 1, 0])
target_redshifts = target_redshifts[::-1]
run='Illustris-1'
basedir = '/n/ghernquist/Illustris/Runs/'
dir = basedir+run+'/output/'
little_h = 0.704
# ============================================== #
scalefactors = illustris.load_scalefactors()
redshifts = 1.0/scalefactors - 1.0
snapshots = np.arange(redshifts.shape[0])
# ============================================== #


start_time = time.time()              
end_time = time.time()
count = 1             
 
min_sf_density = 1000.0

for redshift in target_redshifts:
    diff = np.abs( redshift - redshifts )
    snap = snapshots[ diff == diff.min() ]
    if this_task == 0:
        cat = readsubfHDF5.subfind_catalog(dir, snap[0], keysel=['SubhaloMassType', 'SubhaloSFR'])
    else:
	cat = None
    cat = comm.bcast(cat, root=0)

    stellar_mass   = cat.SubhaloMassType[:,4] * 1e10 / little_h
    total_gas_mass = cat.SubhaloMassType[:,0] * 1e10 / little_h
    sfr = cat.SubhaloSFR    
    n_subs = cat.nsubs

    cold_gas_mass 		= np.zeros( n_subs )		# local array
    cold_gas_mass_fraction 	= np.zeros( n_subs )		# local array
    global_cold_gas_mass               = np.zeros( n_subs )		# global array
    global_cold_gas_mass_fraction      = np.zeros( n_subs )    	# global array

    readhaloHDF5.reset()

    for this_sub_num in np.arange(n_subs):
      if ((this_sub_num % n_tasks) == this_task):

	if ((this_sub_num % (n_tasks*100)) == this_task):
	    end_time = time.time()
	    delta_t = end_time - start_time
            print "processing subhalo "+str(this_sub_num)+" out of "+str(n_subs)+" on task"+str(this_task)+"  (avg t = "+str(delta_t/count)+" secs)"
	    start_time = time.time()
	    count = 1

	if sfr[this_sub_num] > 0:
	    count += 1
            gas_particle_masses     = readhaloHDF5.readhalo(dir, 'snap',snap[0], 'MASS', 0, -1, this_sub_num, run=run)
            gas_particle_densities  = readhaloHDF5.readhalo(dir, 'snap',snap[0], 'RHO ', 0, -1, this_sub_num, run=run)
	    gas_particle_sfr	    = readhaloHDF5.readhalo(dir, 'snap',snap[0], 'SFR ', 0, -1, this_sub_num, run=run)

	    if np.sum(gas_particle_sfr) > 0:
	        sf_gas_rho = gas_particle_densities[ gas_particle_sfr > 0 ]
		#cold_gas_mass[this_sub_num] = np.sum( gas_particle_masses[ gas_particle_densities > 0.000849366 ] ) * 1e10 / little_h
		cold_gas_mass[this_sub_num] = np.sum( gas_particle_masses[ gas_particle_sfr       > 0.0         ] ) * 1e10 / little_h
		if stellar_mass[this_sub_num] > 0:
		    cold_gas_mass_fraction[this_sub_num] = cold_gas_mass[this_sub_num] / stellar_mass[this_sub_num]
		else:
		    cold_gas_mass_fraction[this_sub_num] = 0

    comm.Barrier()
    comm.Allreduce(cold_gas_mass, 		global_cold_gas_mass, 			op=MPI.SUM)
    comm.Allreduce(cold_gas_mass_fraction,   	global_cold_gas_mass_fraction,		op=MPI.SUM)

    if this_task==0:
        f = tables.openFile("cold_gas_masses_z"+str(redshift)+".hdf5", mode="w")
        f.createArray(f.root, "ColdGasMass",             global_cold_gas_mass)
        f.createArray(f.root, "ColdGasMassFraction",     global_cold_gas_mass_fraction)
        f.close()



