import simread.readsubfHDF5 as subf
import simread.readsnapHDF5 as ws
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import torrey_cmf
tc = torrey_cmf.number_density()

little_h = 0.703

min_bin_mass = 1e7
max_bin_mass = 1e13
n_bins = 30

min_bin_array = 10.0**( (0.0 + np.arange( n_bins ) ) / ( 1.0*n_bins ) * (np.log10( max_bin_mass ) - np.log10( min_bin_mass ) ) + np.log10( min_bin_mass ) )
max_bin_array = 10.0**( (1.0 + np.arange( n_bins ) ) / ( 1.0*n_bins ) * (np.log10( max_bin_mass ) - np.log10( min_bin_mass ) ) + np.log10( min_bin_mass ) )
mid_bin_array = 10.0**( (0.0 + np.arange( n_bins ) ) / ( 1.0*n_bins ) * (np.log10( max_bin_mass ) - np.log10( min_bin_mass ) ) + np.log10( min_bin_mass ) )

NUM_COLORS = n_bins
cm        = plt.get_cmap('nipy_spectral')
cNorm     = colors.Normalize(vmin=0, vmax=n_bins)


fig,ax=plt.subplots()
for run in ['explicit_feedback_fiducial']:      #'dm_runs', 'Illustris-Dark-1']:

    if run=='dm_runs': snaparray = np.array([6])	#range(7)
    if run=='Illustris-Dark-1': snaparray = np.array([135])
    if run=='explicit_feedback_fiducial': snaparray = range(12)

    for snapnum in snaparray:
        header = ws.snapshot_header( '../'+run+'/output/snapdir_'+str(snapnum).zfill(3)+'/snapshot_'+str(snapnum).zfill(3)+'.0.hdf5' )
        cat = subf.subfind_catalog(  '../'+run+'/', snapnum, keysel=['Group_M_Crit200', 'GroupMass'] )

        if run == 'dm_runs' or run=='explicit_feedback_fiducial':
            volume = (header.boxsize/little_h)**3.0
        if run == 'Illustris-Dark-1':
            volume = (header.boxsize/little_h/1000.0)**3.0

        print header.boxsize

        masses = np.array( cat.Group_M_Crit200[:] ) * 1e10 / little_h
        #masses = np.array( cat.GroupMass[:] ) * 1e10 / little_h
        
        gal_count = np.zeros( n_bins )
        halo_mass_func = np.zeros( n_bins )
        
        for bin_index in range( n_bins ):
            in_bin_bool = (masses > min_bin_array[bin_index]) & (masses < max_bin_array[bin_index])
            bin_size_dex = np.log10( max_bin_array[bin_index] ) - np.log10( min_bin_array[bin_index] )
            halo_mass_func[ bin_index ] = np.sum(in_bin_bool)*1.0 / (volume * bin_size_dex) 
        
        illustris_dmmf = np.zeros( n_bins )
        z=header.redshift
        for bin_index in range(n_bins):
            illustris_dmmf[bin_index] = 10.0**tc.dm_mf_fit( np.log10( mid_bin_array[bin_index] ), z) 
            
        plt_index = halo_mass_func > 0
        ax.plot( mid_bin_array[plt_index], halo_mass_func[plt_index], color=cm( snapnum /7.0 ) )
       # ax.plot( mid_bin_array, illustris_dmmf, ls='--',  color=cm( snapnum /7.0 ) )
     
    ax.set_ylim([1e-4, 1e2])   
    ax.set_yscale('log')
    ax.set_xscale('log')
        
    fig.savefig('./plots/'+run+'_dm_mass_function.png')



