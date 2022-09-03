import simread.readsubfHDF5 as subf
import simread.readsnapHDF5 as ws
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import torrey_cmf
tc = torrey_cmf.number_density()

little_h = 0.703

min_bin_mass = 1e6
max_bin_mass = 1e12
n_bins = 20

min_bin_array = 10.0**( (0.0 + np.arange( n_bins ) ) / ( 1.0*n_bins ) * (np.log10( max_bin_mass ) - np.log10( min_bin_mass ) ) + np.log10( min_bin_mass ) )
max_bin_array = 10.0**( (1.0 + np.arange( n_bins ) ) / ( 1.0*n_bins ) * (np.log10( max_bin_mass ) - np.log10( min_bin_mass ) ) + np.log10( min_bin_mass ) )
mid_bin_array = 10.0**( (0.0 + np.arange( n_bins ) ) / ( 1.0*n_bins ) * (np.log10( max_bin_mass ) - np.log10( min_bin_mass ) ) + np.log10( min_bin_mass ) )

NUM_COLORS = n_bins
cm        = plt.get_cmap('nipy_spectral')
cNorm     = colors.Normalize(vmin=0, vmax=n_bins)


fig,ax=plt.subplots( figsize=(6,5))


for run in [     'explicit_feedback_256', \
                 'explicit_feedback_256_soft', \
                 'explicit_feedback_512', \
                 'explicit_feedback_512_soft_amd']:
    snaplist = glob.glob( '../'+run_base+'/output/snap*')
    n_snaps = len(snaplist)
    print snaplist
    snaparray = range( len( snaplist ) ) 

    for snapnum in snaparray:
        header = ws.snapshot_header( '../'+run+'/output/snapdir_'+str(snapnum).zfill(3)+'/snapshot_'+str(snapnum).zfill(3)+'.0.hdf5' )
        cat = subf.subfind_catalog(  '../'+run+'/', snapnum, keysel=['SubhaloMassInRadType'] )

        if run == 'dm_runs' or run=='explicit_feedback_fiducial':
            volume = (header.boxsize/little_h)**3.0
        if run == 'Illustris-Dark-1':
            volume = (header.boxsize/little_h/1000.0)**3.0

        z_str = "{:.1f}".format( header.redshift )
        z_str = r"$\mathrm{z="+z_str+"}$"

        print header.boxsize

        masses = np.array( cat.SubhaloMassInRadType[:,4] ) * 1e10 / little_h
        #masses = np.array( cat.GroupMass[:] ) * 1e10 / little_h

        new_masses = np.sort( masses )[::-1]
        print new_masses
        
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
        if np.sum(plt_index) > 4:
            ax.plot( mid_bin_array[plt_index], halo_mass_func[plt_index], color=cm( snapnum /15.0 ), label=z_str )
       # ax.plot( mid_bin_array, illustris_dmmf, ls='--',  color=cm( snapnum /7.0 ) )
     
    ax.set_ylim([1e-5, 5e-1])   
    ax.set_yscale('log')
    ax.set_xscale('log')

    ax.set_ylabel(r'$\mathrm{Stellar\;Mass\;Function\;(counts/Mpc^3/dex)}$')
    ax.set_xlabel(r'$\mathrm{Stellar\;Mass\;(M_\odot)}$')

    fig.subplots_adjust( left=0.13, bottom=0.12, top=0.99, right=0.96)

    data = np.loadtxt( '/n/home01/ptorrey/ObservationalData/Baldry2008_GSMF/GSMF_0005_reform.ascii.dat')
    print data[:,0]
    print data[:,1]
    ax.plot( 10.0**data[:,0], 10.0**data[:,1], c='k', lw=3 , label=r'$\mathrm{Baldry}$'   )
    ax.legend( frameon=False)
        
    fig.savefig('./plots/'+run+'_stellar_mass_function.pdf')



