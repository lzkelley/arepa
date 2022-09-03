import simread.readsubfHDF5 as subf
import simread.readsnapHDF5 as ws
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import torrey_cmf
from matplotlib.colors import LogNorm
import glob

tc = torrey_cmf.number_density()

little_h = 0.703

min_bin_mass = 1e8
max_bin_mass = 1e13
n_bins = 15

min_bin_array = 10.0**( (0.0 + np.arange( n_bins ) ) / ( 1.0*n_bins ) * (np.log10( max_bin_mass ) - np.log10( min_bin_mass ) ) + np.log10( min_bin_mass ) )
max_bin_array = 10.0**( (1.0 + np.arange( n_bins ) ) / ( 1.0*n_bins ) * (np.log10( max_bin_mass ) - np.log10( min_bin_mass ) ) + np.log10( min_bin_mass ) )
mid_bin_array = 10.0**( (0.0 + np.arange( n_bins ) ) / ( 1.0*n_bins ) * (np.log10( max_bin_mass ) - np.log10( min_bin_mass ) ) + np.log10( min_bin_mass ) )

NUM_COLORS = n_bins
cm        = plt.get_cmap('nipy_spectral')
cNorm     = colors.Normalize(vmin=0, vmax=n_bins)

fig,ax=plt.subplots( figsize=(6,5))

for run in ['explicit_feedback_256', 'explicit_feedback_256_soft', 'explicit_feedback_512', 'explicit_feedback_512_soft_amd']:

    snaplist = glob.glob( '../'+run+'/output/snapdir_*' )
    print snaplist
    snaparray = range( len( snaplist) ) 
    print snaparray
#    sys.exit()
#    snaparray = range( 13 ) 

    for snapnum in snaparray:

        fig_,ax_=plt.subplots( figsize=(6,5))
        header = ws.snapshot_header( '../'+run+'/output/snapdir_'+str(snapnum).zfill(3)+'/snapshot_'+str(snapnum).zfill(3)+'.0.hdf5' )
        cat = subf.subfind_catalog(  '../'+run+'/', snapnum, keysel=['SubhaloMassInRadType', 'SubhaloMassType', 'SubhaloGrNr', 'GroupFirstSub' ] )

        subhalo_number = range( cat.SubhaloGrNr.shape[0] )

        is_central = subhalo_number == cat.GroupFirstSub[cat.SubhaloGrNr] #  cat.SubhaloGrNr == 


        stellar_mass = np.array( cat.SubhaloMassInRadType[is_central,4] ) * 1e10 / little_h
        halo_mass    = np.array( cat.SubhaloMassType[is_central,1] )  * 1e10 / little_h
        smhm = stellar_mass / halo_mass 

        z_str = "{:.1f}".format( header.redshift )
        z_str = r"$\mathrm{z="+z_str+"}$"

#        gal_count = np.zeros( n_bins )
#        halo_mass_func = np.zeros( n_bins )
        med_smhm  = np.zeros( n_bins )

        for bin_index in range( n_bins ):
            in_bin_bool = (halo_mass > min_bin_array[bin_index]) & (halo_mass < max_bin_array[bin_index]) & (stellar_mass > 0)
            print smhm[ in_bin_bool ]
            if np.sum(in_bin_bool) > 3:
                med_smhm[ bin_index ] = np.median( smhm[ in_bin_bool ] )


        plt_index = (med_smhm > 0) & np.isfinite( med_smhm )
        if np.sum(plt_index) > 3:
            ax.plot( mid_bin_array[plt_index], med_smhm[plt_index], color=cm( snapnum /15.0 ), label=z_str )
            ax_.plot( np.log10( mid_bin_array[plt_index] ), np.log10( med_smhm[plt_index]) , color=cm( snapnum /15.0 ), label=z_str )

        plt_index = (halo_mass > 0) & (smhm>0) & (np.isfinite( smhm) )
        ax_.hist2d( np.log10(  halo_mass[plt_index] ), np.log10( smhm[plt_index] ), norm=LogNorm() , bins=25, range=[[9,14],[-3,-0.5]], cmap='Greys' )

        data = np.loadtxt( '/n/home01/ptorrey/ObservationalData/Behroozi2012_SMHM/c_smmr_z0.10_red_all_smf_m1p1s1_bolshoi_fullcosmos_ms.dat')
        ax_.plot( data[:,0], data[:,1], c='k', lw=3 , label=r'$\mathrm{Behroozi\;SMHM}$'   )
        ax_.text( 9.5, -0.75, z_str)

        fig_.savefig( './plots/'+run+'_stellar_mass_halo_mass_'+str(snapnum).zfill(3)+'.pdf')
#    ax.set_ylim([1e-5, 5e-1])   
    ax.set_yscale('log')
    ax.set_xscale('log')

    ax.set_ylabel(r'$\mathrm{Stellar\;Mass\;Function\;(counts/Mpc^3/dex)}$')
    ax.set_xlabel(r'$\mathrm{Stellar\;Mass\;(M_\odot)}$')

    fig.subplots_adjust( left=0.13, bottom=0.12, top=0.99, right=0.96)

    data = np.loadtxt( '/n/home01/ptorrey/ObservationalData/Behroozi2012_SMHM/c_smmr_z0.10_red_all_smf_m1p1s1_bolshoi_fullcosmos_ms.dat')
    print data[:,0]
    print data[:,1]
    ax.plot( 10.0**data[:,0], 10.0**data[:,1], c='k', lw=3 , label=r'$\mathrm{Behroozi\;SMHM}$'   )
    ax.legend( frameon=False)
        
    fig.savefig('./plots/'+run+'_stellar_mass_halo_mass.pdf')



