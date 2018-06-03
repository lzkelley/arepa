import illustris
import readsubfHDF5
import numpy as np
import matplotlib.pyplot as plt

#each quantity at z=0, 1, 2, 4
#-stellar mass function
#-mean or median stellar mass as function of halo mass, for central galaxies
#-mean or median cold gas fraction (m_cold/m_star) as function of stellar mass
#-mean or median stellar and gas phase metallicity as function of stellar mass
#-mean or median SFR as function of stellar mass

# ============================================== #
target_redshifts = np.array([0, 1, 2, 4])
run='Illustris-1'
basedir = '/n/ghernquist/Illustris/Runs/'
dir = basedir+run+'/output/'
little_h = 0.704
# ============================================== #
min_mass = 1e10
max_mass = 1e14
n_bins   = 12
min_halo_mass_bins = 10.0**( (0.0 + np.arange(n_bins) ) / ( 1.0 * n_bins ) * (np.log10(max_mass) - np.log10(min_mass)) + np.log10(min_mass) )
mid_halo_mass_bins = 10.0**( (0.5 + np.arange(n_bins) ) / ( 1.0 * n_bins ) * (np.log10(max_mass) - np.log10(min_mass)) + np.log10(min_mass) )
max_halo_mass_bins = 10.0**( (1.0 + np.arange(n_bins) ) / ( 1.0 * n_bins ) * (np.log10(max_mass) - np.log10(min_mass)) + np.log10(min_mass) )
log_bin_size = (np.log10(max_mass) - np.log10(min_mass))/(1.0*n_bins)
smhm = np.zeros( n_bins ) 
# ============================================== #
scalefactors = illustris.load_scalefactors()
redshifts = 1.0/scalefactors - 1.0
snapshots = np.arange(redshifts.shape[0])
# ============================================== #

fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(1, 1, 1)

for redshift in target_redshifts:
    diff = np.abs( redshift - redshifts )
    snap = snapshots[ diff == diff.min() ]
    cat = readsubfHDF5.subfind_catalog(dir, snap[0], keysel=['SubhaloMassType', 'SubhaloGrNr', 'GroupFirstSub', 'Group_M_Crit200'])
    
    stellar_masses = cat.SubhaloMassType[:,4] * 1e10 / little_h
    halo_masses    = cat.SubhaloMassType[:,1] * 1e10 / little_h
    m200c_halo_masses = cat.Group_M_Crit200   * 1e10 / little_h

    gr_first_sub = cat.GroupFirstSub					# the first subhalo ID for each group
    sub_grnr     = cat.SubhaloGrNr					# the group ID for each subhalo
    sub_nr 	 = np.arange( stellar_masses.shape[0] )			# the subhalo ID

    sub_first_in_group = gr_first_sub[   sub_grnr  ]			# the first subhalo for the group to which this subhalo belongs
    is_central = sub_first_in_group == sub_nr				# compare the  "" "" against this subhalo ID.  If matches, is a central

    stellar_masses = stellar_masses[is_central]
    halo_masses = m200c_halo_masses[  sub_grnr[is_central]     ]  

    output = open('./data/stellar_mass_halo_mass_z'+str(redshift)+'.txt', 'w')
    output.write('# stellar mass halo mass relation data  \n')
    output.write('# col1 = M200_c (in solar masses) \n')
    output.write('# col2 = stellar mass (in solar masses)\n')
    output.write('\n')

    for bin_index,this_min_val in enumerate(min_halo_mass_bins):
	this_max_val = max_halo_mass_bins[bin_index]
	if np.sum( np.array( (halo_masses > this_min_val) & (halo_masses < this_max_val)  ) ) > 0:
	    smhm[bin_index] = np.median( stellar_masses[ (halo_masses > this_min_val) & (halo_masses < this_max_val)  ]  )
            str1 = "{0:.5e}".format(mid_halo_mass_bins[bin_index])
            str2 = "{0:.5e}".format(smhm[bin_index])
            output.write(str1+"    "+str2+"\n")

    output.close()
    ax.plot(mid_halo_mass_bins, smhm / mid_halo_mass_bins )
 
ax.set_ylabel(r'M${}_*$/M${}_{200,c}$')
ax.set_xlabel(r'M${}_{200,c}$ (M${}_\odot$)')
   
ax.set_yscale('log')
ax.set_xscale('log')
fig.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15)
fig.savefig('./plots/smhm.pdf')


