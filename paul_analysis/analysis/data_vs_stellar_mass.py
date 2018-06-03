import illustris
import readsubfHDF5
import numpy as np
import matplotlib.pyplot as plt
import hdf5lib

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
min_mass = 1e8
max_mass = 2e12
n_bins   = 20
min_mass_bins = 10.0**( (0.0 + np.arange(n_bins) ) / ( 1.0 * n_bins ) * (np.log10(max_mass) - np.log10(min_mass)) + np.log10(min_mass) )
mid_mass_bins = 10.0**( (0.5 + np.arange(n_bins) ) / ( 1.0 * n_bins ) * (np.log10(max_mass) - np.log10(min_mass)) + np.log10(min_mass) )
max_mass_bins = 10.0**( (1.0 + np.arange(n_bins) ) / ( 1.0 * n_bins ) * (np.log10(max_mass) - np.log10(min_mass)) + np.log10(min_mass) )
log_bin_size = (np.log10(max_mass) - np.log10(min_mass))/(1.0*n_bins)
sfr_array 	= np.zeros( n_bins ) 
m_cold_array 	= np.zeros( n_bins )
z_gas_array	= np.zeros( n_bins )
z_gas_sfr_array = np.zeros( n_bins )
z_stars_array	= np.zeros( n_bins )
cgm_array       = np.zeros( n_bins )            # cold gas mass
cgmf_array      = np.zeros( n_bins )		# cold gas mass fraction
# ============================================== #
scalefactors = illustris.load_scalefactors()
redshifts = 1.0/scalefactors - 1.0
snapshots = np.arange(redshifts.shape[0])
# ============================================== #

fig_sfr 	= plt.figure(figsize=(5,5))
fig_gas_z 	= plt.figure(figsize=(5,5))
fig_star_z      = plt.figure(figsize=(5,5))
fig_cgmf	= plt.figure(figsize=(5,5))

ax_sfr 		= fig_sfr.add_subplot(1, 1, 1)
ax_gas_z	= fig_gas_z.add_subplot(1, 1, 1)
ax_star_z       = fig_star_z.add_subplot(1, 1, 1)
ax_cgmf		= fig_cgmf.add_subplot(1, 1, 1)

for redshift in target_redshifts:
    diff = np.abs( redshift - redshifts )
    snap = snapshots[ diff == diff.min() ]
    cat = readsubfHDF5.subfind_catalog(dir, snap[0], keysel=['SubhaloMassType', 'SubhaloSFR', 'SubhaloGasMetallicity', 'SubhaloGasMetallicitySfr'])
    
    stellar_masses 	= cat.SubhaloMassType[:,4] * 1e10 / little_h		# units of M_solar
    sfr			= cat.SubhaloSFR					# units of M_solar / yr
    gas_z		= cat.SubhaloGasMetallicity				# unitless metallicity
    gas_sfr_z		= cat.SubhaloGasMetallicitySfr				# unitless metallicity
    star_z		= readsubfHDF5.subhalo_stellar_metallicities(snap=snap[0])

    # read this guy in manually... 
    file='/n/home01/ptorrey/ReviewPlots/cold_gas_masses_z'+str(redshift)+'.hdf5'
    f=hdf5lib.OpenFile(file, mode='r')
    cold_gas_mass=np.array(hdf5lib.GetData(f, "ColdGasMass")[:])		# should already be in units of M_solar (per cold_gas_catalog.py)
    f.close()

    output_cgm = open('./data/cold_gas_mass_vs_stellar_mass_z'+str(redshift)+'.txt', 'w')
    output_cgm.write('# Dense gas mass vs stellar mass (note, I used "cold" to describe this data throughout,\n')
    output_cgm.write('#    but its really only a density cut!!! per Romeels email)  \n')
    output_cgm.write('# col1 = stellar mass (in solar masses) \n')
    output_cgm.write('# col2 = dense gas mass (in solar masses)\n')
    output_cgm.write('# col3 = dense gas mass / stellar mass (unitless) \n')
    output_cgm.write('\n')

    output_sfr = open('./data/sfr_vs_stellar_mass_z'+str(redshift)+'.txt', 'w')
    output_sfr.write('# SFMS data  \n')
    output_sfr.write('# col1 = stellar mass (in solar masses) \n')
    output_sfr.write('# col2 = star formation rate (in solar masses / year)\n')
    output_sfr.write('\n')

    output_gas_z = open('./data/gas_z_vs_stellar_mass_z'+str(redshift)+'.txt', 'w')
    output_gas_z.write('# gas phase mass metallicity relation data  \n')
    output_gas_z.write('# col1 = stellar mass (in solar masses) \n')
    output_gas_z.write('# col2 = avg metallicity of *all gas* in the subhalo (unitless) \n')
    output_gas_z.write('# col3 = avg metallicity of *star forming gas* in the subhalo (unitless)  \n')
    output_gas_z.write('\n')

    output_star_z = open('./data/stellar_z_vs_stellar_mass_z'+str(redshift)+'.txt', 'w')
    output_star_z.write('# stellar mass metallicity relation data  \n')
    output_star_z.write('# col1 = stellar mass (in solar masses) \n')
    output_star_z.write('# col2 = avg stellar metallicity of all stars in the subhalo \n')
    output_star_z.write('\n')

    for bin_index,this_min_val in enumerate(min_mass_bins):
	this_max_val = max_mass_bins[bin_index]
	if np.sum( np.array( (stellar_masses > this_min_val) & (stellar_masses < this_max_val)  ) ) > 0:
	    this_bin_indicies = (stellar_masses > this_min_val) & (stellar_masses < this_max_val)

	    sfr_array[bin_index] 	= np.median(    sfr[     (stellar_masses > this_min_val) & (stellar_masses < this_max_val)  ]  )
            str1 = "{0:.5e}".format(mid_mass_bins[bin_index])
            str2 = "{0:.5e}".format(sfr_array[bin_index])
            output_sfr.write(str1+"    "+str2+"\n")

	    if star_z != None:
	        z_stars_array[bin_index]    = np.median(  star_z[    (stellar_masses > this_min_val) & (stellar_masses < this_max_val)  ]  )
		str1 = "{0:.5e}".format(mid_mass_bins[bin_index])
            	str2 = "{0:.5e}".format(z_stars_array[bin_index])
            	output_star_z.write(str1+"    "+str2+"\n")

	    cgm_array[bin_index]  = np.median( cold_gas_mass[this_bin_indicies] )
	    cgmf_array[bin_index] = np.median( cold_gas_mass[this_bin_indicies]/stellar_masses[this_bin_indicies] )
	    str1 = "{0:.5e}".format(mid_mass_bins[bin_index])
	    str2 = "{0:.5e}".format(cgm_array[bin_index])
            str3 = "{0:.5e}".format(cgmf_array[bin_index])
	    output_cgm.write(str1+"    "+str2+"    "+str3+"\n")

        if np.sum( np.array( (stellar_masses > this_min_val) & (stellar_masses < this_max_val) & (sfr > 0) ) ) > 0:
            this_bin_indicies = (stellar_masses > this_min_val) & (stellar_masses < this_max_val) & (sfr > 0)

            z_gas_array[bin_index]      = np.median(  gas_z[     this_bin_indicies  ]  )
            z_gas_sfr_array[bin_index]  = np.median(  gas_sfr_z[ this_bin_indicies  ]  )
            str1 = "{0:.5e}".format(mid_mass_bins[bin_index])
            str2 = "{0:.5e}".format(z_gas_array[bin_index])
            str3 = "{0:.5e}".format(z_gas_sfr_array[bin_index])
            output_gas_z.write(str1+"    "+str2+"    "+str3+"\n")

    output_sfr.close()
    output_gas_z.close()
    output_star_z.close()
    output_cgm.close()

    ax_sfr.plot(mid_mass_bins, sfr_array )
    ax_gas_z.plot( mid_mass_bins, z_gas_array )    
    ax_cgmf.plot(mid_mass_bins,  cgmf_array)



ax_sfr.set_ylabel(r'SFR (M${}_\odot$/yr)')
ax_sfr.set_xlabel(r'M${}_{*}$ (M${}_\odot$)')
ax_sfr.set_yscale('log')
ax_sfr.set_xscale('log')
fig_sfr.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15)
fig_sfr.savefig('./plots/sfr.pdf')

ax_gas_z.set_ylabel(r'Z gas')
ax_gas_z.set_xlabel(r'M${}_{*}$ (M${}_\odot$)')
ax_gas_z.set_yscale('log')
ax_gas_z.set_xscale('log')
fig_gas_z.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15)
fig_gas_z.savefig('./plots/gas_z.pdf')

ax_cgmf.set_ylabel(r'Cold Gas Mass Fraction')
ax_cgmf.set_xlabel(r'M${}_{*}$ (M${}_\odot$)')
ax_cgmf.set_yscale('log')
ax_cgmf.set_xscale('log')
fig_cgmf.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15)
fig_cgmf.savefig('./plots/cgmf.pdf')

