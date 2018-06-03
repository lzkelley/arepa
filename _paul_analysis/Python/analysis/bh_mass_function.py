import illustris
import readsubfHDF5
import numpy as np
import matplotlib.pyplot as plt
import time

#each quantity at z=0, 1, 2, 4
#-stellar mass function
#-mean or median stellar mass as function of halo mass, for central galaxies
#-mean or median cold gas fraction (m_cold/m_star) as function of stellar mass
#-mean or median stellar and gas phase metallicity as function of stellar mass
#-mean or median SFR as function of stellar mass

author = 'P. Torrey (ptorrey@mit.edu)'
date = (time.strftime("%d/%m/%Y"))

# ============================================== #
target_redshifts = np.array([0, 1, 2, 4])
run='Illustris-1'
basedir = '/n/ghernquist/Illustris/Runs/'
dir = basedir+run+'/output/'
little_h = 0.704
# ============================================== #
min_mass = 1e6
max_mass = 1e10
n_bins   = 20
min_stellar_mass_bins = 10.0**( (0.0 + np.arange(n_bins) ) / ( 1.0 * n_bins ) * (np.log10(max_mass) - np.log10(min_mass)) + np.log10(min_mass) )
mid_stellar_mass_bins = 10.0**( (0.5 + np.arange(n_bins) ) / ( 1.0 * n_bins ) * (np.log10(max_mass) - np.log10(min_mass)) + np.log10(min_mass) )
max_stellar_mass_bins = 10.0**( (1.0 + np.arange(n_bins) ) / ( 1.0 * n_bins ) * (np.log10(max_mass) - np.log10(min_mass)) + np.log10(min_mass) )
log_stellar_mass_bin_size = (np.log10(max_mass) - np.log10(min_mass))/(1.0*n_bins)
smf = np.zeros( n_bins ) 

# ============================================== #
scalefactors = illustris.load_scalefactors()
redshifts = 1.0/scalefactors - 1.0
snapshots = np.arange(redshifts.shape[0])
volume=75.0 ** 3 /(little_h * little_h * little_h)				#  volume in Mpc^3
# ============================================== #

fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(1, 1, 1)

for redshift in target_redshifts:
    diff = np.abs( redshift - redshifts )
    snap = snapshots[ diff == diff.min() ]
    cat = readsubfHDF5.subfind_catalog(dir, snap[0], keysel=['SubhaloMassType'])
    stellar_masses = cat.SubhaloMassType[:,5] * 1e10 / little_h


    output = open('./data/stellar_mass_function_z'+str(redshift)+'.txt', 'w')
    output.write('# bh mass function data.  '+author+'  '+date+'  \n')
    output.write('# col1 = bh mass (in solar masses) \n')
    output.write('# col2 = number density (in #/Mpc^3/dex)\n')

    output.write('\n')

    for bin_index,this_min_val in enumerate(min_stellar_mass_bins):
	this_max_val = max_stellar_mass_bins[bin_index]
	count = np.sum( np.array( (stellar_masses > this_min_val) & (stellar_masses < this_max_val)) )
	smf[bin_index] = ( 1.0 * count )/(volume * (np.log10( this_max_val )-np.log10( this_min_val )) )
	str1 = "{0:.5e}".format(mid_stellar_mass_bins[bin_index])
        str2 = "{0:.5e}".format(smf[bin_index])


        output.write(str1+"    "+str2+"  \n")

    ax.plot(mid_stellar_mass_bins, smf)
    
    output.close()

ax.set_ylabel(r'$\Phi$ (#/Mpc${}^3$/dex)')
ax.set_xlabel(r'M${}_{BH}$ (M${}_\odot$)')
ax.set_yscale('log')
ax.set_xscale('log')
fig.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15)
fig.savefig('./plots/bhmf.pdf')


