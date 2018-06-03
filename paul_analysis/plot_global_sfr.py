import numpy as np

import matplotlib.pyplot as plt


fig,ax=plt.subplots(figsize=(5,4))

ls_array = ['-', '--', ':', ':', ':', ':', ':' ]
label_list = [r'$\mathrm{256}$',       r'$\mathrm{512}$' , 
              r'$\mathrm{256\;Soft}$', r'$\mathrm{512\;Soft}$', \
              r'$\mathrm{256\;Soft\;Old\;Dens}$', \
              r'$\mathrm{256\;Metal\;Opacity}$', \
              r'$\mathrm{256\;Metal\;Opacity\;Softer}$', \
              ]

#adrwxr-xr-x 3 ptorrey hernquist_lab 553 Apr 18 16:22 explicit_feedback_512
#drwxr-xr-x 3 ptorrey hernquist_lab 301 Apr 19 01:37 explicit_feedback_256
#drwxr-xr-x 3 ptorrey hernquist_lab 253 Apr 19 09:28 explicit_feedback_512_speed_test
#drwxr-xr-x 3 ptorrey hernquist_lab 397 Apr 19 09:35 explicit_feedback_512_soft
#drwxr-xr-x 3 ptorrey hernquist_lab 239 Apr 19 16:37 explicit_feedback_256_soft
#drwxr-xr-x 3 ptorrey hernquist_lab 275 Apr 19 16:42 explicit_feedback_512_soft_amd
#drwxr-xr-x 3 ptorrey hernquist_lab 239 Apr 19 16:58 explicit_feedback_256_soft_old_dens


for run_index,run in enumerate(['explicit_feedback_256', 'explicit_feedback_512', 'explicit_feedback_256_soft', \
                                'explicit_feedback_512_soft_amd', 'explicit_feedback_256_soft_old_dens', 'explicit_feedback_256_tau_normal', 'explicit_feedback_256_tau_normal_softer']):
    data = np.loadtxt( '../'+run+'/output/sfr.txt')
    redshift = 1.0 / ( data[:,0] ) - 1.0
    csfrd = data[:,2] / ( (10.0 / 0.704)**3.0  )

    ax.plot(np.log10( 1.0 + redshift), csfrd, ls=ls_array[run_index], label=label_list[run_index] )
    print run
    print np.log10( 1.0 + redshift )
    print csfrd
    print np.max(csfrd)
    print " "

b_data = np.loadtxt( '/n/home01/ptorrey/ObservationalData/Behroozi2012_CSFR/csfrs_new.dat' )

ax.scatter( np.log10( b_data[:,0] + 1.0), 10.0**b_data[:,1], color='k', s=1.0,label=r'$\mathrm{Behroozi}$' )


ax.legend(frameon=False, loc=3, fontsize=14)


ax.set_xlim([0,1.3])
ax.set_ylim([1e-5, 1e0])        #[1e-4, 1e0])

ax.set_xlabel(r'$\mathrm{Log(1+}z\mathrm{)}$', fontsize=14)
ax.set_ylabel(r'$\mathrm{Star\;Formation\;Rate\;(M_\odot/yr/Mpc^3)}$', fontsize=14)
ax.set_yscale('log')

fig.subplots_adjust(left=0.16, bottom=0.14, top=0.965, right=0.96 )

fig.savefig('cosmic_sfr.pdf')



