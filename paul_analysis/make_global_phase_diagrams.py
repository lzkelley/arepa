import simread.readsnapHDF5 as ws
import glob
import units.springel_units as units
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import no_repo_illustris_python as illustris_python

scalefactors = np.loadtxt( '/n/home01/ptorrey/Runs/BugTests/output_list.txt')   #/n/home01/ptorrey/IllustrisTNGAuxFiles/IllustrisTNG_scalefactors.txt' )
scalefactors = scalefactors[:,0]
print scalefactors

rundir='/n/home01/ptorrey/Runs/ISM_Boxes/'
snapnr = 50
little_h = 0.6774

print_phase_fractions = True

for snapnr in [11]:
    this_scalefactor = scalefactors[snapnr]
    this_redshift    = 1.0 / this_scalefactor - 1.0
    for runname in ['explicit_feedback_256_soft']:       #[455, 910, 1820]:
        base=rundir+runname+'/'
        path=base+'output/snapdir_'+str(snapnr).zfill(3)+'/snapshot_'+str(snapnr).zfill(3)
    
        print "density loading"
        rho = np.array( ws.read_block(path, 'RHO ', 0)).astype( np.float32 ) / 1000.0 **3.0
        print "Internal Energy loading"
        u   = np.array( ws.read_block(path, 'U   ', 0)).astype( np.float32 )
        print "Electron Abundance loading"
        ne  = np.array( ws.read_block(path, 'NE  ', 0)).astype( np.float32 )
        print "Temperature Loading"
        t   = np.array( units.gas_code_to_temperature( u, ne)).astype( np.float32 )
        del u
        del ne
        print "Mass Loading"
        m   = np.array( ws.read_block(path, 'MASS', 0)).astype( np.float32 )
        ms  = np.array( ws.read_block(path, 'MASS', 4)).astype( np.float32 )
    
        rho = np.array( units.gas_code_to_cgs_density( rho )).astype( np.float32 ) * little_h * little_h / this_scalefactor**3
        print "Loading done.  Starting to plot."
    
        for type_index,type in enumerate(['gas']):  #, 'metals']):
            fig,ax = plt.subplots(figsize=(5,5))
            if type=='gas':
                str_label=r'$\mathrm{Gas\/\/Mass}$'
                weight = m
            elif type=='metals':
                str_label=r'$\mathrm{Metal\/\/Mass}$'
                weight = m * z
            z_str   ='{:.1f}'.format( this_redshift )
            z_label = r'$\mathrm{z='+z_str+'}$'


            crit_density = units.critical_density( 1.0 )    # crit density in units of cm^{-3}

            if print_phase_fractions:
                denominator = np.sum( weight ) # either sum of mass, or sum of metal mass
                hot_index = (np.log10( t ) > 7 ) & ( np.isfinite( np.log10( t) ) )
                whim_index= (np.log10( t ) < 7 ) & ( np.log10( t ) > 5 ) & ( np.isfinite( np.log10( t) ) )
                diff_index= (np.log10( t ) < 5 ) & ( rho < 1000.0*crit_density * 0.0486) & ( np.isfinite( np.log10( t ) ) ) & (np.isfinite( np.log10( rho ) ) )
                cond_index= (np.log10( t ) < 5 ) & ( rho > 1000.0*crit_density * 0.0486) & ( np.isfinite( np.log10( t ) ) ) & (np.isfinite( np.log10( rho ) ) )
    
                hot_frac = np.sum( weight[ hot_index ] ) / denominator
                whim_frac = np.sum( weight[ whim_index]) / denominator
                diff_frac = np.sum( weight[ diff_index]) / denominator
                cond_frac = np.sum( weight[ cond_index]) / denominator



            if type == 'metals':
                rho_gas_mass = np.sum( m ) * 1e10 / little_h / (75000.0/little_h)**3.0        # total mass in solar masses / kpc^3
                rho_gas_mass = rho_gas_mass * 1.99e33 / 3.086e21                                            # g / kpc^3
                rho_gas_mass = rho_gas_mass / (3.086e21 * 1.67e-24) / (3.086e21)                                    # g / cm^3

                rho_metal_mass = np.sum( m * z ) * 1e10 / little_h / (75000.0/little_h)**3.0        # total mass in solar masses / kpc^3
                rho_metal_mass = rho_metal_mass * 1.99e33 / 3.086e21                                           # g / kpc^3
                rho_metal_mass = rho_metal_mass / (3.086e21 * 1.67e-24) / (3.086e21)                                    # g / cm^3

                gas_mass_str="{:.3f}".format( rho_gas_mass / (  0.0486 *  crit_density)   )
                gas_mass_str=r"$\rho_{\mathrm{gas}}/\Omega_{\mathrm{b}}\rho_{\mathrm{crit}}="+gas_mass_str+"$"
                ax.text( -7.65, 4.5, gas_mass_str, fontsize=12  )
                print "gas mass density is {:.8f}".format(  rho_gas_mass / (  0.0486 *  crit_density)   )
    
                gas_mass_str="{:.5f}".format( rho_metal_mass / (  0.0486 *  crit_density)   )
                gas_mass_str=r"$\rho_{\mathrm{Z,gas}}/\Omega_{\mathrm{b}}\rho_{\mathrm{crit}}="+gas_mass_str+"$"
                ax.text( -7.65, 4.25, gas_mass_str, fontsize=12  )
    
                print "gas metal mass density is {:.8f}".format(  rho_metal_mass / (  0.0486 *  crit_density)   )
    
                print rho_gas_mass
                print rho_gas_mass / crit_density
    
    
                rho_star_mass = np.sum( ms ) * 1e10 / little_h / (75000.0/little_h)**3.0        # total mass in solar masses / kpc^3
                rho_star_mass = rho_star_mass * 1.99e33 / 3.086e21                                            # g / kpc^3
                rho_star_mass = rho_star_mass / (3.086e21 * 1.67e-24) / (3.086e21)                                    # g / cm^3
    
                gas_mass_str="{:.3f}".format( rho_star_mass / (  0.0486 *  crit_density)   )
                gas_mass_str=r"$\rho_{\mathrm{*}}/\Omega_{\mathrm{b}}\rho_{\mathrm{crit}}="+gas_mass_str+"$"
                ax.text( -7.65, 3.75, gas_mass_str, fontsize=12  )
    
                print "star mass density is {:.8f}".format(  rho_star_mass / (  0.0486 *  crit_density)   )
    
    
                rho_star_metal_mass = np.sum( ms * zs ) * 1e10 / little_h / (75000.0/little_h)**3.0        # total mass in solar masses / kpc^3
                rho_star_metal_mass = rho_star_metal_mass * 1.99e33 / 3.086e21                                           # g / kpc^3
                rho_star_metal_mass = rho_star_metal_mass / (3.086e21 * 1.67e-24) / (3.086e21)                                    # g / cm^3
    
                gas_mass_str="{:.5f}".format( rho_star_metal_mass / (  0.0486 *  crit_density)   )
                gas_mass_str=r"$\rho_{\mathrm{Z,*}}/\Omega_{\mathrm{b}}\rho_{\mathrm{crit}}="+gas_mass_str+"$"
                ax.text( -7.65, 3.5, gas_mass_str, fontsize=12  )
    
                print "star metal mass density is {:.8f}".format(  rho_star_metal_mass / (  0.0486 *  crit_density)   )
    
    
            if type=='gas':
                ax.hist2d( np.log10( rho), np.log10( t), bins=256, cmap='afmhot_r', norm=LogNorm(), weights=weight )
            elif type=='metals':
                vals, xe, ye, im = ax.hist2d( np.log10( rho), np.log10( t), bins=256, cmap='afmhot_r', norm=LogNorm(), weights=weight, cmin=0.0001 )

            ax.set_xlim(  [-8,6]  )
            ax.set_ylim(  [1, 8.5]  )
            ax.set_xlabel(r'$\mathrm{Log(\rho [cm^{-3}])}$')
            ax.set_ylabel(r'$\mathrm{Log(T [K])}$')
            ax.text(  -2, 8.15, str_label, ha='center', va='center'  )
            ax.text(  -6.0, 8.0, z_label,   ha='center', va='center'  )


            ax.plot( [np.log10(1000.0*crit_density * 0.0486), np.log10( 1000.0 * crit_density * 0.0486)], [-100, 5], ls='-', lw=1, c='k' )
            ax.plot( [-100, 100], [5, 5], ls='-', lw=1, c='k')
            ax.plot( [-100, 100], [7, 7], ls='-', lw=1, c='k')

            if print_phase_fractions:
                tmp_str = "{:.1f}".format(hot_frac * 100.0)
                ax.text( -1.25, 7.5, r"$\mathrm{Hot: "+tmp_str+"\%}$", ha='left', va='center'  )
    
                tmp_str = "{:.1f}".format(whim_frac * 100.0)
                ax.text( -1.25, 6, r"$\mathrm{WHIM: "+tmp_str+"\%}$", ha='left', va='center'  )
    
                tmp_str = "{:.1f}".format(cond_frac * 100.0)
                ax.text( -1.25, 1.5, r"$\mathrm{Cond: "+tmp_str+"\%}$", ha='left', va='center'  )
    
                tmp_str = "{:.1f}".format(diff_frac * 100.0)
                ax.text( -6.0, 1.5, r"$\mathrm{Diffuse: "+tmp_str+"\%}$", ha='center', va='center'  )
    
            fig.subplots_adjust(left=0.10, bottom=0.115, top=0.99, right=0.98)
            if type=='gas':
                fig.savefig('./plots/global_phase_diagram_snr_'+str(snapnr).zfill(3)+'_n'+str(runname)+'.pdf')
            else:
                fig.savefig('./plots/global_metals_phase_diagram_snr_'+str(snapnr).zfill(3)+'_n'+str(runname)+'.pdf')
    
            plt.close('all')



