import numpy as np
import glob 
import sys
import os 

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm     as mplcm
import matplotlib.colors as colors
from matplotlib.backends.backend_pdf import PdfPages

from cycler import cycler


from scipy.stats import binned_statistic

import h5py

from scipy import interpolate

fontsize=14

# should be able to get rid of this in place of 
datablocks = [ 	'time',
		'sep',
		'sfr',
		'sfr_sm',
		'sfr_sm2',	
		'gas_mass',
		'z',
		'molec_gas_mass',
		'atomic_gas_mass',
		'warm_gas_mass',
		'hot_gas_mass'
	]

class merger_plotter:
    def __init__(self, saved_data_path, min_run=0, max_run=None, run_list=None, orbits=False, masses=False, legend_fontsize=5, legend_ncol=2):
	# figure out what data is available.  Sort and label it.
        if orbits:
            all_files = np.array( glob.glob(saved_data_path+"/saved_data/*oGG*hdf5")  )	# old
            all_files = np.array( glob.glob(saved_data_path+"/saved_data/orbit_*.hdf5")  ) # new
#G2G3_e_orbit_1_make_sfr_saved_data.dat
        else:
    	    all_files = np.array( glob.glob(saved_data_path+"/saved_data/")  )
            all_files = np.array( glob.glob(saved_data_path+"/saved_data/mass_ratio_*.hdf5"))


	tag_nr = np.zeros( all_files.shape[0] )
	for index, folder in enumerate(all_files):
            if orbits:
                tmp = os.path.basename( folder )
                #tag_nr[index] = int( tmp[tmp.index('oGG')+4:-5]  )			# old
                tag_nr[index] = int( tmp[tmp.index('orbit_')+6:tmp.index('.hdf5')]  )	# new
            else:
                #tag_nr[index] = int( os.path.basename( folder )[4:-5] )
                print os.path.basename( folder )
                tag_nr[index] = int( os.path.basename( folder )[11:-5] )

        self.orbits=orbits
        self.masses=masses
        self.legend_fontsize = legend_fontsize
        self.legend_ncol = legend_ncol



	order = np.argsort( tag_nr )
	new_files = [0] * all_files.shape[0]
	self.merger_int   = [0] * all_files.shape[0]
	for index,folder in enumerate( order ):
	    new_files[index] = all_files[order[index]]
	    self.merger_int[index]   = tag_nr[order[index]]
 
	all_files = new_files
	
	# remove any mergers we do not which to load/plot
	if max_run == None:  max_run = order.shape[0]
	all_files = all_files[min_run:max_run]
	self.merger_int   = self.merger_int[min_run:max_run]

	# prepare to open and package data from hdf5 files
        all_files = np.array( all_files  )
        self.all_data = {}
        self.all_tags = np.chararray( all_files.shape[0], itemsize=15 )
        for index, file in enumerate(all_files):
            print file
            f = h5py.File(file)
            tag = os.path.basename( file )[:-5]
            self.all_tags[index] = tag

            for field_name in f['data']:
                self.all_data[field_name+'_'+tag] = f['data'][field_name].value
            self.all_data['z_'+tag] /= 0.0127
            self.all_data['z_'+tag] = np.log10( 0.0004 * self.all_data['z_'+tag] ) + 12

	    good_index = self.all_data['time_'+tag] > 0
	    for field_name in f['data']:
		if np.sum( self.all_data[field_name+'_'+tag].shape )  > np.sum( np.array( good_index).shape ):
		    self.all_data[field_name+'_'+tag] = self.all_data[field_name+'_'+tag][:,:,good_index]
                    for iii in range(2):
                        for jjj in range(5):
                            if field_name in ['sfr', 'sfr_sm', 'sfr_sm2', 'm_atomic_gas', 'm_hot_gas', 'm_molec_gas', 'm_warm_gas']:
                                self.all_data[field_name+'_'+tag][iii,jjj,:] = smooth( self.all_data[field_name+'_'+tag][iii,jjj,:], window_len=5 )
			    if field_name in ['z']:
                                self.all_data[field_name+'_'+tag][iii,jjj,:] = smooth( self.all_data[field_name+'_'+tag][iii,jjj,:], window_len=5 )
		else: 
		    self.all_data[field_name+'_'+tag] = self.all_data[field_name+'_'+tag][good_index]
            f.close()

        # prepare to open and package data from hdf5 control files
        if self.orbits:
	    control_files = np.array( glob.glob(saved_data_path+"/saved_data/G?.hdf5")  )
        elif self.masses:
            control_files = np.array(glob.glob(saved_data_path+"/saved_data/G?.hdf5"))
        else:
            print "WARNING:  No merger suite type specified"
            control_files = np.array(glob.glob(saved_data_path+"/saved_data/G?.hdf5"))
 
        control_files.sort()
        self.control_data = {}
        self.control_tags = np.chararray( control_files.shape[0], itemsize=10 )
        for index, file in enumerate(control_files):
            f = h5py.File(file)
            tag = os.path.basename( file )[:-5]
#	    print "\n\n\n"
            print file, tag
#            print "\n\n\n"
            self.control_tags[index] = tag
            for field_name in f['data']:
                self.control_data[field_name+'_'+tag] = f['data'][field_name].value

            good_index = self.control_data['time_'+tag] > 0
            for field_name in f['data']:
                if np.sum( self.control_data[field_name+'_'+tag].shape )  > np.sum( np.array( good_index).shape ):
                    self.control_data[field_name+'_'+tag] = self.control_data[field_name+'_'+tag][:,:,good_index]
		    for iii in range(2):
			for jjj in range(5):
			    if field_name in ['sfr', 'sfr_sm', 'sfr_sm2', 'm_atomic_gas', 'm_hot_gas', 'm_molec_gas', 'm_warm_gas']:
			        self.control_data[field_name+'_'+tag][iii,jjj,:] = smooth(  self.control_data[field_name+'_'+tag][iii,jjj,:], window_len=5  )
			    if field_name in ['z']:
                                self.control_data[field_name+'_'+tag][iii,jjj,:] = smooth(  self.control_data[field_name+'_'+tag][iii,jjj,:], window_len=5 )
                else:
                    self.control_data[field_name+'_'+tag] = self.control_data[field_name+'_'+tag][good_index]

#            self.control_data['delta_sep_'+tag]  = self.control_data['sep_'+tag]
            self.control_data['z_'+tag] /= 0.0127
            self.control_data['z_'+tag] = np.log10( 0.0004 * self.control_data['z_'+tag] ) + 12
            f.close()
#            print 'time_'+tag
#            print self.control_data['time_'+tag]
            #print 'm_atomic_gas_'+tag
            #print self.control_data['m_atomic_gas_'+tag]
#            print "\n\n\n"



	n_delta = 200	# number of points for "delta" values
        f = h5py.File(control_files[0])	# just open the first control.  Any file will do.

	for index, file in enumerate(all_files):	# for all mergers
            if self.orbits:
                control_tag_1 = 'G2'
                control_tag_2 = 'G3'
            else:
                if all_files.shape[0] > 16:
                    control_tag_1 = 'Control_G1'
                    control_tag_2 = 'Control_G2'
                else:
                    control_tag_1 = 'G'+str(int(np.floor((self.merger_int[index]-1.0)/4.0))+1)
                    control_tag_2 = 'G'+str(int(np.round((self.merger_int[index]-1.0) % 4))+1)

            print control_tag_1, control_tag_2
            tag = self.all_tags[index]
	    #n_common = np.min( [self.all_data['time_'     +tag].shape[0] , self.control_data['time_'+control_tag_1].shape[0], self.control_data['time_'+control_tag_2].shape[0] ] )
	    min_time = np.max( [ np.min(self.all_data['time_'     +tag]), np.min(self.control_data['time_'+control_tag_1]), np.min(self.control_data['time_'+control_tag_2])  ]   )
            max_time = np.min( [ np.max(self.all_data['time_'     +tag]), np.max(self.control_data['time_'+control_tag_1]), np.max(self.control_data['time_'+control_tag_2])  ]   )
	    int_time = np.linspace( min_time, max_time, num=n_delta )

	    for field_name in f['data']:
		if field_name in ['time', 'sep']:
                    func = interpolate.interp1d(self.all_data['time_'+tag] ,  
                                                   self.all_data[field_name+'_'+tag] )
	            self.all_data['delta_'+field_name+'_'+tag] = func( int_time ) 
		elif field_name in ['z']:
                    self.all_data['delta_'+field_name+'_'+tag] =  np.zeros( (2, 5, n_delta) )

		    for iii in range(2):
                        if iii==0: control_tag=control_tag_1
                        else:      control_tag=control_tag_2
			for jjj in range(5):
			    tmp = np.zeros( n_delta )
			    c_func = interpolate.interp1d( self.control_data['time_'+control_tag] , 
							   self.control_data[field_name+'_'+control_tag][iii,jjj,:] )
			    m_func = interpolate.interp1d( self.all_data['time_'+tag] ,
                                                           self.all_data[field_name+'_'+tag][iii,jjj,:] )
		            
                            self.all_data['delta_'+field_name+'_'+tag][iii,jjj,:] = m_func( int_time ) - c_func( int_time )

		elif field_name in ['bh_acc_rate']:
                    #print "skip"
		    dummy = 1
		else:
		    tmp = np.zeros( (2, 5, n_delta) )
		    self.all_data['delta_'+field_name+'_'+tag] = np.zeros( (2, 5, n_delta) )
		    for iii in range(2):
			if iii==0: control_tag=control_tag_1
		        else:      control_tag=control_tag_2
		        for jjj in range(5):
		            c1_func = interpolate.interp1d(self.control_data['time_'+control_tag][:] , 
		            			   self.control_data[field_name+'_'+control_tag][0,jjj,:] )
		            g1_func = interpolate.interp1d(  self.all_data['time_'+tag][:], 
		            				 self.all_data[field_name+'_'+tag][iii,jjj,:] )

			    #print self.control_data['time_'+control_tag][:].min(), self.control_data['time_'+control_tag][:].max()				
                            #print self.all_data['time_'+tag][:].min(), self.all_data['time_'+tag][:].max()
			    #print int_time.min(), int_time.max()
			    #print " "
#			    print int_time
		            self.all_data['delta_'+field_name+'_'+tag][iii,jjj,:] = g1_func( int_time )/c1_func( int_time )
                            if(0):
				print 'time_'+control_tag
				print field_name+'_'+control_tag
                                print self.control_data['time_'+control_tag][:]
				print self.control_data[field_name+'_'+control_tag][0,jjj,:]
                                print 'delta_'+field_name+'_'+tag+'\n\n\n'
                                print self.all_data['delta_'+field_name+'_'+tag][0,jjj,:]
                                sys.exit()


#		    if "sfr_sm2" in field_name:
#			print 'delta_'+field_name+'_'+tag
#			print control_tag_1, control_tag_2
			#print self.all_data['delta_'+field_name+'_'+tag]

	    self.all_data['molec_gas_frac_'+tag]	= self.all_data['m_molec_gas_'+tag] / self.all_data['mgas_' +tag]

	f.close()



    def plot_dsfr_vs_psep_contour( self, n_proj=4, min_index=0, max_index=None, savetag='', label=None, rad_bin=4,
				  min_sep=0, max_sep=150.0, min_delta=-1, max_delta=1, 
				  **kwargs):
	fig, ax = self._setup_val_vs_sep_fig( min_sep=min_sep, max_sep=max_sep, min_y=min_delta, max_y=max_delta, ylog=False, ylabel=r'$\Delta$ SFR')
        if max_index==None: max_index=self.all_tags.shape[0]

        plot_x, plot_y = self.loop_and_store( 'delta_sep_', 'delta_sfr_sm2_', min_index, max_index, n_proj, rad_bin)
	plot_y = np.log10(plot_y)

        ax.hexbin( plot_x ,plot_y, gridsize=(40,20), cmap='cubehelix_r',  linewidth=0, bins='log',  extent=(min_sep, max_sep, min_delta, max_delta))

        result, edges, number = binned_statistic( plot_x[ np.isfinite(plot_y)   ], plot_y[np.isfinite(plot_y) ], statistic='median')
        ax.plot( edges, np.zeros_like( edges ), c='k', lw=2, ls='-')
        ax.plot( edges[:-1], result, c='r', lw=2, ls='--', label='median')

        ax.legend(loc=4)

        if label != None:
            ax.text( 75, 7.0/8.0*max_delta, label, ha='center', va='center', fontsize=18 )
        fig.savefig('./plots/'+savetag+'dsfr_vs_psep.pdf')
        plt.close()


    def plot_dz_vs_psep_contour( self, n_proj=4, min_index=0, max_index=None, savetag='', label=None, rad_bin=4, off_nuclear=False,
                                          min_sep=0, max_sep=150.0, min_delta=-0.1, max_delta=0.1,
                                          **kwargs):
        fig, ax = self._setup_val_vs_sep_fig( min_sep=min_sep, max_sep=max_sep, min_y=min_delta, max_y=max_delta, ylog=False, ylabel=r'$\Delta$Log(O/H)')
        if max_index==None: max_index=self.all_tags.shape[0]

        plot_x, plot_y = self.loop_and_store( 'delta_sep_', 'delta_z_', min_index, max_index, n_proj, rad_bin)
        ax.hexbin( plot_x ,plot_y,  gridsize=(40,25),  cmap='cubehelix_r', linewidth=0, bins='log', extent=(min_sep, max_sep, min_delta, max_delta) )

        result, edges, number = binned_statistic( plot_x[ np.isfinite(plot_y)   ], plot_y[np.isfinite(plot_y) ], statistic='median')
        ax.plot( edges[:-1], np.zeros_like( edges[:-1] ), c='b', lw=2, ls='-')
        ax.plot( edges[:-1], result, c='r', lw=2, ls='--', label='median')

        ax.legend(loc=4)
        if label != None:
            ax.text( 75, 7.0/8.0 * max_delta, label, ha='center', va='center', fontsize=18 )
	if rad_bin==0:
            fig.savefig('./plots/'+savetag+'dnz_vs_psep.pdf')
        else:
            fig.savefig('./plots/'+savetag+'dz_vs_psep.pdf')

        plt.close()


    def plot_dgas_vs_psep_contour( self, n_proj=4, min_index=0, max_index=None, savetag='', label=None, rad_bin=4, off_nuclear=False,
                                          min_sep=0, max_sep=150.0, min_delta=-0.9, max_delta=0.9, type=None,
                                          **kwargs):
        fig, ax = self._setup_val_vs_sep_fig( min_sep=min_sep, max_sep=max_sep, min_y=min_delta, max_y=max_delta, ylog=False, ylabel=r'$\Delta \mathrm{Log(M}_{\mathrm{'+type+'}}\mathrm{)}$')
        if max_index==None: max_index=self.all_tags.shape[0]

	if type not in ['molec', 'atomic', 'hot', 'warm']: sys.exit()
        plot_x, plot_y = self.loop_and_store( 'delta_sep_', 'delta_m_'+type+'_gas_', min_index, max_index, n_proj, rad_bin)
	plot_y = np.log10(plot_y) 

        ax.hexbin( plot_x ,plot_y,  gridsize=(40,25),  cmap='cubehelix_r', linewidth=0, bins='log', extent=(min_sep, max_sep, min_delta, max_delta) )

        result, edges, number = binned_statistic( plot_x[ np.isfinite(plot_y)   ], plot_y[np.isfinite(plot_y) ], statistic='median')
        ax.plot( edges[:-1], np.zeros_like( edges[:-1] ), c='b', lw=2, ls='-')
        ax.plot( edges[:-1], result, c='r', lw=2, ls='--', label='median')

        ax.legend(loc=4)
        if label != None:
            ax.text( 75, 7.0/8.0 * max_delta, label, ha='center', va='center', fontsize=18 )
        fig.savefig('./plots/'+savetag+'d'+type+'_vs_psep.pdf')
        plt.close()


    def loop_and_store( self, x_tag, y_tag, min_index, max_index, n_proj, rad_bin):
        count = 0
	for index,tag in enumerate(self.all_tags[min_index:max_index]):
	    if tag != "control":
	        count += self.all_data[x_tag+tag].shape[0] * 2 

	plot_x = np.zeros( count * n_proj )
	plot_y = np.zeros( count * n_proj )

	count = 0
	for index,tag in enumerate(self.all_tags[min_index:max_index]):
	    if tag != "control" and tag != "Control_G0" and tag != "Control_G1":
	      for iii in np.arange( n_proj ):
		phi = np.arcsin( np.random.uniform( size = self.all_data[x_tag      +tag].shape[0] ) )
                plot_x[count : count + self.all_data[x_tag      +tag].shape[0] ] = self.all_data[x_tag+tag]  * np.abs( np.cos( phi ) )
                plot_y[count : count + self.all_data[x_tag      +tag].shape[0] ] = self.all_data[y_tag+tag][0, rad_bin, :]
                count += self.all_data[x_tag      +tag].shape[0]

                plot_x[count : count + self.all_data[x_tag      +tag].shape[0] ] = self.all_data[x_tag+tag]  * np.abs( np.cos( phi ) )
                plot_y[count : count + self.all_data[x_tag      +tag].shape[0] ] = self.all_data[y_tag+tag][1, rad_bin, :]
                count += self.all_data[x_tag      +tag].shape[0]

	cull_index = plot_x > 1.0
	plot_x = plot_x [ cull_index] 
	plot_y = plot_y [ cull_index]

        return plot_x, plot_y



    def plot_all_sep_vs_time(self, name='all_sep_vs_time.pdf', **kwargs):
        self.x_axis_tag = 'time_'
        self.y_axis_tag = 'sep_'
        fig, ax = self._setup_val_vs_time_fig( ylog=False, min_y=0, max_y=200, ylabel=r'Sep (kpc)', **kwargs)
	self.loop_and_plot( ax, **kwargs )
        fig.savefig('./plots/'+name) 

    def plot_all_dsfr_vs_time(self, name='all_dsfr_vs_time.pdf',  rad_bin=0, **kwargs):
        self.x_axis_tag = 'delta_time_'
        self.y_axis_tag = 'delta_sfr_sm2_'
        fig, ax = self._setup_val_vs_time_fig( ylog=False, min_y=-1, max_y=1, ylabel=r'Log($\Delta$ SFR)', **kwargs)
        self.loop_and_plot( ax, galaxy_to_plot='gal1', take_y_log=True, rad_bin=4, **kwargs )
        self.loop_and_plot( ax, galaxy_to_plot='gal2', take_y_log=True, rad_bin=4, **kwargs )
        fig.savefig('./plots/'+name)

    def plot_all_sfr_vs_time(self, rad_bin=0, name='all_sfr_vs_time.pdf', min_y=1e-1, max_y=1e1, **kwargs):
        self.x_axis_tag = 'time_'
        self.y_axis_tag = 'sfr_sm2_'
        fig, ax = self._setup_val_vs_time_fig( ylog=True, min_y=min_y, max_y=max_y, ylabel=r'SFR ($M_\odot$)', **kwargs)

	self.loop_and_plot( ax, galaxy_to_plot='gal1', rad_bin=4, **kwargs) 
        self.loop_and_plot( ax, galaxy_to_plot='gal2', rad_bin=4, **kwargs)
        #self.loop_and_plot( ax, galaxy_to_plot='control', rad_bin=0, **kwargs)
        #self.loop_and_plot( ax, galaxy_to_plot='control', rad_bin=1, **kwargs)
        #self.loop_and_plot( ax, galaxy_to_plot='control', rad_bin=2, **kwargs)
        #self.loop_and_plot( ax, galaxy_to_plot='control', rad_bin=3, **kwargs)

        self.loop_and_plot( ax, galaxy_to_plot='control', rad_bin=4, **kwargs)

        fig.savefig('./plots/'+name)

    def plot_all_nz_vs_time(self, rad_bin=0, name='all_nz_vs_time.pdf',  **kwargs):
        self.x_axis_tag = 'time_'
        self.y_axis_tag = 'z_'
        fig, ax = self._setup_val_vs_time_fig( ylog=False, min_y=8.4, max_y=9.4, ylabel=r'Log(O/H)+12', **kwargs)
	self.loop_and_plot( ax, galaxy_to_plot='gal1',    rad_bin=0, **kwargs) 
        self.loop_and_plot( ax, galaxy_to_plot='gal2',    rad_bin=0, **kwargs)
        self.loop_and_plot( ax, galaxy_to_plot='control', rad_bin=0, **kwargs)


        fig.savefig('./plots/'+name)

    def plot_all_dnz_vs_time(self, rad_bin=0,  name='all_dnz_vs_time.pdf', **kwargs):
        self.x_axis_tag = 'delta_time_'
        self.y_axis_tag = 'delta_z_'
        fig, ax = self._setup_val_vs_time_fig( ylog=False, min_y=-0.2, max_y=0.1, ylabel=r'$\Delta$Log(O/H)', **kwargs)
	self.loop_and_plot( ax, galaxy_to_plot='gal1', rad_bin=0, **kwargs) 
        self.loop_and_plot( ax, galaxy_to_plot='gal2', rad_bin=0, **kwargs)
        fig.savefig('./plots/'+name)


    def plot_all_gas_vs_time(self, rad_bin=0, name='all_molec_gas_vs_time.pdf', type=None, min_y=1e6, max_y=2e8,  **kwargs):
	if type not in ['molec', 'atomic', 'warm', 'hot']: sys.exit()
	self.x_axis_tag = 'time_'
	self.y_axis_tag = 'm_'+type+'_gas_'
        fig, ax = self._setup_val_vs_time_fig( ylog=True, min_y=min_y, max_y=max_y, ylabel=r'M${}_{\mathrm{'+type+'}}$ ($M_\odot$)', **kwargs)
	self.loop_and_plot( ax, galaxy_to_plot='gal1',    rad_bin=4, y_factor=1e10,  **kwargs) 
        self.loop_and_plot( ax, galaxy_to_plot='gal2',    rad_bin=4, y_factor=1e10,  **kwargs)
        self.loop_and_plot( ax, galaxy_to_plot='control', rad_bin=4, y_factor=1e10,  **kwargs)
        fig.savefig('./plots/'+name)


    def plot_all_dgas_vs_time(self, rad_bin=0,  name='all_dmolec_gas_vs_time.pdf', type=None, min_y=-0.9, max_y=0.3, **kwargs):
        if type not in ['molec', 'atomic', 'warm', 'hot']: sys.exit()
        self.x_axis_tag = 'delta_time_'
        self.y_axis_tag = 'delta_m_'+type+'_gas_'
        fig, ax = self._setup_val_vs_time_fig( ylog=False, min_y=min_y, max_y=max_y, ylabel=r'Log($\Delta$M${}_{\mathrm{'+type+'}}$/$M_\odot$)', **kwargs)
        self.loop_and_plot( ax, galaxy_to_plot='gal1', rad_bin=4,  take_y_log=True,  **kwargs)
        self.loop_and_plot( ax, galaxy_to_plot='gal2', rad_bin=4,  take_y_log=True,  **kwargs)
        fig.savefig('./plots/'+name)



    def plot_all_molec_frac_vs_time(self, rad_bin=0, name='all_molec_frac_vs_time.pdf',  **kwargs):
        self.x_axis_tag = 'time_'
        self.y_axis_tag = 'molec_gas_frac_'
        fig, ax = self._setup_val_vs_time_fig( ylog=True, ylabel=r'Molec Gas Frac', **kwargs)
	self.loop_and_plot( ax, galaxy_to_plot='gal1', rad_bin=4, **kwargs)
        self.loop_and_plot( ax, galaxy_to_plot='gal2', rad_bin=4, **kwargs)
	fig.savefig('./plots/'+name)






    def plot_all_dsfr_vs_sep(self,  name='dsfr_vs_sep.pdf', **kwargs):
	self.x_axis_tag = 'delta_sep_'
	self.y_axis_tag = 'delta_sfr_sm2_'
        fig, ax = self._setup_val_vs_sep_fig( ylog=False, ylabel=r'Log($\Delta$ SFR)', **kwargs)
	self.loop_and_plot( ax, galaxy_to_plot='gal1', rad_bin=4, take_y_log=True, **kwargs) 
        self.loop_and_plot( ax, galaxy_to_plot='gal2', rad_bin=4, take_y_log=True, **kwargs)                        

	fig.savefig('./plots/'+name)

    def plot_all_dnz_vs_sep(self, rad_bin=0,  name='all_dnz_vs_sep.pdf', **kwargs):
        self.x_axis_tag = 'delta_sep_'
        self.y_axis_tag = 'delta_z_'
        fig, ax = self._setup_val_vs_sep_fig( ylog=False, ylabel=r'$\Delta$Log(O/H)', **kwargs)
        self.loop_and_plot( ax, galaxy_to_plot='gal1', rad_bin=0, **kwargs)
        self.loop_and_plot( ax, galaxy_to_plot='gal2', rad_bin=0, **kwargs)

        fig.savefig('./plots/'+name)


    def plot_all_dgas_vs_sep(self, rad_bin=0,  name='all_dmolec_gas_vs_sep.pdf', type=None, min_y=-0.9, max_y=0.3,  **kwargs):
        if type not in ['molec', 'atomic', 'warm', 'hot']: sys.exit()
        self.x_axis_tag = 'delta_sep_'
        self.y_axis_tag = 'delta_m_'+type+'_gas_'
        fig, ax = self._setup_val_vs_sep_fig( ylog=False, min_y=min_y, max_y=max_y, ylabel=r'Log($\Delta$M${}_{\mathrm{'+type+'}}$/$M_\odot$)', **kwargs)
        self.loop_and_plot( ax, galaxy_to_plot='gal1', rad_bin=4, take_y_log=True, **kwargs)
        self.loop_and_plot( ax, galaxy_to_plot='gal2', rad_bin=4, take_y_log=True, **kwargs)
        fig.savefig('./plots/'+name)




    def loop_and_plot(self, ax,  min_index=0, max_index=None, merg_ind=None, label_ind=None, cmap='nipy_spectral', rad_bin=None, multiple_rad_bins=None, legend=False, 
				take_y_log=False, y_factor=None, smooth=False, galaxy_to_plot='sum', **kwargs ):
        if galaxy_to_plot not in ['gal1', 'gal2', 'both', 'sum', 'control']:
	    print "bad choice for galaxy_to_plot"

        cm        = plt.get_cmap(cmap)

        if max_index==None: max_index=self.all_tags.shape[0]
        if merg_ind==None:  merg_ind=np.arange(self.all_tags.shape[0]) + min_index
	if label_ind==None:  label_ind=merg_ind

	if galaxy_to_plot=='control':
            this_tag_list = self.control_tags
            this_data     = self.control_data
        else:
            this_tag_list = self.all_tags[min_index:max_index]
            this_data     = self.all_data

        for index,tag in enumerate(this_tag_list):		#self.all_tags[min_index:max_index]):
            this_x = this_data[self.x_axis_tag+tag]
            this_y = this_data[self.y_axis_tag+tag]

	    if this_y.ndim==3:
	        if galaxy_to_plot=='sum':
                    this_y = this_y[0,:,:] + this_y[1,:,:]
	        elif galaxy_to_plot=='gal1':
                    this_y = this_y[0,:,:]
                    #print "gal1"
		    #print this_y
	        elif galaxy_to_plot=='gal2':
                    this_y = this_y[1,:,:]
                    #print "gal2"
                    #print this_y
                elif galaxy_to_plot=='control':
                    this_y = this_y[0,:,:]

            if rad_bin != None:  this_y = this_y[rad_bin,:]

	    if y_factor:         this_y *= y_factor
	    if take_y_log:       this_y = np.log10( this_y )	
	    if smooth:           this_y = smooth(this_y)

            color_val=index/(1.0*len(this_tag_list) )		#(5.0 - (index % 5.0))/7.5 + 2.0/7.5
            alpha_val=(np.floor((index-1)/5.0)+5)  /10.0  

	    if index in merg_ind or galaxy_to_plot=='control':
                if galaxy_to_plot=='control':
                    ax.plot( this_x, this_y , lw=2, ls='-', c='k', alpha =1.0, label='control' )
                else:
                    labels=['1','2','4','6','8']
                    if galaxy_to_plot=='control':
                        ax.plot( this_x, this_y , lw=2, ls='-', c='k' )
                    else:
		        if index in label_ind:
                            ax.plot( this_x, this_y , lw=2, ls='--', c=cm(color_val), alpha = alpha_val,  label=tag)
#'b='+labels[int((index-1) % 5)] )
		        else:
                            ax.plot( this_x, this_y , lw=2, ls='--', c=cm(color_val), alpha = alpha_val )

	if legend != False:
	    if type(legend) is int:
	        ax.legend(loc=legend, prop={'size':self.legend_fontsize},  ncol=self.legend_ncol, frameon=False )
	    else:
		ax.legend(prop={'size':self.legend_fontsize}, ncol=self.legend_ncol, frameon=False)


#    def _setup_dsfr_vs_psep_fig(self, num=1):
#        fig = plt.figure(num, figsize=(5,5))
#        plt.clf()
#
#        ax = setup_generic_axis( fig, self.all_tags.shape[0] )
#        ax.set_xlim([0, 150.0])
#        ax.set_ylim([-1, 1])       #[1e-6, 3e-1])
#        ax.set_xlabel(r'Proj Sep (kpc)', fontsize=14)
#        ax.set_ylabel(r'$\Delta$ SFR', fontsize=14)
#
#        return fig,ax



    def _setup_val_vs_sep_fig(self, num=1, min_sep=0, max_sep=150.0, min_y=1e6, max_y=1e11, ylog=False, ylabel=None, **kwargs):
        fig = plt.figure(num, figsize=(5,5))
        plt.clf()
        ax = setup_generic_axis( fig, self.all_tags.shape[0] )
        ax.set_xlim([min_sep, max_sep])
        ax.set_ylim([min_y, max_y])
        ax.set_xlabel(r'Sep (kpc)', fontsize=14)
        if ylog:    ax.set_yscale('log')
        if ylabel != None:    ax.set_ylabel(ylabel, fontsize=14)
        return fig,ax


    def _setup_val_vs_time_fig(self, num=1, min_time=0, max_time=5.0, min_y=1e6, max_y=1e11, ylog=False, ylabel=None, **kwargs):
        fig = plt.figure(num, figsize=(5,5))
        plt.clf()
        ax = setup_generic_axis( fig, self.all_tags.shape[0] )
        ax.set_xlim([min_time, max_time])
        ax.set_ylim([min_y, max_y])
        ax.set_xlabel(r't (Gyrs)', fontsize=14)
        if ylog:    ax.set_yscale('log')
        if ylabel != None:    ax.set_ylabel(ylabel, fontsize=14)
        return fig,ax


def setup_generic_axis( fig, n_colors ):
    ax = fig.add_subplot(1,1,1)
    fig.subplots_adjust(left=0.19, right=0.97, top=0.97, bottom=0.12)
    for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(fontsize)

    cm        = plt.get_cmap('nipy_spectral')
    cNorm     = colors.Normalize(vmin=0, vmax=n_colors)
    scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)
    ax.set_prop_cycle( cycler('color', [scalarMap.to_rgba(i) for i in range(n_colors)] ) )

#['c', 'm', 'y', 'k']) )
#set_color_cycle([scalarMap.to_rgba(i) for i in range(n_colors)])

    return ax



def smooth(x, **kwargs):
    try:
	new_x = smooth_1d(x, **kwargs)
    except:
	new_x = np.zeros_like(x)
	for iii in np.arange(x.shape[1]):
	    new_x[:,iii] = smooth_1d(x[:,iii])

    return new_x


def smooth_1d(x, window_len=6, window='hanning', **kwargs):
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

    if (window_len % 2)==0:
	window_len += 1

    if np.isnan(x[0]): 
        x[0] = x[1]

    s = np.zeros( x.shape[0] + window_len - 1 )
    for iii in np.arange( (window_len-1)/2.0 ):
        s[iii] = x[0]
	s[-(iii+1)] = x[-1]
    start_index = (window_len-1)/2.0
    s[start_index : start_index + x.shape[0]] = x

    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y


