import no_repo_illustris_python as illustris_python
import numpy as np
from util import hdf5lib 
import h5py
import sys

scalefactor_list = np.loadtxt('/n/home01/ptorrey/IllustrisAuxFiles/Illustris_scalefactors.txt')
redshift_list    = 1.0/scalefactor_list - 1.0


class aux_catalogs:
    def __init__(self, base, run, snapnum, subbox=None, **kwargs):

#        self.run = 'Illustris-'+str(run)+'/'
        self.run = run
        self.base = base
        self.basePath = base+self.run+'/output/'
        print self.basePath
        self.boxsize=75000.0
	self.little_h = 0.704
        self.min_n_stars = 500
        self.min_n_gas = 500
        self.min_n_pca = 5000
        self.snapnum = snapnum
        self.subbox = subbox

        if 'TNG' in self.run:
            all_scalefactors = np.loadtxt('/n/home01/ptorrey/IllustrisTNGAuxFiles/IllustrisTNG_scalefactors.txt')
            self.little_h = 0.6774
        else:
	    all_scalefactors = np.loadtxt('/n/home01/ptorrey/IllustrisAuxFiles/Illustris_scalefactors.txt')
            self.little_h = 0.704

        if subbox!= None:
            all_scalefactors = np.loadtxt('/n/hernquistfs3/ptorrey/IllustrisTNG/SubboxTreeConstruction/walk_trees/subbox_redshift_list.txt')
            all_scalefactors = all_scalefactors[:,1]
            import simread.readhaloHDF5 as readhaloHDF5
            full_path = self.base+'/output_subbox'+str(self.subbox)
            self.halo_reader = readhaloHDF5.HaloReader( full_path, self.snapnum, snapbase='snap-groupordered' )


        self.scalefactor = all_scalefactors[snapnum]

        if subbox != None:
            full_path = base+'/output_subbox'+str(subbox)
        else:
            full_path = base+self.run+'/output'

        self.subcat = illustris_python.cat.loadSubhalos(full_path, snapnum)	#snapnum, run=self.run, fields=fields, base=base)
        self.nsubs  = self.subcat['count']

    def make_sz_aux(self, **kwargs):
        fields=[u'Coordinates', u'Masses', u'GFM_StellarFormationTime']
        partType=4
        stored_core_1     = np.zeros( self.nsubs ) - 1.0
        stored_core_2     = np.zeros( self.nsubs ) - 1.0

        for iii in np.arange(self.nsubs):
          if self.subcat['SubhaloLenType'][iii, 4] > self.min_n_stars:
            if iii % 1000 == 0:
                print "processing subhalo "+str(iii)+" out of "+str(self.nsubs)
            data =  illustris_python.snapshot.loadSubhalo(self.basePath, self.snapnum,iii,partType,fields=fields)
            for xyz in np.arange(3):
                data['Coordinates'][:,xyz] -= self.subcat['SubhaloPos'][iii,xyz]
                data['Coordinates'][ (data['Coordinates'][:,xyz] < -self.boxsize/2.0) ,xyz]  += self.boxsize
                data['Coordinates'][ (data['Coordinates'][:,xyz] >  self.boxsize/2.0) ,xyz]  -= self.boxsize

            r  = np.sqrt( data['Coordinates'][:,0]**2 + data['Coordinates'][:,1]**2 + data['Coordinates'][:,2]**2 ) * self.scalefactor / self.little_h
            age= np.array( data['GFM_StellarFormationTime'] )

            index = (age > 0)
            if np.sum(index) > 0:
                core_1_index = r < 1.0
                if np.sum(core_1_index)>0:
                    stored_core_1[iii] = np.sum(data['Masses'][core_1_index]) * 1e10 / self.little_h
                core_2_index = r < 2.0
                if np.sum(core_2_index)>0:
                    stored_core_2[iii] = np.sum(data['Masses'][core_2_index]) * 1e10 / self.little_h

        new_file_name = './catalogs/'+self.run+'/stellar_core_mass/stellar_core_mass_'+str(self.snapnum)+'.hdf5'
        print new_file_name
        f_write = hdf5lib.OpenFile(new_file_name, mode="w")
        subhalo_group = hdf5lib.CreateGroup(f_write, "Subhalo")
        hdf5lib.CreateArray(f_write, subhalo_group, 'StellarCoreMass1'  , stored_core_1)
        hdf5lib.CreateArray(f_write, subhalo_group, 'StellarCoreMass2'  , stored_core_2)
        f_write.close()

    def make_rotation_curve_aux(self, **kwargs):        
        fields=[u'Coordinates', u'Velocities', u'GFM_StellarFormationTime',u'Masses']
        stored_80_light_rad  = np.zeros( self.nsubs ) - 1.0 
        stored_v80_from_pot  = np.zeros( self.nsubs ) - 1.0
        stored_v80_mean_star = np.zeros( self.nsubs ) - 1.0 
	stored_v80_med_star  = np.zeros( self.nsubs ) - 1.0
        with h5py.File(illustris_python.snapshot.snapPath(self.basePath,self.snapnum),'r') as f:
            header = dict( f['Header'].attrs.items() )
#	print header[u'MassTable']

        for iii in np.arange(self.nsubs):
            if self.subcat['SubhaloLenType'][iii, 4] > 50:
                if iii % 100 == 0:  print "processing subhalo "+str(iii)+" out of "+str(self.nsubs)
		data = {}
		for partType in [0, 1, 4, 5]:
                    if self.subcat['SubhaloLenType'][iii,partType] > 0:
			if partType==1:
                            data[str(partType)] = illustris_python.snapshot.loadSubhalo(self.basePath, self.snapnum,iii,partType,fields=[u'Coordinates', u'Velocities'])
			else:
                            data[str(partType)] = illustris_python.snapshot.loadSubhalo(self.basePath, self.snapnum,iii,partType,fields=[u'Masses', u'Coordinates', u'Velocities'])

                        for xyz in np.arange(3):
                            data[str(partType)]['Coordinates'][:,xyz] -= self.subcat['SubhaloPos'][iii,xyz]
                            data[str(partType)]['Coordinates'][ (data[str(partType)]['Coordinates'][:,xyz] < -self.boxsize/2.0) ,xyz]  += self.boxsize
                            data[str(partType)]['Coordinates'][ (data[str(partType)]['Coordinates'][:,xyz] >  self.boxsize/2.0) ,xyz]  -= self.boxsize
                            data[str(partType)]['Velocities'][:,xyz] -= np.median( data[str(partType)]['Velocities'][:,xyz] )

		    

                r = np.sqrt( data['4']['Coordinates'][:,0]**2 + data['4']['Coordinates'][:,1]**2 + data['4']['Coordinates'][:,2]**2 )   
                order = np.argsort( r )
                r_ = r[order]
	        m_ = data['4']['Masses'][order]
		prof = np.abs(  np.cumsum( m_ ) / np.sum( m_ ) - 0.8 )
		index = np.where(prof == np.min(prof) )[0][0]
                stored_80_light_rad[iii] = r_[index]
	
	        menc = 0
		for partType in [0, 1, 4, 5]:
		    if self.subcat['SubhaloLenType'][iii,partType] > 0:
		        r = np.sqrt( data[str(partType)]['Coordinates'][:,0]**2 + data[str(partType)]['Coordinates'][:,1]**2 + data[str(partType)]['Coordinates'][:,2]**2 )
                        if partType==1:
                            menc += np.sum( r<stored_80_light_rad[iii] ) * header[u'MassTable'][1]
                        else:
                            menc += np.sum( data[str(partType)]['Masses'][r<stored_80_light_rad[iii]] )
		menc *= 1e10 / 0.704 * 1.989e30 			# mass enclosed in kg 
		renc = stored_80_light_rad[iii] / 0.704 * 3.086e19	# radius in m
		grav_const = 6.674e-11 					# m^3 / kg / s^2
		vcirc = np.sqrt( grav_const * menc / renc ) 
		stored_v80_from_pot[iii] = vcirc 

		r = np.sqrt( data['4']['Coordinates'][:,0]**2 + data['4']['Coordinates'][:,1]**2 + data['4']['Coordinates'][:,2]**2 ) + 1e-3
                v = np.sqrt( data['4']['Velocities'][:,0]**2  + data['4']['Velocities'][:,1]**2  + data['4']['Velocities'][:,2]**2 )
                vr_star =(data['4']['Coordinates'][:,0]*data['4']['Velocities'][:,0] + 
                          data['4']['Coordinates'][:,1]*data['4']['Velocities'][:,1] + 
                          data['4']['Coordinates'][:,2]*data['4']['Velocities'][:,2] )/r
		vphi_star = np.sqrt( v**2 - vr_star**2 )
		index = (r > 0.8 * stored_80_light_rad[iii]) & (r < 1.2 * stored_80_light_rad[iii])
		stored_v80_mean_star[iii] = np.mean( vphi_star[index] )
		stored_v80_med_star[iii]  = np.median( vphi_star[index] )


        new_file_name = './catalogs/'+self.run+'/rotation_curve/rotation_curve_'+str(self.snapnum)+'.hdf5'
        print new_file_name
        f_write = hdf5lib.OpenFile(new_file_name, mode="w")
        subhalo_group = hdf5lib.CreateGroup(f_write, "Subhalo")
        hdf5lib.CreateArray(f_write, subhalo_group, 'r80_from_mass'       , stored_80_light_rad )
        hdf5lib.CreateArray(f_write, subhalo_group, 'v80_from_pot'        , stored_v80_from_pot )
        hdf5lib.CreateArray(f_write, subhalo_group, 'v80_from_mean_star'  , stored_v80_mean_star)
        hdf5lib.CreateArray(f_write, subhalo_group, 'v80_from_median_star', stored_v80_med_star )

        f_write.close()
        

    def make_velocity_dispersion_aux(self, **kwargs):
        fields=[u'Coordinates', u'Velocities', u'GFM_StellarFormationTime']
        partType=4
        stored_vd_all        = np.zeros( self.nsubs ) - 1.0
        stored_vd_half_rad   = np.zeros( self.nsubs ) - 1.0
        stored_vd_fixed_rad1 = np.zeros( self.nsubs ) - 1.0
        stored_vd_fixed_rad5 = np.zeros( self.nsubs ) - 1.0
	stored_vd_outer	     = np.zeros( self.nsubs ) - 1.0

        for iii in np.arange(self.nsubs):
          if self.subcat['SubhaloLenType'][iii, 4] > 50: 	#self.min_n_stars:
            if iii % 1000 == 0:
                print "processing subhalo "+str(iii)+" out of "+str(self.nsubs)
            data =  illustris_python.snapshot.loadSubhalo(self.basePath, self.snapnum,iii,partType,fields=fields)
            for xyz in np.arange(3):
                data['Coordinates'][:,xyz] -= self.subcat['SubhaloPos'][iii,xyz]
                data['Coordinates'][ (data['Coordinates'][:,xyz] < -self.boxsize/2.0) ,xyz]  += self.boxsize
                data['Coordinates'][ (data['Coordinates'][:,xyz] >  self.boxsize/2.0) ,xyz]  -= self.boxsize
                data['Velocities'][:,xyz] -= np.median( data['Velocities'][:,xyz] )
         

            r = np.sqrt( data['Coordinates'][:,0]**2 + data['Coordinates'][:,1]**2 + data['Coordinates'][:,2]**2 )
	    vx = data['Velocities'][:,0]
            vy = data['Velocities'][:,1]
            vz = data['Velocities'][:,2]
	    age= np.array( data['GFM_StellarFormationTime'] )

	    index = (age > 0)
	    if np.sum(index) > 0:
                stored_vd_all[iii] = np.sqrt( np.std(vx[index])**2 + np.std(vy[index])**2 +np.std(vz[index])**2 )
    
            index = (r < self.subcat['SubhaloHalfmassRadType'][iii,partType]) & (age > 0)
            if np.sum(index) > 0:
                stored_vd_half_rad[iii] = np.sqrt( np.std(vx[index])**2 + np.std(vy[index])**2 +np.std(vz[index])**2 )

            index = (r > self.subcat['SubhaloHalfmassRadType'][iii,partType]) & (age > 0)
            if np.sum(index) > 0:
                stored_vd_outer[iii] = np.sqrt( np.std(vx[index])**2 + np.std(vy[index])**2 +np.std(vz[index])**2 )

            index = (r*scalefactor_list[self.snapnum] < 1.0) & (age > 0)
            if np.sum(index) > 0:
                stored_vd_fixed_rad1[iii] = np.sqrt( np.std(vx[index])**2 + np.std(vy[index])**2 +np.std(vz[index])**2 )

            index = (r*scalefactor_list[self.snapnum] < 5.0) & (age > 0)
            if np.sum(index) > 0:
                stored_vd_fixed_rad5[iii] = np.sqrt( np.std(vx[index])**2 + np.std(vy[index])**2 +np.std(vz[index])**2 )



        new_file_name = './catalogs/'+self.run+'/stellar_vel_disp/stellar_vel_disp_'+str(self.snapnum)+'.hdf5'
        print new_file_name
        f_write = hdf5lib.OpenFile(new_file_name, mode="w")
        subhalo_group = hdf5lib.CreateGroup(f_write, "Subhalo")
        hdf5lib.CreateArray(f_write, subhalo_group, 'StellarVelDisp_All'      , stored_vd_all)
        hdf5lib.CreateArray(f_write, subhalo_group, 'StellarVelDisp_HalfMassRad'  , stored_vd_half_rad)
        hdf5lib.CreateArray(f_write, subhalo_group, 'StellarVelDisp_1kpc_phys', stored_vd_fixed_rad1)
        hdf5lib.CreateArray(f_write, subhalo_group, 'StellarVelDisp_5kpc_phys', stored_vd_fixed_rad5)
	hdf5lib.CreateArray(f_write, subhalo_group, 'StellarVelDisp_outer', stored_vd_outer)        


        f_write.close()


    def make_stellar_core_mass_aux(self, **kwargs):
        fields=[u'Coordinates', u'Masses', u'GFM_StellarFormationTime']
        partType=4
        stored_core_1     = np.zeros( self.nsubs ) - 1.0
        stored_core_2     = np.zeros( self.nsubs ) - 1.0

        for iii in np.arange(self.nsubs):
          if self.subcat['SubhaloLenType'][iii, 4] > self.min_n_stars:
            if iii % 1000 == 0:
                print "processing subhalo "+str(iii)+" out of "+str(self.nsubs)
            data =  illustris_python.snapshot.loadSubhalo(self.basePath, self.snapnum,iii,partType,fields=fields)
            for xyz in np.arange(3):
                data['Coordinates'][:,xyz] -= self.subcat['SubhaloPos'][iii,xyz]
                data['Coordinates'][ (data['Coordinates'][:,xyz] < -self.boxsize/2.0) ,xyz]  += self.boxsize
                data['Coordinates'][ (data['Coordinates'][:,xyz] >  self.boxsize/2.0) ,xyz]  -= self.boxsize

            r  = np.sqrt( data['Coordinates'][:,0]**2 + data['Coordinates'][:,1]**2 + data['Coordinates'][:,2]**2 ) * self.scalefactor / self.little_h
            age= np.array( data['GFM_StellarFormationTime'] )

            index = (age > 0)
            if np.sum(index) > 0:
		core_1_index = r < 1.0
		if np.sum(core_1_index)>0:
		    stored_core_1[iii] = np.sum(data['Masses'][core_1_index]) * 1e10 / self.little_h
                core_2_index = r < 2.0
                if np.sum(core_2_index)>0:
                    stored_core_2[iii] = np.sum(data['Masses'][core_2_index]) * 1e10 / self.little_h

        new_file_name = './catalogs/'+self.run+'/stellar_core_mass/stellar_core_mass_'+str(self.snapnum)+'.hdf5'
        print new_file_name
        f_write = hdf5lib.OpenFile(new_file_name, mode="w")
        subhalo_group = hdf5lib.CreateGroup(f_write, "Subhalo")
        hdf5lib.CreateArray(f_write, subhalo_group, 'StellarCoreMass1'  , stored_core_1)
        hdf5lib.CreateArray(f_write, subhalo_group, 'StellarCoreMass2'	, stored_core_2)
        f_write.close()



    def make_ang_mom_aux(self, **kwargs):
        fields=[u'Coordinates', u'Velocities', u'Masses']
        partType=4
        stored_ang_mom_all        = np.zeros( (3,self.nsubs) ) 

        for iii in np.arange(self.nsubs):
          if self.subcat['SubhaloLenType'][iii, 4] > self.min_n_stars:
            if iii % 1000 == 0:
                print "processing subhalo "+str(iii)+" out of "+str(self.nsubs)
            data =  illustris_python.snapshot.loadSubhalo(self.basePath, self.snapnum,iii,partType,fields=fields)
            for xyz in np.arange(3):
                data['Coordinates'][:,xyz] -= self.subcat['SubhaloPos'][iii,xyz]
                data['Coordinates'][ (data['Coordinates'][:,xyz] < -self.boxsize/2.0) ,xyz]  += self.boxsize
                data['Coordinates'][ (data['Coordinates'][:,xyz] >  self.boxsize/2.0) ,xyz]  -= self.boxsize
                data['Velocities'][:,xyz] -= np.median( data['Velocities'][:,xyz] )


            lx  = np.sum( data['Masses'] * 1e10 / 0.704 * (data['Coordinates'][:,1] * data['Velocities'][:,2] - data['Coordinates'][:,2] * data['Velocities'][:,1] ) / 0.704 )
            ly  = np.sum( data['Masses'] * 1e10 / 0.704 * (data['Coordinates'][:,2] * data['Velocities'][:,0] - data['Coordinates'][:,0] * data['Velocities'][:,2] ) / 0.704 )
            lz  = np.sum( data['Masses'] * 1e10 / 0.704 * (data['Coordinates'][:,0] * data['Velocities'][:,1] - data['Coordinates'][:,1] * data['Velocities'][:,0] ) / 0.704 )

            stored_ang_mom_all[0,iii] = lx
            stored_ang_mom_all[1,iii] = ly
            stored_ang_mom_all[2,iii] = lz


        new_file_name = './catalogs/'+self.run+'/stellar_angular_momentum_vectors/stellar_angular_momentum_vectors_'+str(self.snapnum)+'.hdf5'
        print new_file_name
        f_write = hdf5lib.OpenFile(new_file_name, mode="w")
        subhalo_group = hdf5lib.CreateGroup(f_write, "Subhalo")
        hdf5lib.CreateArray(f_write, subhalo_group, 'StellarAngularMom'  , stored_ang_mom_all)
        f_write.close()








    def make_pca_axis_aux(self, **kwargs):
        fields=[u'Coordinates', u'Masses']
        eig_val_vec    = np.zeros( (6, 3, self.nsubs) ) 
        eig_val_vec[:] = np.nan
        pc_val_vec     = np.zeros( (6, 9, self.nsubs) )
        pc_val_vec[:]  = np.nan


        for partType in [4]:	#[0,1,4]:
            for iii in range(self.nsubs):	#[86206, 204399, 283832]:		#np.arange(1000)+1:		#self.nsubs):
              if self.subcat['SubhaloLenType'][iii, partType] > self.min_n_pca:
                if iii % 1000 == 0:
                    print "processing subhalo "+str(iii)+" out of "+str(self.nsubs)
                if (partType==0) or (partType==4):
                    fields=[u'Coordinates', u'Masses']
                    data =  illustris_python.snapshot.loadSubhalo(self.basePath, self.snapnum,iii,partType,fields=fields)
                else:
                    fields=[u'Coordinates', u'Velocities']
                    data =  illustris_python.snapshot.loadSubhalo(self.basePath, self.snapnum,iii,partType,fields=fields)


                for xyz in range(3):
                   data['Coordinates'][:,xyz] -= self.subcat['SubhaloPos'][iii,xyz]
                   data['Coordinates'][ (data['Coordinates'][:,xyz] < -self.boxsize/2.0) ,xyz]  += self.boxsize
                   data['Coordinates'][ (data['Coordinates'][:,xyz] >  self.boxsize/2.0) ,xyz]  -= self.boxsize

                print "data loaded, starting cov matrix calculation"

                sample_matrix = np.array( [ data['Coordinates'][:,xyz]  for xyz in range(3) ] )
                cov_matrix = np.cov( sample_matrix )
                #Calculating the eigenvalues and eigenvector for sample_matrix and the covariance matrix
                #The Eigenvector point in the direction of the principal component.

                print "cov matrix calculated, starting eigenvector calculation"

                eig_val, eig_vec = np.linalg.eig(cov_matrix)

                print "Eigenvectors calculated for subhalo {:d}".format( iii )

                p1 = []
                p2 = []
                p3 = []

                for i in range(len(eig_val)):
                    if eig_val[i] == np.amax(eig_val):
                        p1 = eig_vec[i]
                        e1 = eig_val[i]

                    if eig_val[i] != np.amax(eig_val) and eig_val[i] != np.amin(eig_val):
                        p2 = eig_vec[i]
                        e2 = eig_val[i]

                    if eig_val[i] == np.amin(eig_val):
                        p3 = eig_vec[i]
                        e3 = eig_val[i]

                #Storing the principal components
                pc = np.concatenate( (p1[:], p2[:], p3[:]) )
                EIG_VAL = [e1, e2, e3]

                eig_val_vec[partType,:,iii] = EIG_VAL
                pc_val_vec[ partType,:,iii]  = pc
                
        new_file_name = './catalogs/'+self.run+'/pca_axes/pca_axes_'+str(self.snapnum)+'.hdf5'
        print new_file_name
        f_write = hdf5lib.OpenFile(new_file_name, mode="w")
        subhalo_group = hdf5lib.CreateGroup(f_write, "Subhalo")
        hdf5lib.CreateArray(f_write, subhalo_group, 'eigenvalues'  , eig_val_vec)
        hdf5lib.CreateArray(f_write, subhalo_group, 'pc_val_vec'  , pc_val_vec)
        f_write.close()




    def make_gas_metallicity_gradient(self, **kwargs):
        fields=[u'Coordinates', u'GFM_Metallicity', u'StarFormationRate']
        partType=0
        wdata = {}
        for name in ['z0_5', 'z0_10', 'z0_20', 'z0_rh', 'dz_5', 'dz_10', 'dz_20', 'dz_rh']:
            wdata[name] = np.zeros( self.nsubs )
            wdata[name][:] = np.nan

        for name in ['z0_5_ism', 'z0_10_ism', 'z0_20_ism', 'z0_rh_ism', 'dz_5_ism', 'dz_10_ism', 'dz_20_ism', 'dz_rh_ism']:
            wdata[name] = np.zeros( self.nsubs )
            wdata[name][:] = np.nan

        for name in ['nuclear_metallicity_1kpc', 'nuclear_metallicity_5kpc', 'nuclear_metallicity_1kpc_ism', 'nuclear_metallicity_5kpc_ism']:
            wdata[name] = np.zeros( self.nsubs )
            wdata[name][:] = np.nan

#        z0_5  = np.zeros( self.nsubs )   - 1.0
#        z0_10 = np.zeros( self.nsubs )  - 1.0
#        z0_20 = np.zeros( self.nsubs )  - 1.0
#        dz_5  = np.zeros( self.nsubs )  - 99.0
#        dz_10 = np.zeros( self.nsubs ) - 99.0
#        dz_20 = np.zeros( self.nsubs ) - 99.0
      
        for iii in np.arange(self.nsubs):
            if self.subcat['SubhaloLenType'][iii, 0] > self.min_n_gas:
                if iii % 1000 == 0:
                    print "processing subhalo "+str(iii)+" out of "+str(self.nsubs)
                data =  illustris_python.snapshot.loadSubhalo(self.basePath, self.snapnum,iii,partType,fields=fields)
                for xyz in np.arange(3):
                    data['Coordinates'][:,xyz] -= self.subcat['SubhaloPos'][iii,xyz]
                    data['Coordinates'][ (data['Coordinates'][:,xyz] < -self.boxsize/2.0) ,xyz]  += self.boxsize
                    data['Coordinates'][ (data['Coordinates'][:,xyz] >  self.boxsize/2.0) ,xyz]  -= self.boxsize
        
                r = np.sqrt( data['Coordinates'][:,0]**2 + data['Coordinates'][:,1]**2 + data['Coordinates'][:,2]**2 )
                met = np.log10(data['GFM_Metallicity'])
                sfr = data['StarFormationRate']

                index = (r<5)
                if (np.sum(index) > 30 ): result_5 = np.polyfit(r[r<5],  met[r<5],  1)
                else: result_5 = np.array([0,0])
                index = (r<10)
                if (np.sum(index) > 30 ): result_10= np.polyfit(r[r<10], met[r<10], 1)
                else: result_10 = np.array([0,0])
                index = (r<20)
                if (np.sum(index) > 30): result_20= np.polyfit(r[r<20], met[r<20], 1)
                else: result_20 = np.array([0,0])
                index = (r<self.subcat['SubhaloHalfmassRadType'][iii,4])
                if (np.sum(index) > 30): result_rh= np.polyfit(r[index], met[index], 1)
                else: result_rh = np.array([0,0])
                

                wdata['z0_5'][iii]  = result_5[1]
                wdata['z0_10'][iii] = result_10[1]
                wdata['z0_20'][iii] = result_20[1]
                wdata['z0_rh'][iii] = result_rh[1]

                wdata['dz_5'][iii]  = result_5[0]
                wdata['dz_10'][iii] = result_10[0]
                wdata['dz_20'][iii] = result_20[0]
                wdata['dz_rh'][iii] = result_rh[0]

                index = (r<1)
                if np.sum(index) > 10:
                    wdata['nuclear_metallicity_1kpc'] = np.median( met[index] )
                index = (r<5)
                if np.sum(index) > 10:
                    wdata['nuclear_metallicity_5kpc'] = np.median( met[index] )


                r = r[ sfr > 0 ]
                met = met[ sfr > 0]
                sfr = sfr[ sfr > 0]
                index = (r<5)
                if (np.sum(index) > 30 ): result_5 = np.polyfit(r[r<5],  met[r<5],  1)
                else: result_5 = np.array([0,0])
                index = (r<10)
                if (np.sum(index) > 30 ): result_10= np.polyfit(r[r<10], met[r<10], 1)
                else: result_10 = np.array([0,0])
                index = (r<20)
                if (np.sum(index) > 30): result_20= np.polyfit(r[r<20], met[r<20], 1)
                else: result_20 = np.array([0,0])
                index = (r<self.subcat['SubhaloHalfmassRadType'][iii,4])
                if (np.sum(index) > 30): result_rh= np.polyfit(r[index], met[index], 1)
                else: result_rh = np.array([0,0])

                wdata['z0_5_ism'][iii]  = result_5[1]
                wdata['z0_10_ism'][iii] = result_10[1]
                wdata['z0_20_ism'][iii] = result_20[1]
                wdata['z0_rh_ism'][iii] = result_rh[1]

                wdata['dz_5_ism'][iii]  = result_5[0]
                wdata['dz_10_ism'][iii] = result_10[0]
                wdata['dz_20_ism'][iii] = result_20[0]
                wdata['dz_rh_ism'][iii] = result_rh[0]

                index = (r<1)
                if np.sum(index) > 10:
                    wdata['nuclear_metallicity_1kpc_ism'] = np.median( met[index] )
                index = (r<5)
                if np.sum(index) > 10:
                    wdata['nuclear_metallicity_5kpc_ism'] = np.median( met[index] )


        new_file_name = './catalogs/'+self.run+'/gas_metallicity/gas_metallicity_info_'+str(self.snapnum)+'.hdf5'

        print new_file_name
        f = hdf5lib.OpenFile(new_file_name, mode="w")
        subhalo_group = hdf5lib.CreateGroup(f, "Subhalo")
        hdf5lib.CreateArray(f, subhalo_group, 'CentralMetallicity_5',  wdata['z0_5'])
        hdf5lib.CreateArray(f, subhalo_group, 'CentralMetallicity_10', wdata['z0_10'])
        hdf5lib.CreateArray(f, subhalo_group, 'CentralMetallicity_20', wdata['z0_20'])
        hdf5lib.CreateArray(f, subhalo_group, 'CentralMetallicity_rh', wdata['z0_rh'])

        hdf5lib.CreateArray(f, subhalo_group, 'GradMetallicity_5',  wdata['dz_5'])
        hdf5lib.CreateArray(f, subhalo_group, 'GradMetallicity_10', wdata['dz_10'])
        hdf5lib.CreateArray(f, subhalo_group, 'GradMetallicity_20', wdata['dz_20'])
        hdf5lib.CreateArray(f, subhalo_group, 'GradMetallicity_rh', wdata['dz_rh'])
        

        hdf5lib.CreateArray(f, subhalo_group, 'CentralMetallicity_5_ISM',  wdata['z0_5_ism'])
        hdf5lib.CreateArray(f, subhalo_group, 'CentralMetallicity_10_ISM', wdata['z0_10_ism'])
        hdf5lib.CreateArray(f, subhalo_group, 'CentralMetallicity_20_ISM', wdata['z0_20_ism'])
        hdf5lib.CreateArray(f, subhalo_group, 'CentralMetallicity_rh_ISM', wdata['z0_rh_ism'])

        hdf5lib.CreateArray(f, subhalo_group, 'GradMetallicity_5_ISM',  wdata['dz_5_ism'])
        hdf5lib.CreateArray(f, subhalo_group, 'GradMetallicity_10_ISM', wdata['dz_10_ism'])
        hdf5lib.CreateArray(f, subhalo_group, 'GradMetallicity_20_ISM', wdata['dz_20_ism'])
        hdf5lib.CreateArray(f, subhalo_group, 'GradMetallicity_rh_ISM', wdata['dz_rh_ism'])

        hdf5lib.CreateArray(f, subhalo_group, 'NuclearMetallicity_1kpc', wdata['nuclear_metallicity_1kpc'])
        hdf5lib.CreateArray(f, subhalo_group, 'NuclearMetallicity_5kpc', wdata['nuclear_metallicity_5kpc'])
        hdf5lib.CreateArray(f, subhalo_group, 'NuclearMetallicity_1kpc_ISM', wdata['nuclear_metallicity_1kpc_ism'])
        hdf5lib.CreateArray(f, subhalo_group, 'NuclearMetallicity_5kpc_ISM', wdata['nuclear_metallicity_5kpc_ism'])

#nuclear_metallicity_5kpc_ism

        f.close()



    def make_ism_mass_and_size(self, h=0.6774, **kwargs):
        fields=[u'Coordinates', u'Masses', u'StarFormationRate', ]
        partType=0
        ism_mass = np.zeros( self.nsubs ) 
        ism_hmr  = np.zeros( self.nsubs )
        ism_maxr = np.zeros( self.nsubs )

        ism_dyn_time = np.zeros( self.nsubs )           # dynamical time at hmr
        ism_dyn_time2 = np.zeros( self.nsubs )          # dynamical time at maxr
        ism_dyn_time3 = np.zeros( self.nsubs )          # dynamical time at twice maxr
        ism_dyn_time4 = np.zeros( self.nsubs )          # dynamical time at twice maxr
        ism_dyn_time5 = np.zeros( self.nsubs )          # dynamical time at twice maxr
        ism_dyn_time6 = np.zeros( self.nsubs )          # dynamical time at twice maxr
        ism_dyn_time7 = np.zeros( self.nsubs )          # dynamical time at twice maxr
        ism_dyn_time8 = np.zeros( self.nsubs )          # dynamical time at twice maxr
        ism_dyn_time9 = np.zeros( self.nsubs )          # dynamical time at twice maxr

        ism_mass[:] = np.NAN
        ism_hmr[:]  = np.NAN
        ism_maxr[:] = np.NAN

        ism_dyn_time[:] = np.NAN
        ism_dyn_time2[:] = np.NAN
        ism_dyn_time3[:] = np.NAN
        ism_dyn_time4[:] = np.NAN
        ism_dyn_time5[:] = np.NAN
        ism_dyn_time6[:] = np.NAN
        ism_dyn_time7[:] = np.NAN
        ism_dyn_time8[:] = np.NAN
        ism_dyn_time9[:] = np.NAN


        if self.nsubs > 0:
            mask = (self.subcat['SubhaloLenType'][:, 0] > 10) & (self.subcat['SubhaloSFR'][:] > 0) 
            subnr= np.arange( len(mask) )
            mask = (self.subcat['SubhaloLenType'][:, 0] > 10) & (self.subcat['SubhaloSFR'][:] > 0)  & (self.subcat['SubhaloMassType'][:,4] > 0) & (self.subcat['SubhaloLenType'][:,4] > 10) 
            # & (subnr > 500000) 
#            print len(mask), np.sum(mask)
            #sys.exit()
            tmp = self.subcat['SubhaloMassType'][mask,4] * 1e10 / self.little_h
#            print np.min(tmp), np.min(tmp[tmp>0]), np.max(tmp)

            iii_list = subnr[ mask ]
            # map(myfunc, iii_list)
            for tmp_iii,iii in enumerate(iii_list):            #np.arange(self.nsubs):
                if True:            #self.subcat['SubhaloLenType'][iii, 0] > self.min_n_gas:
                    if tmp_iii % 100 == 0:
                        print "processing subhalo "+str(tmp_iii)+" out of "+str(len(iii_list))
                        print "    loading data..."

                    if self.subbox != None:
                        #full_path = self.base+'/output_subbox'+str(self.subbox)
                        data={}
                        data['Coordinates']       = self.halo_reader.read('POS ', partType, -1, iii)
                        data['Masses']            = self.halo_reader.read('MASS', partType, -1, iii)
                        data['StarFormationRate'] = self.halo_reader.read('SFR ', partType, -1, iii)

                        data['AllCoords'] = np.array([])
                        data['AllMasses'] = np.array([])

                        for tmp_partType in range(6):
                            data['AllCoords'] = np.append( data['AllCoords'],  self.halo_reader.read('POS ', tmp_partType, -1, iii)   )
                            data['AllMasses'] = np.append( data['AllMasses'],  self.halo_reader.read('MASS', tmp_partType, -1, iii)   )

                    else:
                        data =  illustris_python.snapshot.loadSubhalo(self.basePath, self.snapnum,iii,partType,fields=fields)

                        tmp_all_x = np.array([])
                        tmp_all_y = np.array([])
                        tmp_all_z = np.array([])
                        data['AllMasses'] = np.array([])
                        for tmp_partType in range(6):
                            tmp_data = illustris_python.snapshot.loadSubhalo(self.basePath, self.snapnum,iii,partType,fields=[r'Coordinates', r'Masses'])
                            tmp_all_x = np.append( tmp_all_x ,  tmp_data['Coordinates'][:,0]   )
                            tmp_all_y = np.append( tmp_all_y ,  tmp_data['Coordinates'][:,1]   )
                            tmp_all_z = np.append( tmp_all_z ,  tmp_data['Coordinates'][:,2]   )
                            data['AllMasses'] = np.append( data['AllMasses'], tmp_data['Masses'] )
                        data['AllCoords'] = np.column_stack( (tmp_all_x, tmp_all_y, tmp_all_z)  )

                    if tmp_iii % 100 == 0:
                        print "    extracting coordinates and calculating radii..."
                    for xyz in np.arange(3):
                        data['Coordinates'][:,xyz] -= self.subcat['SubhaloPos'][iii,xyz]
                        data['Coordinates'][ (data['Coordinates'][:,xyz] < -self.boxsize/2.0) ,xyz]  += self.boxsize
                        data['Coordinates'][ (data['Coordinates'][:,xyz] >  self.boxsize/2.0) ,xyz]  -= self.boxsize

                        data['AllCoords'][:,xyz] -= self.subcat['SubhaloPos'][iii,xyz]
                        data['AllCoords'][ (data['AllCoords'][:,xyz] < -self.boxsize/2.0) ,xyz]  += self.boxsize
                        data['AllCoords'][ (data['AllCoords'][:,xyz] >  self.boxsize/2.0) ,xyz]  -= self.boxsize


                    r = np.sqrt( data['Coordinates'][:,0]**2 + data['Coordinates'][:,1]**2 + data['Coordinates'][:,2]**2 )
                    m = data['Masses'] * 1e10 / self.little_h
                    sfr = data['StarFormationRate']


                    rall = np.sqrt( data['AllCoords'][:,0]**2 + data['AllCoords'][:,1]**2 + data['AllCoords'][:,2]**2 )
                    mall = data['AllMasses'] * 1e10 / self.little_h


                    index = sfr > 0
                    if np.sum(index) > 5:
                        r = r[index]
                        m = m[index]
                        sfr = sfr[index]

                        ism_mass[iii] = np.sum( m )

                        tmp = self.subcat['SubhaloMassType'][iii,4] * 1e10 / self.little_h
#                        print "Summary:  Galaxy {:d} has a (log) stellar mass of {:.2f} and a (log) ISM mass of {:.2f}; gf={:.2f}".format( iii, np.log10(tmp), np.log10( ism_mass[iii]  ) , ism_mass[iii] / tmp    )

                        if tmp_iii % 100 == 0:
                            print "    sorting radii..."
                        order = np.argsort( r ) 
                        r = r[order]
                        m = np.cumsum( m[order] ) / np.sum(m)

                        m = np.abs( m - 0.5 )
                        index_of_hmr = np.where( m == np.min( m ) )[0][0]
                        ism_hmr[iii] = r[index_of_hmr]
                        ism_maxr[iii] = r[-1]

                        grav_g = 6.67e-11 * (1.99e30 * 3.24e-20 * 3.24e-20) * 3.24e-20 * 3.15e16 * 3.15e16

#                        order = np.argsort( rall )
#                        rall = rall[order]

#                        tmp_diff = np.abs( rall -  ism_hmr[iii] )
                        do_dyn_time_calcs = False
                        if do_dyn_time_calcs:
                            this_r = ism_hmr[iii]
                            dyn_time_index  = rall < this_r              #anp.where(  tmp_diff == np.min(tmp_diff)   )[0][0]
                            ism_dyn_time[iii] = np.sqrt(  (self.scalefactor * this_r)**3  / (grav_g * np.sum(mall[dyn_time_index])  )   )


                            #tmp_diff = np.abs( rall -  ism_maxr[iii] )

                            this_r = ism_maxr[iii]
                            dyn_time_index  = rall < this_r         #np.where(  tmp_diff == np.min(tmp_diff)   )[0][0]
                            ism_dyn_time2[iii] = np.sqrt(  (self.scalefactor * this_r)**3 / (grav_g * np.sum(mall[dyn_time_index]) )   )



#                            tmp_diff = np.abs( rall -  2.0*ism_maxr[iii] )
                            this_r = 2.0*ism_maxr[iii]
                            dyn_time_index  = rall < this_r         #np.where(  tmp_diff == np.min(tmp_diff)   )[0][0]
                            ism_dyn_time3[iii] = np.sqrt(  (self.scalefactor * this_r)**3 / (grav_g * np.sum(mall[dyn_time_index]) )   ) 


#                            tmp_diff = np.abs( rall -  5.0*ism_maxr[iii] )
                            this_r = 5.0*ism_maxr[iii]
                            dyn_time_index  = rall < this_r         #np.where(  tmp_diff == np.min(tmp_diff)   )[0][0]
                            ism_dyn_time4[iii] = np.sqrt(  (self.scalefactor * this_r)**3 / (grav_g * np.sum(mall[dyn_time_index]) )   ) 



#                            tmp_diff = np.abs( rall -  self.subcat['SubhaloHalfmassRadType'][iii,4] )
#                            dyn_time_index  = np.where(  tmp_diff == np.min(tmp_diff)   )[0][0]

                            this_rad = self.subcat['SubhaloHalfmassRadType'][iii,4]
                            dyn_time_indicies = rall < this_rad
                            ism_dyn_time5[iii] = np.sqrt(  (self.scalefactor * this_rad )**3 / (grav_g * np.sum(mall[dyn_time_indicies])  )   ) 
#                            ism_dyn_time5[iii] = np.sqrt(  (self.scalefactor * rall[dyn_time_index])**3 / (grav_g * mall[dyn_time_index] * 1e10 / self.little_h )   ) 

#                            tmp_diff = np.abs( rall -  2.0*self.subcat['SubhaloHalfmassRadType'][iii,4] )
                            this_r = 2.0*self.subcat['SubhaloHalfmassRadType'][iii,4]
                            dyn_time_index  = rall < this_r         #np.where(  tmp_diff == np.min(tmp_diff)   )[0][0]
                            ism_dyn_time6[iii] = np.sqrt(  (self.scalefactor * this_r)**3 / (grav_g * np.sum(mall[dyn_time_index]) )   ) 

#                            tmp_diff = np.abs( rall -  5.0*self.subcat['SubhaloHalfmassRadType'][iii,4] )
                            this_r = 5.0*self.subcat['SubhaloHalfmassRadType'][iii,4]
                            dyn_time_index  = rall < this_r #np.where(  tmp_diff == np.min(tmp_diff)   )[0][0]
                            ism_dyn_time7[iii] = np.sqrt(  (self.scalefactor * this_r)**3 / (grav_g * np.sum(mall[dyn_time_index])  )  ) 

                            tmp_diff = np.abs( rall -  10.0*self.subcat['SubhaloHalfmassRadType'][iii,4] )
                            this_r = 10.0*self.subcat['SubhaloHalfmassRadType'][iii,4]
                            dyn_time_index  = rall < this_r #np.where(  tmp_diff == np.min(tmp_diff)   )[0][0]
                            ism_dyn_time8[iii] = np.sqrt(  (self.scalefactor * this_r)**3 / (grav_g * np.sum(mall[dyn_time_index]) )   ) 

                            tmp_diff = np.abs( rall -  self.subcat['SubhaloHalfmassRad'][iii] )
                            this_r = self.subcat['SubhaloHalfmassRad'][iii]
                            dyn_time_index  = rall < this_r # np.where(  tmp_diff == np.min(tmp_diff)   )[0][0]
                            ism_dyn_time9[iii] = np.sqrt(  (self.scalefactor * this_r)**3 / (grav_g * np.sum(mall[dyn_time_index]) )   ) 
                        #print iii, this_r , np.sum( mall[dyn_time_index] )
                        #print np.sqrt(  (self.scalefactor * this_r)**3 / (grav_g * np.sum(mall[dyn_time_index]) )   )
                        #print ism_dyn_time
                        #print ism_dyn_time9



#self.subcat['SubhaloPos'][iii,xyz]



                        if tmp_iii % 100 == 0:
                            print "    done..."
                    else:
                        ism_mass[iii] = np.NAN
                        ism_hmr[iii] = np.NAN
                        ism_maxr[iii] = np.NAN
                        ism_dyn_time[iii] = np.NAN
                        ism_dyn_time2[iii] = np.NAN
                        ism_dyn_time3[iii] = np.NAN
                        ism_dyn_time4[iii] = np.NAN
                        ism_dyn_time5[iii] = np.NAN
                        ism_dyn_time6[iii] = np.NAN
                        ism_dyn_time7[iii] = np.NAN
                        ism_dyn_time8[iii] = np.NAN
                        ism_dyn_time9[iii] = np.NAN

        new_file_name = './catalogs/'+self.run+'/ism_mass_and_size/ism_mass_and_size_info_'+str(self.snapnum)+'.hdf5'

        try:
            print self.subcat['SubhaloMassType'][:, 4] 
            print np.min(ism_mass[ np.isfinite(ism_mass) ]), np.max( ism_mass[ np.isfinite(ism_mass) ] )
        except:
            print "failed at some print statements."

        print new_file_name
        f = hdf5lib.OpenFile(new_file_name, mode="w")
        subhalo_group = hdf5lib.CreateGroup(f, "Subhalo")
        hdf5lib.CreateArray(f, subhalo_group, 'ISM_mass',  ism_mass)
        hdf5lib.CreateArray(f, subhalo_group, 'ISM_HalfMassRadius', ism_hmr)
        hdf5lib.CreateArray(f, subhalo_group, 'ISM_MaxRadius', ism_maxr)


        hdf5lib.CreateArray(f, subhalo_group, 'ISM_HalfMassRadius_DynTime', ism_dyn_time)
        hdf5lib.CreateArray(f, subhalo_group, 'ISM_MaxRadius_DynTime',      ism_dyn_time2)
        hdf5lib.CreateArray(f, subhalo_group, 'ISM_TwiceMaxRadius_DynTime', ism_dyn_time3)
        hdf5lib.CreateArray(f, subhalo_group, 'ISM_FiveMaxRadius_DynTime', ism_dyn_time4)

        hdf5lib.CreateArray(f, subhalo_group, 'ISM_StellarHalfMassRadius_DynTime', ism_dyn_time5)
        hdf5lib.CreateArray(f, subhalo_group, 'ISM_TwiceStellarHalfMassRadius_DynTime', ism_dyn_time6)
        hdf5lib.CreateArray(f, subhalo_group, 'ISM_FiveStellarHalfMassRadius_DynTime', ism_dyn_time7)
        hdf5lib.CreateArray(f, subhalo_group, 'ISM_TenStellarHalfMassRadius_DynTime', ism_dyn_time8)

        hdf5lib.CreateArray(f, subhalo_group, 'ISM_TotalHalfMassRadius_DynTime', ism_dyn_time9)

        f.close()






    def make_timescales(self, h=0.6774, **kwargs):
        fields=[u'Coordinates', u'Masses', u'StarFormationRate', ]
        partType=0

        ism_dyn_time = np.zeros( self.nsubs )           # dynamical time at hmr
        ism_dyn_time2 = np.zeros( self.nsubs )          # dynamical time at maxr
        ism_dyn_time3 = np.zeros( self.nsubs )          # dynamical time at twice maxr
        ism_dyn_time4 = np.zeros( self.nsubs )          # dynamical time at twice maxr
        ism_dyn_time5 = np.zeros( self.nsubs )          # dynamical time at twice maxr
        ism_dyn_time6 = np.zeros( self.nsubs )          # dynamical time at twice maxr
        ism_dyn_time7 = np.zeros( self.nsubs )          # dynamical time at twice maxr
        ism_dyn_time8 = np.zeros( self.nsubs )          # dynamical time at twice maxr
        ism_dyn_time9 = np.zeros( self.nsubs )          # dynamical time at twice maxr

        ism_dyn_time[:] = np.NAN
        ism_dyn_time2[:] = np.NAN
        ism_dyn_time3[:] = np.NAN
        ism_dyn_time4[:] = np.NAN
        ism_dyn_time5[:] = np.NAN
        ism_dyn_time6[:] = np.NAN
        ism_dyn_time7[:] = np.NAN
        ism_dyn_time8[:] = np.NAN
        ism_dyn_time9[:] = np.NAN


        if self.nsubs > 0:
            subnr= np.arange( len(self.subcat['SubhaloSFR'][:]) )
            mask = (self.subcat['SubhaloLenType'][:, 0] > 10) & (self.subcat['SubhaloSFR'][:] > 0)  & (self.subcat['SubhaloMassType'][:,4] > 0) & (self.subcat['SubhaloLenType'][:,4] > 10) 
            tmp = self.subcat['SubhaloMassType'][mask,4] * 1e10 / self.little_h

            iii_list = subnr[ mask ]
            # map(myfunc, iii_list)
            for tmp_iii,iii in enumerate(iii_list):            #np.arange(self.nsubs):
                if True:            #self.subcat['SubhaloLenType'][iii, 0] > self.min_n_gas:
                    if tmp_iii % 100 == 0:
                        print "processing subhalo "+str(tmp_iii)+" out of "+str(len(iii_list))
                        print "    loading data..."

                    if self.subbox != None:
                        #full_path = self.base+'/output_subbox'+str(self.subbox)
                        data={}
                        data['Coordinates']       = self.halo_reader.read('POS ', partType, -1, iii)
                        data['Masses']            = self.halo_reader.read('MASS', partType, -1, iii)
                        data['StarFormationRate'] = self.halo_reader.read('SFR ', partType, -1, iii)

                        data['AllCoords'] = np.array([])
                        data['AllMasses'] = np.array([])

                        for tmp_partType in range(6):
                            data['AllCoords'] = np.append( data['AllCoords'],  self.halo_reader.read('POS ', tmp_partType, -1, iii)   )
                            data['AllMasses'] = np.append( data['AllMasses'],  self.halo_reader.read('MASS', tmp_partType, -1, iii)   )

                    else:
                        data =  illustris_python.snapshot.loadSubhalo(self.basePath, self.snapnum,iii,partType,fields=fields)

                        tmp_all_x = np.array([])
                        tmp_all_y = np.array([])
                        tmp_all_z = np.array([])
                        data['AllMasses'] = np.array([])
                        for tmp_partType in range(6):
                            tmp_data = illustris_python.snapshot.loadSubhalo(self.basePath, self.snapnum,iii,partType,fields=[r'Coordinates', r'Masses'])
                            tmp_all_x = np.append( tmp_all_x ,  tmp_data['Coordinates'][:,0]   )
                            tmp_all_y = np.append( tmp_all_y ,  tmp_data['Coordinates'][:,1]   )
                            tmp_all_z = np.append( tmp_all_z ,  tmp_data['Coordinates'][:,2]   )
                            data['AllMasses'] = np.append( data['AllMasses'], tmp_data['Masses'] )
                        data['AllCoords'] = np.column_stack( (tmp_all_x, tmp_all_y, tmp_all_z)  )

                    if tmp_iii % 100 == 0:
                        print "    extracting coordinates and calculating radii..."
                    for xyz in np.arange(3):
                        data['Coordinates'][:,xyz] -= self.subcat['SubhaloPos'][iii,xyz]
                        data['Coordinates'][ (data['Coordinates'][:,xyz] < -self.boxsize/2.0) ,xyz]  += self.boxsize
                        data['Coordinates'][ (data['Coordinates'][:,xyz] >  self.boxsize/2.0) ,xyz]  -= self.boxsize

                        data['AllCoords'][:,xyz] -= self.subcat['SubhaloPos'][iii,xyz]
                        data['AllCoords'][ (data['AllCoords'][:,xyz] < -self.boxsize/2.0) ,xyz]  += self.boxsize
                        data['AllCoords'][ (data['AllCoords'][:,xyz] >  self.boxsize/2.0) ,xyz]  -= self.boxsize


                    r = np.sqrt( data['Coordinates'][:,0]**2 + data['Coordinates'][:,1]**2 + data['Coordinates'][:,2]**2 )
                    m = data['Masses'] * 1e10 / self.little_h
                    sfr = data['StarFormationRate']


                    rall = np.sqrt( data['AllCoords'][:,0]**2 + data['AllCoords'][:,1]**2 + data['AllCoords'][:,2]**2 )
                    mall = data['AllMasses'] * 1e10 / self.little_h


                    index = sfr > 0
                    if np.sum(index) > 5:
                        r = r[index]
                        m = m[index]
                        sfr = sfr[index]

                        tmp = self.subcat['SubhaloMassType'][iii,4] * 1e10 / self.little_h

                        if tmp_iii % 100 == 0:
                            print "    sorting radii..."
                        order = np.argsort( r ) 
                        r = r[order]
                        m = np.cumsum( m[order] ) / np.sum(m)

                        m = np.abs( m - 0.5 )
                        index_of_hmr = np.where( m == np.min( m ) )[0][0]
                        hmr = r[index_of_hmr]

                        grav_g = 6.67e-11 * (1.99e30 * 3.24e-20 * 3.24e-20) * 3.24e-20 * 3.15e16 * 3.15e16

                        do_dyn_time_calcs = True
                        if do_dyn_time_calcs:
                            this_r = hmr
                            dyn_time_index  = rall < this_r              #anp.where(  tmp_diff == np.min(tmp_diff)   )[0][0]
                            ism_dyn_time[iii] = np.sqrt(  (self.scalefactor * this_r)**3  / (grav_g * np.sum(mall[dyn_time_index])  )   )

                            this_r = np.max(r)
                            dyn_time_index  = rall < this_r         #np.where(  tmp_diff == np.min(tmp_diff)   )[0][0]
                            ism_dyn_time2[iii] = np.sqrt(  (self.scalefactor * this_r)**3 / (grav_g * np.sum(mall[dyn_time_index]) )   )

                            this_r = 2.0*np.max(r)      #ism_maxr[iii]
                            dyn_time_index  = rall < this_r         #np.where(  tmp_diff == np.min(tmp_diff)   )[0][0]
                            ism_dyn_time3[iii] = np.sqrt(  (self.scalefactor * this_r)**3 / (grav_g * np.sum(mall[dyn_time_index]) )   ) 

                            this_r = 5.0*np.max(r)      #ism_maxr[iii]
                            dyn_time_index  = rall < this_r         #np.where(  tmp_diff == np.min(tmp_diff)   )[0][0]
                            ism_dyn_time4[iii] = np.sqrt(  (self.scalefactor * this_r)**3 / (grav_g * np.sum(mall[dyn_time_index]) )   ) 

                            this_rad = self.subcat['SubhaloHalfmassRadType'][iii,4]
                            dyn_time_indicies = rall < this_rad
                            ism_dyn_time5[iii] = np.sqrt(  (self.scalefactor * this_rad )**3 / (grav_g * np.sum(mall[dyn_time_indicies])  )   ) 

                            this_r = 2.0*self.subcat['SubhaloHalfmassRadType'][iii,4]
                            dyn_time_index  = rall < this_r         #np.where(  tmp_diff == np.min(tmp_diff)   )[0][0]
                            ism_dyn_time6[iii] = np.sqrt(  (self.scalefactor * this_r)**3 / (grav_g * np.sum(mall[dyn_time_index]) )   ) 

                            this_r = 5.0*self.subcat['SubhaloHalfmassRadType'][iii,4]
                            dyn_time_index  = rall < this_r #np.where(  tmp_diff == np.min(tmp_diff)   )[0][0]
                            ism_dyn_time7[iii] = np.sqrt(  (self.scalefactor * this_r)**3 / (grav_g * np.sum(mall[dyn_time_index])  )  ) 

                            tmp_diff = np.abs( rall -  10.0*self.subcat['SubhaloHalfmassRadType'][iii,4] )
                            this_r = 10.0*self.subcat['SubhaloHalfmassRadType'][iii,4]
                            dyn_time_index  = rall < this_r #np.where(  tmp_diff == np.min(tmp_diff)   )[0][0]
                            ism_dyn_time8[iii] = np.sqrt(  (self.scalefactor * this_r)**3 / (grav_g * np.sum(mall[dyn_time_index]) )   ) 

                            tmp_diff = np.abs( rall -  self.subcat['SubhaloHalfmassRad'][iii] )
                            this_r = self.subcat['SubhaloHalfmassRad'][iii]
                            dyn_time_index  = rall < this_r # np.where(  tmp_diff == np.min(tmp_diff)   )[0][0]
                            ism_dyn_time9[iii] = np.sqrt(  (self.scalefactor * this_r)**3 / (grav_g * np.sum(mall[dyn_time_index]) )   ) 

                        if tmp_iii % 100 == 0:
                            print "    done..."
                    else:
                        ism_dyn_time[iii] = np.NAN
                        ism_dyn_time2[iii] = np.NAN
                        ism_dyn_time3[iii] = np.NAN
                        ism_dyn_time4[iii] = np.NAN
                        ism_dyn_time5[iii] = np.NAN
                        ism_dyn_time6[iii] = np.NAN
                        ism_dyn_time7[iii] = np.NAN
                        ism_dyn_time8[iii] = np.NAN
                        ism_dyn_time9[iii] = np.NAN

        new_file_name = './catalogs/'+self.run+'/timescales/timescales_'+str(self.snapnum)+'.hdf5'

        try:
            print self.subcat['SubhaloMassType'][:, 4] 
            print np.min(ism_mass[ np.isfinite(ism_mass) ]), np.max( ism_mass[ np.isfinite(ism_mass) ] )
        except:
            print "failed at some print statements."

        f = hdf5lib.OpenFile(new_file_name, mode="w")
        subhalo_group = hdf5lib.CreateGroup(f, "Subhalo")
        hdf5lib.CreateArray(f, subhalo_group, 'ISM_HalfMassRadius_DynTime', ism_dyn_time)
        hdf5lib.CreateArray(f, subhalo_group, 'ISM_MaxRadius_DynTime',      ism_dyn_time2)
        hdf5lib.CreateArray(f, subhalo_group, 'ISM_TwiceMaxRadius_DynTime', ism_dyn_time3)
        hdf5lib.CreateArray(f, subhalo_group, 'ISM_FiveMaxRadius_DynTime', ism_dyn_time4)

        hdf5lib.CreateArray(f, subhalo_group, 'ISM_StellarHalfMassRadius_DynTime', ism_dyn_time5)
        hdf5lib.CreateArray(f, subhalo_group, 'ISM_TwiceStellarHalfMassRadius_DynTime', ism_dyn_time6)
        hdf5lib.CreateArray(f, subhalo_group, 'ISM_FiveStellarHalfMassRadius_DynTime', ism_dyn_time7)
        hdf5lib.CreateArray(f, subhalo_group, 'ISM_TenStellarHalfMassRadius_DynTime', ism_dyn_time8)

        hdf5lib.CreateArray(f, subhalo_group, 'ISM_TotalHalfMassRadius_DynTime', ism_dyn_time9)

        f.close()





    def make_gas_disk_scale_length(self, **kwargs):
        fields=[u'Coordinates', u'Masses', u'Density']
        partType=0
        sigma0_5 = np.zeros( self.nsubs )   - 1.0	# central surface density
        sigma0_10 = np.zeros( self.nsubs )  - 1.0
        sigma0_20 = np.zeros( self.nsubs )  - 1.0
        sl_5 = np.zeros( self.nsubs )  - 99.0		# disk scale length (ln -- exponential) 
        sl_10 = np.zeros( self.nsubs ) - 99.0
        sl_20 = np.zeros( self.nsubs ) - 99.0

        for iii in np.arange(self.nsubs):
            if self.subcat['SubhaloLenType'][iii, 0] > self.min_n_gas:
                if iii % 1000 == 0:
                    print "processing subhalo "+str(iii)+" out of "+str(self.nsubs)
                data =  illustris_python.snapshot.loadSubhalo(self.basePath, self.snapnum,iii,partType,fields=fields)
                for xyz in np.arange(3):
                    data['Coordinates'][:,xyz] -= self.subcat['SubhaloPos'][iii,xyz]
                    data['Coordinates'][ (data['Coordinates'][:,xyz] < -self.boxsize/2.0) ,xyz]  += self.boxsize
                    data['Coordinates'][ (data['Coordinates'][:,xyz] >  self.boxsize/2.0) ,xyz]  -= self.boxsize

                r = np.sqrt( data['Coordinates'][:,0]**2 + data['Coordinates'][:,1]**2 + data['Coordinates'][:,2]**2 )
                m = data['Masses'] * 1e10 / 0.704 
                surface_density_profile, edges = np.histogram( r, weights=m, range=[0,20], bins=100  )
                surface_density_profile /= 3.14159 * ( edges[1:]*edges[1:] - edges[:-1]*edges[:-1] )
                surface_density_profile  = np.log( surface_density_profile )
                inner_edges = edges[:-1]
 
                index = (inner_edges<5)
                if np.sum( np.isfinite( surface_density_profile[index] ) ) > 10:  result_5 = np.polyfit( inner_edges[index], surface_density_profile[index], 1)
                else: result_5 = np.array([0,0])

                index = (inner_edges<10)
                if np.sum( np.isfinite( surface_density_profile[index] ) ) > 10:  result_10 = np.polyfit( inner_edges[index], surface_density_profile[index], 1)
                else: result_10 = np.array([0,0])

                index = (inner_edges<20)
                if np.sum( np.isfinite( surface_density_profile[index] ) ) > 10:  result_20 = np.polyfit( inner_edges[index], surface_density_profile[index], 1)
                else: result_20 = np.array([0,0])

                #print result_5
                if np.sum(result_5)==0:
                    sigma0_5[iii]  = np.nan
                    sl_5[iii]      = np.nan
                else:
                    sigma0_5[iii]  = result_5[1]
                    sl_5[iii]  = -1.0/result_5[0]

                if np.sum(result_10)==0:
                    sigma0_10[iii]  = np.nan
                    sl_10[iii]      = np.nan
                else:
                    sigma0_10[iii] = result_10[1]
                    sl_10[iii] = -1.0/result_10[0]

                if np.sum(result_20)==0:
                    sigma0_20[iii] = np.nan
                    sl_20[iii]     = np.nan
                else:
                    sigma0_20[iii] = result_20[1]
                    sl_20[iii] = -1.0/result_20[0]

        new_file_name = './catalogs/'+self.run+'/gas_disk_properties/gas_disk_properties_info_'+str(self.snapnum)+'.hdf5'

        print new_file_name
        f = hdf5lib.OpenFile(new_file_name, mode="w")
        subhalo_group = hdf5lib.CreateGroup(f, "Subhalo")
        hdf5lib.CreateArray(f, subhalo_group, 'CentralDensity_5',  sigma0_5)
        hdf5lib.CreateArray(f, subhalo_group, 'CentralDensity_10', sigma0_10)
        hdf5lib.CreateArray(f, subhalo_group, 'CentralDensity_20', sigma0_20)
        hdf5lib.CreateArray(f, subhalo_group, 'ScaleLength_5',  sl_5)
        hdf5lib.CreateArray(f, subhalo_group, 'ScaleLength_10', sl_10)
        hdf5lib.CreateArray(f, subhalo_group, 'ScaleLength_20', sl_20)

        f.close()




    def make_gas_kinematics(self, **kwargs):
        fields=[u'Coordinates', u'Velocities', u'Masses']
        partType=0

        v_1 = np.zeros( self.nsubs )
        v_2 = np.zeros( self.nsubs )
        v_3 = np.zeros( self.nsubs )
        v_4 = np.zeros( self.nsubs )
        v_5 = np.zeros( self.nsubs )
        s_1 = np.zeros( self.nsubs )
        s_2 = np.zeros( self.nsubs )
        s_3 = np.zeros( self.nsubs )
        s_4 = np.zeros( self.nsubs )
        s_5 = np.zeros( self.nsubs )
        
        
        for iii in np.arange(self.nsubs):
          if self.subcat['SubhaloLenType'][iii, 0] > self.min_n_gas:
            if iii % 1000 == 0:
                print "processing subhalo "+str(iii)+" out of "+str(self.nsubs)
            data =  illustris_python.snapshot.loadSubhalo(self.basePath, self.snapnum,iii,partType,fields=fields)
            for xyz in np.arange(3):
                data['Coordinates'][:,xyz] -= self.subcat['SubhaloPos'][iii,xyz]
                data['Coordinates'][ (data['Coordinates'][:,xyz] < -self.boxsize/2.0) ,xyz]  += self.boxsize
                data['Coordinates'][ (data['Coordinates'][:,xyz] >  self.boxsize/2.0) ,xyz]  -= self.boxsize
                data['Velocities'][:,xyz] -= np.median( data['Velocities'][:,xyz] )

            x  = data['Coordinates'][:,0]
            y  = data['Coordinates'][:,1]
            z  = data['Coordinates'][:,2]
            vx = data['Velocities'][:,0]
            vy = data['Velocities'][:,1]
            vz = data['Velocities'][:,2]
            mass = data['Masses']
            r = np.sqrt( data['Coordinates'][:,0]**2 + data['Coordinates'][:,1]**2 + data['Coordinates'][:,2]**2 )

            index = (r<20)
            Lx = np.sum(mass[index] * (y[index]*vz[index] - z[index]*vy[index]))
            Ly = np.sum(mass[index] * (z[index]*vx[index] - x[index]*vz[index]))
            Lz = np.sum(mass[index] * (x[index]*vy[index] - y[index]*vx[index]))
            phi   = np.arctan2(Ly,Lx)
            theta = np.arctan2(np.sqrt(Lx**2+Ly**2),Lz)

            x_ = -z  * np.sin(theta)   + (x *  np.cos(phi) + y * np.sin(phi)) * np.cos(theta)
            y_ = -x  * np.sin(phi)     +  y  * np.cos(phi)
            z_ =  z  * np.cos(theta)   + (x *  np.cos(phi) + y * np.sin(phi)) * np.sin(theta)
            vx_ = -vz  * np.sin(theta) + (vx * np.cos(phi) + vy *np.sin(phi)) * np.cos(theta)
            vy_ = -vx  * np.sin(phi)   +  vy * np.cos(phi)
            vz_ =  vz  * np.cos(theta) + (vx * np.cos(phi) + vy *np.sin(phi)) * np.sin(theta)

            x = x_
            y = y_
            z = z_
            vx = vx_
            vy = vy_
            vz = vz_

            r2 = np.sqrt( x**2 + y**2)
            vr = ((x * vy ) - (y * vx)) / r2

            nbins=5
            dummy_vr= np.zeros(nbins)
            dummy_sigma = np.zeros(nbins)

            for ijk in np.arange(nbins):
                min_r = (0.0+1.0*ijk)/(nbins) * 25.0
                max_r = (1.0+1.0*ijk)/(nbins) * 25.0
                index = (r2 > min_r) & (r2 < max_r) & (z > -20) & (z < 20 )
                if np.sum(index) > 5:
                    dummy_vr[ijk] = np.median( vr[index] )
                    dummy_sigma[ijk] = np.std( vr[index] )

            v_1[iii] = dummy_vr[0]
            v_2[iii] = dummy_vr[1]
            v_3[iii] = dummy_vr[2]
            v_4[iii] = dummy_vr[3]
            v_5[iii] = dummy_vr[4]

            s_1[iii] = dummy_sigma[0]
            s_2[iii] = dummy_sigma[1]
            s_3[iii] = dummy_sigma[2]
            s_4[iii] = dummy_sigma[3]
            s_5[iii] = dummy_sigma[4]


        new_file_name = './catalogs/'+self.run+'/gas_kinematics/gas_kinematic_info_'+str(self.snapnum)+'.hdf5'
        f = hdf5lib.OpenFile(new_file_name, mode="w")
        subhalo_group = hdf5lib.CreateGroup(f, "Subhalo")
        hdf5lib.CreateArray(f, subhalo_group, 'v_5',  v_1)
        hdf5lib.CreateArray(f, subhalo_group, 'v_10', v_2)
        hdf5lib.CreateArray(f, subhalo_group, 'v_15', v_3)
        hdf5lib.CreateArray(f, subhalo_group, 'v_20', v_4)
        hdf5lib.CreateArray(f, subhalo_group, 'v_25', v_5)
        hdf5lib.CreateArray(f, subhalo_group, 's_5',  s_1)
        hdf5lib.CreateArray(f, subhalo_group, 's_10', s_2)
        hdf5lib.CreateArray(f, subhalo_group, 's_15', s_3)
        hdf5lib.CreateArray(f, subhalo_group, 's_20', s_4)
        hdf5lib.CreateArray(f, subhalo_group, 's_25', s_5)

        
        f.close()

