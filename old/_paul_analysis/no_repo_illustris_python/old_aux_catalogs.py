import illustris_python
import numpy as np
import util.hdf5lib as hdf5lib
import h5py
import sys

base='/n/ghernquist/Illustris/Runs/'
#base='/n/hernquistfs3/IllustrisTNG/Runs/'
#'L75n1820TNG'
scalefactor_list = np.loadtxt('/n/home01/ptorrey/IllustrisAuxFiles/Illustris_scalefactors.txt')
redshift_list    = 1.0/scalefactor_list - 1.0


class aux_catalogs:
    def __init__(self, snapnum, run, **kwargs):

        self.run = 'Illustris-'+str(run)+'/'
#        self.run = run
        self.basePath = base+self.run+'/output/'
        print self.basePath
        self.boxsize=75000.0
	self.little_h = 0.704
        self.min_n_stars = 500
        self.min_n_gas = 500
        self.snapnum = snapnum

	all_scalefactors = np.loadtxt('/n/home01/ptorrey/IllustrisAuxFiles/Illustris_scalefactors.txt')
        self.scalefactor = all_scalefactors[snapnum]

        fields      = [u'SubhaloPos', u'SubhaloVel', u'SubhaloHalfmassRadType', u'SubhaloLenType']
        #if "TNG" in self.run:
        self.subcat = illustris_python.cat.loadSubhalos(snapnum, run=self.run, fields=fields)
        #else:
        #    self.subcat = illustris_python.cat.loadSubhalos(snapnum, run=self.run, fields=fields)

        self.nsubs  = self.subcat['count']

    def make_sz_aux(self, **kwargs):
        fields=[u'Coordinates', u'Masses', u'GFM_StellarFormationTime']
        partType=4
        stored_core_1     = np.zeros( self.nsubs ) - 1.0
        stored_core_2     = np.zeros( self.nsubs ) - 1.0

        for iii in np.arange(self.nsubs):
          if self.subcat['SubhaloLenType'][iii, 4] > self.min_n_stars:
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

        new_file_name = './catalogs/'+self.run+'stellar_core_mass/stellar_core_mass_'+str(self.snapnum)+'.hdf5'
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


        new_file_name = './catalogs/'+self.run+'rotation_curve/rotation_curve_'+str(self.snapnum)+'.hdf5'
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



        new_file_name = './catalogs/'+self.run+'stellar_vel_disp/stellar_vel_disp_'+str(self.snapnum)+'.hdf5'
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

        new_file_name = './catalogs/'+self.run+'stellar_core_mass/stellar_core_mass_'+str(self.snapnum)+'.hdf5'
        print new_file_name
        f_write = hdf5lib.OpenFile(new_file_name, mode="w")
        subhalo_group = hdf5lib.CreateGroup(f_write, "Subhalo")
        hdf5lib.CreateArray(f_write, subhalo_group, 'StellarCoreMass1'  , stored_core_1)
        hdf5lib.CreateArray(f_write, subhalo_group, 'StellarCoreMass2'	, stored_core_2)
        f_write.close()



    def make_gas_metallicity_gradient(self, **kwargs):
        fields=[u'Coordinates', u'GFM_Metallicity']
        partType=0
        z0_5 = np.zeros( self.nsubs )   - 1.0
        z0_10 = np.zeros( self.nsubs )  - 1.0
        z0_20 = np.zeros( self.nsubs )  - 1.0
        dz_5 = np.zeros( self.nsubs )  - 99.0
        dz_10 = np.zeros( self.nsubs ) - 99.0
        dz_20 = np.zeros( self.nsubs ) - 99.0
        
        for iii in np.arange(self.nsubs):
            if self.subcat['SubhaloLenType'][iii, 0] > self.min_n_gas:
                print "processing subhalo "+str(iii)+" out of "+str(self.nsubs)
                data =  illustris_python.snapshot.loadSubhalo(self.basePath, self.snapnum,iii,partType,fields=fields)
                for xyz in np.arange(3):
                    data['Coordinates'][:,xyz] -= self.subcat['SubhaloPos'][iii,xyz]
                    data['Coordinates'][ (data['Coordinates'][:,xyz] < -self.boxsize/2.0) ,xyz]  += self.boxsize
                    data['Coordinates'][ (data['Coordinates'][:,xyz] >  self.boxsize/2.0) ,xyz]  -= self.boxsize
        
                r = np.sqrt( data['Coordinates'][:,0]**2 + data['Coordinates'][:,1]**2 + data['Coordinates'][:,2]**2 )
                met = data['GFM_Metallicity']

                index = (r<5)
                if (np.sum(index) > 30 ): result_5 = np.polyfit(r[r<5],  met[r<5],  1)
                else: result_5 = np.array([0,0])
                index = (r<10)
                if (np.sum(index) > 30 ): result_10= np.polyfit(r[r<10], met[r<10], 1)
                else: result_10 = np.array([0,0])
                index = (r<20)
                if (np.sum(index) > 30): result_20= np.polyfit(r[r<20], met[r<20], 1)
                else: result_20 = np.array([0,0])

                z0_5[iii]  = result_5[1]
                z0_10[iii] = result_10[1]
                z0_20[iii] = result_20[1]

                dz_5[iii]  = result_5[0]
                dz_10[iii] = result_10[0]
                dz_20[iii] = result_20[0]

        new_file_name = './catalogs/'+self.run+'gas_metallicity/gas_metallicity_info_'+str(self.snapnum)+'.hdf5'

        f = hdf5lib.OpenFile(new_file_name, mode="w")
        subhalo_group = hdf5lib.CreateGroup(f, "Subhalo")
        hdf5lib.CreateArray(f, subhalo_group, 'CentralMetallicity_5',  z0_5)
        hdf5lib.CreateArray(f, subhalo_group, 'CentralMetallicity_10', z0_10)
        hdf5lib.CreateArray(f, subhalo_group, 'CentralMetallicity_20', z0_20)
        hdf5lib.CreateArray(f, subhalo_group, 'GradMetallicity_5',  dz_5)
        hdf5lib.CreateArray(f, subhalo_group, 'GradMetallicity_10', dz_10)
        hdf5lib.CreateArray(f, subhalo_group, 'GradMetallicity_20', dz_20)
        
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

            for iii in np.arange(nbins):
                min_r = (0.0+1.0*iii)/(nbins) * 25.0
                max_r = (1.0+1.0*iii)/(nbins) * 25.0
                index = (r2 > min_r) & (r2 < max_r) & (z > -20) & (z < 20 )
                if np.sum(index) > 5:
                    dummy_vr[iii] = np.median( vr[index] )
                    dummy_sigma[iii] = np.std( vr[index] )

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

        new_file_name = './catalogs/'+self.run+'gas_kinematics/gas_kinematic_info_'+str(self.snapnum)+'.hdf5'
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

