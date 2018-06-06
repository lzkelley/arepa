import glob
import numpy as np
import simread.readsnapHDF5 as ws
import units.springel_units as units
import os.path
import util.hdf5lib as hdf5lib


class merger_time_series:
    def __init__(self, path, tag, factor=1, **kwargs):
       
        snap_list = np.array(glob.glob(path+'/snap*_*hdf5'))
        print path, snap_list
        all_snapnrs = np.sort( np.array([int( file[file.index('snapshot_')+9:file.index('.hdf5')] ) for file in snap_list], dtype=int) )
        first_snapnr = all_snapnrs[0]
        new_snap_list = [None] * snap_list.shape[0]
        snap_timing  = np.zeros(snap_list.shape[0])
        for index in np.arange(snap_list.shape[0]):
                this_file = path+"snapshot_"+str(index*factor+int(first_snapnr)).zfill(3)+".hdf5"
                print this_file
                new_snap_list[index] = this_file
		try:
                    head = ws.snapshot_header(this_file[:this_file.index('.hdf5')])
                    snap_timing[index] = head.time
		except:
	            snap_timing[index] = -1.0

        self.first_snapnr	= first_snapnr
        self.savetag            = tag
        self.snap_list          = np.array(new_snap_list)
        self.snap_timing        = snap_timing
        self.n_snaps 		= self.snap_list.shape[0]
        
        print " "
        print tag
        print " "

        self.alloc_loc_arrays()

        this_file = path+"snapshot_"+str(int(first_snapnr)).zfill(3)+".hdf5"
        print this_file
        gas_xyz     = np.array(ws.read_block(this_file,'POS ', parttype=0))
        star_xyz    = np.append(
                        np.array(ws.read_block(this_file,'POS ', parttype=2)),
                        np.array(ws.read_block(this_file,'POS ', parttype=3)),
                                axis=0)
        star_mass   = np.append(
                        np.array(ws.read_block(this_file,'MASS', parttype=2)),
                        np.array(ws.read_block(this_file,'MASS', parttype=3)),
                        axis=0) * 1e10 


        bh_id       = np.array(ws.read_block(this_file,'ID  ', parttype=5))
        bh_xyz      = np.array(ws.read_block(this_file,'POS ', parttype=5))

        bh_order = np.argsort(bh_id)
        bh_xyz = bh_xyz[bh_order,:]


#        bh_xyz      = np.array(ws.read_block(this_file,'POS ', parttype=5))

	# calculate stellar offsets, to determine each galaxy's initial stellar mass (for MZ rleaiton)
        dr = np.zeros( (bh_xyz.shape[0], star_xyz.shape[0]) )
        for iii in np.arange(bh_xyz.shape[0]):
            dx2 = ( star_xyz[:,0] - bh_xyz[iii,0] )**2
            dy2 = ( star_xyz[:,1] - bh_xyz[iii,1] )**2
            dz2 = ( star_xyz[:,2] - bh_xyz[iii,2] )**2
            dr[iii,:]  = np.sqrt( dx2 + dy2 + dz2 )             # offset from this BH; and then you have it.
        min_dr = np.min(dr,axis=0)              # now the min_dr, regardless of n_bh

        print bh_xyz
        print dr

        print dr.shape
        print min_dr.shape

        m_gal_1 = 0
        m_gal_2 = 0
        for iii in np.arange( min_dr.shape[0] ):
            if min_dr[iii]==dr[0,iii]:
                m_gal_1 += star_mass[iii]
            else:
                m_gal_2 += star_mass[iii]
        print m_gal_1, m_gal_2
        print np.log10(m_gal_1), np.log10(m_gal_2)
        x1 = np.log10(m_gal_1)
        x2 =  np.log10(m_gal_2)
        z10_tmp = 27.7911 - 6.94493 * x1 + 0.808097 * x1 * x1 - 0.0301508 * x1 * x1 * x1
        z20_tmp = 27.7911 - 6.94493 * x2 + 0.808097 * x2 * x2 - 0.0301508 * x2 * x2 * x2


        print "  Initial Metallicity Summary: "
        print "  We think galaxy 1 has a mass of {:.3f} and a metallicity of {:.3f}".format( x1, z10_tmp )
        print "  We think galaxy 2 has a mass of {:.3f} and a metallicity of {:.3f}".format( x2, z20_tmp )


	# calculate gas offsets, to determine each galaxy's metallicity gradient	
        dr = np.zeros( (bh_xyz.shape[0], gas_xyz.shape[0]) )
        for iii in np.arange(bh_xyz.shape[0]):
            dx2 = ( gas_xyz[:,0] - bh_xyz[iii,0] )**2
            dy2 = ( gas_xyz[:,1] - bh_xyz[iii,1] )**2
            dz2 = ( gas_xyz[:,2] - bh_xyz[iii,2] )**2
            dr[iii,:]  = np.sqrt( dx2 + dy2 + dz2 )             # offset from this BH; and then you have it.
        min_dr = np.min(dr,axis=0)              # now the min_dr, regardless of n_bh, (right?!)

        tmp_ids = np.array(ws.read_block(this_file,'ID  ',parttype=0))

        self.init_ids         = np.zeros( tmp_ids.max()+1 )
        self.init_metallicity = np.zeros( tmp_ids.max()+1 )
        for iii in np.arange( min_dr.shape[0] ):
            if min_dr[iii]==dr[0, iii]:
                x = np.log10(m_gal_1)
            else:
                x = np.log10(m_gal_2)
            oh_12 = 27.7911 - 6.94493 * x + 0.808097 * x * x - 0.0301508 * x * x * x
            z_zsolar    = 10.0**(oh_12 - 12) / 0.0004	# central metallicity in solar units.

            this_id = tmp_ids[iii]
            this_z  = 0.0127 * z_zsolar * 10.0**( - 0.1 * min_dr[iii] )
            self.init_ids[this_id]		= this_id
            self.init_metallicity[this_id]      = this_z
	    # print "{:.3f}  {:.3f}  {:.2f} \n".format( np.log10( this_z * 0.0004 / 0.0127) + 12 , oh_12, min_dr[iii] )


    def alloc_loc_arrays(self):
        self.rad_bins = np.array( [1.0, 2, 3, 5, 20] )
        
        self.values = {}
        for global_value in ['time', 'sep']:
            self.values[global_value] = np.zeros( self.n_snaps )
        
        for radial_value in ['sfr', 'z', 'mgas', 'sfr_sm', 'sfr_sm2',
                             'stellar_ages', 'mstar', 'm_molec_gas',
                             'm_atomic_gas', 'm_warm_gas', 'm_hot_gas',
                             'bh_acc_rate', 'torque_acc_rate']:
            self.values[radial_value] = np.zeros( (2, 5, self.n_snaps) )


    def write_one_snapshot(self, snap_index, **kwargs):
        this_file   = self.snap_list[snap_index]

        header = ws.snapshot_header( this_file )
        self.values['time'][ snap_index ] = header.time

        print this_file
        gas_xyz     = np.array(ws.read_block(this_file,'POS ', parttype=0))
        gas_sfr     = np.array(ws.read_block(this_file,'SFR ', parttype=0))
        gas_z	    = np.array(ws.read_block(this_file,'Z   ', parttype=0))
        gas_z       = gas_z[:,0]
        gas_mass    = np.array(ws.read_block(this_file,'MASS', parttype=0))
        gas_ids	    = np.array(ws.read_block(this_file,'ID  ',parttype=0))

        gas_u       = np.array(ws.read_block(this_file,'U   ',parttype=0))
        gas_ne      = np.array(ws.read_block(this_file,'NE  ',parttype=0))
        gas_t	    = units.gas_code_to_temperature( gas_u, gas_ne)
        
        


        for index, value in enumerate(gas_ids):	# IDs, add init metallicity
            #print index, value, gas_z.shape, self.init_metallicity.shape
            try:
                gas_z[index] += self.init_metallicity[ value ]              # all that matters here is that init metallicty can be indexed by ID
            except:
		try:
                    gas_z[index] += self.init_metallicity[ value & 2**30-1 ]
		except:
		    gas_z[index] += 0.0	    # confused...
		    sys.exit()
    

        try:
            star_xyz    = np.array(ws.read_block(this_file,'POS ', parttype=4))
            star_mass   = np.array(ws.read_block(this_file,'MASS', parttype=4))
            star_ages   = np.array(ws.read_block(this_file,'AGE ', parttype=4)) 
            star_ages   = header.time - star_ages
        except:
            star_xyz = np.zeros( (10,3) )
            star_mass = np.zeros( 10 )
            star_ages = np.zeros( 10 )

        bh_id       = np.array(ws.read_block(this_file,'ID  ', parttype=5))
        bh_xyz	    = np.array(ws.read_block(this_file,'POS ', parttype=5))

        bh_order = np.argsort(bh_id)
        bh_xyz = bh_xyz[bh_order,:]

        if bh_xyz.shape[0] == 1:
            self.values['sep'][ snap_index ] = 0
        else:
            self.values['sep'][  snap_index ] = np.sqrt(
                                (bh_xyz[0,0] - bh_xyz[1,0])**2 +
                                (bh_xyz[0,1] - bh_xyz[1,1])**2 +
                                (bh_xyz[0,2] - bh_xyz[1,2])**2 )

        dr = np.zeros( (bh_xyz.shape[0], gas_xyz.shape[0]) )
        for iii in np.arange(bh_xyz.shape[0]):
            dx2 = ( gas_xyz[:,0] - bh_xyz[iii,0] )**2
            dy2 = ( gas_xyz[:,1] - bh_xyz[iii,1] )**2
            dz2 = ( gas_xyz[:,2] - bh_xyz[iii,2] )**2
            dr[iii,:]  = np.sqrt( dx2 + dy2 + dz2 )		# offset from this BH; and then you have it.
	
        min_dr = np.min(dr,axis=0)		# now the min_dr, regardless of n_bh, (right?!)

        all_molec_gas_mass = gas_mass[    gas_t < 100 ]
        all_molec_gas_min_dr    = min_dr[ gas_t < 100 ]
        all_molec_gas_dr        = dr[  :, gas_t < 100 ]

        all_atomic_gas_mass = gas_mass[    (gas_t > 100) & (gas_t < 1000) ]
        all_atomic_gas_min_dr    = min_dr[ (gas_t > 100) & (gas_t < 1000) ]
        all_atomic_gas_dr        = dr[  :, (gas_t > 100) & (gas_t < 1000) ]

        all_warm_gas_mass = gas_mass[    (gas_t > 1000) & (gas_t < 1e5) ]
        all_warm_gas_min_dr    = min_dr[ (gas_t > 1000) & (gas_t < 1e5) ]
        all_warm_gas_dr        = dr[  :, (gas_t > 1000) & (gas_t < 1e5) ]

        all_hot_gas_mass = gas_mass[    (gas_t > 1e5) ]
        all_hot_gas_min_dr    = min_dr[ (gas_t > 1e5) ]
        all_hot_gas_dr        = dr[  :, (gas_t > 1e5) ]
 
        for iii in np.arange(bh_xyz.shape[0]):
            for jjj in np.arange(self.rad_bins.shape[0]):
                this_index = (dr[iii,:] < self.rad_bins[jjj]) & (np.equal(dr[iii,:], min_dr))
                self.values['sfr'][iii, jjj, snap_index] = np.sum(     gas_sfr[ this_index] )
                self.values['z'][  iii, jjj, snap_index] = np.median(  gas_z[ this_index]   )
                self.values['mgas'][iii,jjj, snap_index] = np.sum(     gas_mass[ this_index ] )

		
                this_index = (all_molec_gas_dr[iii,:] < self.rad_bins[jjj]) & (np.equal(all_molec_gas_dr[iii,:],  all_molec_gas_min_dr ) )
                self.values['m_molec_gas'][iii, jjj, snap_index] = np.sum( all_molec_gas_mass[this_index] )

                this_index = (all_atomic_gas_dr[iii,:] < self.rad_bins[jjj]) & (np.equal(all_atomic_gas_dr[iii,:],  all_atomic_gas_min_dr ) )
                self.values['m_atomic_gas'][iii, jjj, snap_index] = np.sum( all_atomic_gas_mass[this_index] )

                this_index = (all_warm_gas_dr[iii,:] < self.rad_bins[jjj]) & (np.equal(all_warm_gas_dr[iii,:],  all_warm_gas_min_dr ) )
                self.values['m_warm_gas'][iii, jjj, snap_index] = np.sum( all_warm_gas_mass[this_index] )

                this_index = (all_hot_gas_dr[iii,:] < self.rad_bins[jjj]) & (np.equal(all_hot_gas_dr[iii,:],  all_hot_gas_min_dr ) )
                self.values['m_hot_gas'][iii, jjj, snap_index] = np.sum( all_hot_gas_mass[this_index] )


        dr = np.zeros( (bh_xyz.shape[0], star_xyz.shape[0]) )
        for iii in np.arange(bh_xyz.shape[0]):
            dx2 = ( star_xyz[:,0] - bh_xyz[iii,0] )**2
            dy2 = ( star_xyz[:,1] - bh_xyz[iii,1] )**2
            dz2 = ( star_xyz[:,2] - bh_xyz[iii,2] )**2
            dr[iii,:]  = np.sqrt( dx2 + dy2 + dz2 )             # offset from this BH; and then you have it.

        min_dr = np.min(dr,axis=0)              # now the min_dr, regardless of n_bh, (right?!)
        for iii in np.arange(bh_xyz.shape[0]):
            for jjj in np.arange(self.rad_bins.shape[0]):
                this_index = (dr[iii,:] < self.rad_bins[jjj]) & (np.equal(dr[iii,:], min_dr))
                if np.sum(this_index) > 0:
                    self.values['stellar_ages'][  iii, jjj, snap_index] = np.median(  star_ages[ this_index]   )
                    self.values['mstar'][iii,jjj, snap_index]  = np.sum(     star_mass[ this_index ] )

                this_index = (dr[iii,:] < self.rad_bins[jjj]) & (np.equal(dr[iii,:], min_dr)) & ( star_ages < 0.01 )
                if np.sum(this_index) > 0:
                    self.values['sfr_sm'][iii,jjj, snap_index] = np.sum(     star_mass[this_index] ) / 0.01 * 10.0

                this_index = (dr[iii,:] < self.rad_bins[jjj]) & (np.equal(dr[iii,:], min_dr)) & ( star_ages < 0.1 )
                if np.sum(this_index) > 0:
                    self.values['sfr_sm2'][iii,jjj, snap_index] = np.sum(     star_mass[this_index] ) / 0.1 * 10.0


    def write_all_snapshots(self, rank=0, size=1, **kwargs):
        for index in np.arange(self.n_snaps):
            if (index % size) == rank:
		try:
                    self.write_one_snapshot(index, **kwargs)
		except:
                    print "snapshot "+str(index)+" failed"

    def collect_at_root(self, comm, write=True):
        for name in self.values:
            tmp = np.zeros_like ( self.values[name] )
            comm.Allreduce( self.values[name], tmp)
            self.values[name] = np.copy( tmp )

        if (comm.Get_rank() == 0) and write==True:
            f_write = hdf5lib.OpenFile('./saved_data/'+self.savetag+'.hdf5', mode="w")
            group = hdf5lib.CreateGroup(f_write, "data")
            for name in self.values:
                hdf5lib.CreateArray(f_write, group, name, self.values[name])
            f_write.close()
                
                
#            3self.my_write( )
            
            #            for name in self.values:
            #    if name in ['time', 'sep']:
            #        self.my_write(name, self.values[name], oned=True)
            #    else:
#        self.my_write(name, self.values[name])


    def my_write(self, loc_tag, loc_data, oned=False):
        if not os.path.exists('./saved_data'): os.mkdir('./saved_data')
        if not os.path.exists('./saved_data/'+self.savetag): os.mkdir('./saved_data/'+self.savetag)
        
        
        f = open('./saved_data/'+self.savetag+'/'+loc_tag+'.txt','w')
        for iii in np.arange(self.n_snaps):
            this_row_txt = ""
            if oned:
                this_row_txt = "  {:.8f}  ".format( loc_data[iii] )
            else:
                for jjj in np.arange(self.rad_bins.shape[0]):
                    this_row_txt = this_row_txt+"  {:.8f}  {:.8f}  ".format( loc_data[0,jjj,iii], loc_data[1,jjj,iii])
            f.write( this_row_txt+" \n")
        f.close()




