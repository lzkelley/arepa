import numpy as np
import math
import ctypes
import util.utilities as util
import gadget_lib.gadget as gadget

def checklen(x):
    return len(np.array(x,ndmin=1));

def int_round(x):
    return np.int(np.round(x));

def ok_scan(input,xmax=1.0e10,pos=0):
    if (pos==1):
        return (np.isnan(input)==False) & (np.fabs(input)<=xmax) & (input > 0.);
    else:
        return (np.isnan(input)==False) & (np.fabs(input)<=xmax);

def fcor(x):
    return np.array(x,dtype='f',ndmin=1)
def vfloat(x):
    return x.ctypes.data_as(ctypes.POINTER(ctypes.c_float));



def calculate_zoom_center(sdir,snum,cen=[0.,0.,0.],clip_size=2.e10):
    rgrid=np.array([1.0e10,1000.,700.,500.,300.,200.,100.,70.,50.,30.,20.,10.,5.,2.5,1.]);
    rgrid=rgrid[rgrid <= clip_size];
    Ps=gadget.readsnap(sdir,snum,4,cosmological=1);
    n_new=Ps['m'].shape[0];
    if (n_new > 1):
        pos=Ps['p']; x0s=pos[:,0]; y0s=pos[:,1]; z0s=pos[:,2];
    Pg=gadget.readsnap(sdir,snum,0,cosmological=1);
    rho=np.array(Pg['rho'])*407.5;
    if (rho.shape[0] > 0):
        pos=Pg['p']; x0g=pos[:,0]; y0g=pos[:,1]; z0g=pos[:,2];
    rho_cut=1.0e-5;
    cen=np.array(cen);

    for i_rcut in range(len(rgrid)):
        for j_looper in range(5):
            if (n_new > 1000):
                x=x0s; y=y0s; z=z0s;
            else:
                ok=(rho > rho_cut);
                x=x0g[ok]; y=y0g[ok]; z=z0g[ok];
            x=x-cen[0]; y=y-cen[1]; z=z-cen[2];
            r = np.sqrt(x*x + y*y + z*z);
            ok = (r < rgrid[i_rcut]);
            if (len(r[ok]) > 1000):
                x=x[ok]; y=y[ok]; z=z[ok];
                if (i_rcut <= len(rgrid)-5):
                    cen+=np.array([np.median(x),np.median(y),np.median(z)]);
                else:
                    cen+=np.array([np.mean(x),np.mean(y),np.mean(z)]);
            else:
                if (len(r[ok]) > 200):
                    x=x[ok]; y=y[ok]; z=z[ok];
                    cen+=np.array([np.mean(x),np.mean(y),np.mean(z)]);
                    
    return cen;
    
    
def test():
    import gadget
    sdir='../zooms/m10_dw/'
    snum=440
    cen=calculate_zoom_center(sdir,snum);
    cen=np.array([ 3306.55575601,  3029.35099489,  3235.4757624 ]);
    print cen;

    PPPs=gadget.readsnap(sdir,snum,4,cosmological=1);
    PPP=gadget.readsnap(sdir,snum,0,cosmological=1);
    source_pos=(np.random.rand(3,5)-0.5)*100.
    source_pos[2,:] *= 0.02
    for i in [0,1,2]: source_pos[i,:]+=cen[i];
    source_pos=PPPs['p']
    
    mm=PPP['m'];
    rho=PPP['rho'];
    
    #mm *= PPP['nh'];
    #rho *= PPP['nh'];
    #print np.min(PPP['nh']), np.median(PPP['nh']), np.max(PPP['nh'])
    
    
    
    xlen = 10.
    los_nh,los_nh_hot,los_z = \
      return_columns_to_sources(source_pos,\
      PPP['p'],PPP['u'],rho,PPP['h'],PPP['ne'],PPP['nh'],PPP['z'],mm,\
      xrange=cen[0]+[-xlen,xlen],yrange=cen[1]+[-xlen,xlen],zrange=cen[2]+[-xlen,xlen])

    print los_nh
    print los_nh_hot
    print los_z



##
## return: los_NH_allgas, los_NH_hotphase, los_gas_metallicity 
##
def return_columns_to_sources( source_pos, gas_pos, \
    gas_u, gas_rho, gas_hsml, gas_numh, gas_nume, gas_metallicity, gas_mass, \
    xrange=0, yrange=0, zrange=0, \
    MIN_CELL_SIZE=0.01, OUTER_RANGE_OF_INT=1200., \
    TRIM_PARTICLES=1 ):
    
    ## check the ordering of the position matrices:
    if ((checklen(gas_pos[0,:])==3) & (checklen(gas_pos[:,0]) !=3)): gas_pos=np.transpose(gas_pos);
    if ((checklen(source_pos[0,:])==3) & (checklen(source_pos[:,0]) !=3)): source_pos=np.transpose(source_pos);
    ## and that metallicities are a vector, not a matrix
    if (len(gas_metallicity.shape)>1): gas_metallicity=gas_metallicity[:,0]

    if ((checklen(gas_pos[:,0]) != 3) | (checklen(gas_pos[0,:]) <= 1)):
        print 'ERROR WILL OCCUR :: need pos to be (3,N)'

    x=source_pos[0,:] ; y=source_pos[1,:] ; z=source_pos[2,:]
    if(checklen(xrange)<=1): xrange=[np.min(x),np.max(x)];
    if(checklen(yrange)<=1): yrange=[np.min(y),np.max(y)];
    xr=xrange; yr=yrange;
    if(checklen(zrange)<=1):
        zrr=np.sqrt((xr[1]-xr[0])**2.+(yr[1]-yr[0])**2.)/np.sqrt(2.);
        zmin=np.median(z)-zrr; zmax=np.median(z)+zrr;
        if (np.min(z) > zmin): zmin=np.min(z);
        zrange=[zmin,zmax]; print 'z_range (calc) == ',zrange
    zr=zrange;
    x00=0.5*(xr[1]+xr[0]); y00=0.5*(yr[1]+yr[0]); z00=0.5*(zr[1]+zr[0]); 
    tolfac = 1.0e10;
    if (TRIM_PARTICLES==1):
        tolfac = 0.05; 
        #tolfac = -0.01;
        ## trim down the incoming list to only whats in the range plotted 
        ##   (saves a ton of time and memory overflow crashes)

    dx=(0.5+tolfac)*(xr[1]-xr[0]); dy=(0.5+tolfac)*(yr[1]-yr[0]); dz=(0.5+tolfac)*(zr[1]-zr[0]);
    ok_sources=ok_scan(x-x00,xmax=dx) & ok_scan(y-y00,xmax=dy) & ok_scan(z-z00,xmax=dz);
    x=gas_pos[0,:] ; y=gas_pos[1,:] ; z=gas_pos[2,:]
    gw=gas_rho ; gh=gas_hsml ; gz=gas_metallicity ; gm=gas_mass
    ok_gas=ok_scan(x-x00,xmax=dx) & ok_scan(y-y00,xmax=dy) & ok_scan(z-z00,xmax=dz) & \
        ok_scan(gw,pos=1) & ok_scan(gh,pos=1) & ok_scan(gz,pos=1) & ok_scan(gm,pos=1,xmax=1.0e40);

    Ngas = checklen(gas_mass[ok_gas]);
    Nstars = checklen(source_pos[0,ok_sources]);
    if (Nstars<=1) or (Ngas<=1):
        print ' UH-OH: EXPECT ERROR NOW, there are no valid source/gas particles to send!'
        print 'Ngas=',Ngas,'Nstars=',Nstars,'dx=',dx,'dy=',dy,'dz=',dz,'x00=',x00,'y00=',y00,'z00=',z00
        return -1,-1,-1;

    dzmax=np.max(gas_pos[2,ok_gas])-z00; 
    if(dzmax<OUTER_RANGE_OF_INT): OUTER_RANGE_OF_INT=dzmax;
    print 'PASSING: N_gas=',Ngas,'N_sources=',Nstars,'MaxDist=',OUTER_RANGE_OF_INT,'MinCell=',MIN_CELL_SIZE;
    Nbh=0; theta=1.0e-4; phi=1.0e-4;
  
    ## load the routine we need
    exec_call=util.return_python_routines_cdir()+'/LOS_column_singlePOV/getnh.so'
    NH_routine=ctypes.cdll[exec_call];
    ## cast the variables to store the results
    nh_out_cast=ctypes.c_float*Nstars; 
    los_NH_out=nh_out_cast(); los_NH_hot_out=nh_out_cast(); los_Z_out=nh_out_cast();

    ## ok this is a bit arcane but the routine will read appropriately this block order
    Coord = np.zeros((Ngas+Nstars,10),dtype='f');
    Coord[0:Ngas,0] = gas_pos[0,ok_gas]-x00;
    Coord[0:Ngas,1] = gas_pos[1,ok_gas]-y00;
    Coord[0:Ngas,2] = gas_pos[2,ok_gas]-z00;
    Coord[0:Ngas,3] = gas_u[ok_gas]
    Coord[0:Ngas,4] = gas_rho[ok_gas]
    Coord[0:Ngas,5] = gas_hsml[ok_gas]
    Coord[0:Ngas,6] = gas_numh[ok_gas]
    Coord[0:Ngas,7] = gas_nume[ok_gas]
    Coord[0:Ngas,8] = gas_metallicity[ok_gas]
    Coord[0:Ngas,9] = gas_mass[ok_gas]
    Coord[Ngas:Nstars+Ngas,0] = source_pos[0,ok_sources]-x00;
    Coord[Ngas:Nstars+Ngas,1] = source_pos[1,ok_sources]-y00;
    Coord[Ngas:Nstars+Ngas,2] = source_pos[2,ok_sources]-z00;
    Coord=np.copy(np.transpose(Coord));

    ## main call to the NH-calculation routine
    NH_routine.getnh( ctypes.c_int(Ngas), ctypes.c_int(Nstars), ctypes.c_int(Nbh), \
        ctypes.c_float(theta), ctypes.c_float(phi), \
        vfloat(Coord), \
        ctypes.byref(los_NH_out),  ctypes.byref(los_NH_hot_out),  ctypes.byref(los_Z_out), \
        ctypes.c_float(OUTER_RANGE_OF_INT), ctypes.c_float(MIN_CELL_SIZE) );
    ## now put the output arrays into a useful format 
    print type(los_NH_out), los_NH_out
    los_NH = np.ctypeslib.as_array(np.copy(los_NH_out));
    los_NH_hot = np.ctypeslib.as_array(np.copy(los_NH_hot_out));
    los_Z = np.ctypeslib.as_array(np.copy(los_Z_out));

    # trap for really low NH value and zero metallicity (make it small instead)
    low_NH = 1.0e10;
    los_NH[los_NH<low_NH]=low_NH; los_NH_hot[los_NH_hot<low_NH]=low_NH;
    los_Z[los_Z<=1.0e-5]=1.0e-5;

    ## assign strong attenuation to all 'off-grid' sources, then fill in calc. vals
    Nstarstot=checklen(source_pos[0,:]);
    los_NH_allgas=np.zeros(Nstarstot,dtype='f')+1.0e23;
    los_NH_hotgas=np.zeros(Nstarstot,dtype='f')+1.0e23;
    los_gas_metallicity=np.zeros(Nstarstot,dtype='f')+0.02;
    nok=checklen(los_NH_allgas[ok_sources])
    los_NH_allgas[ok_sources]=fcor(los_NH[0:Nstars]);
    los_NH_hotgas[ok_sources]=fcor(los_NH_hot[0:Nstars]);
    los_gas_metallicity[ok_sources]=fcor(los_Z[0:Nstars]);

    return los_NH_allgas, los_NH_hotgas, los_gas_metallicity;
