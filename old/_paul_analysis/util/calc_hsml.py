import util
import util.cast
import numpy as np
import ctypes


def get_particle_hsml( x, y, z, DesNgb=32, Hmax=0. , **kwargs ):
    x=util.cast.fcor(x); y=util.cast.fcor(y); z=util.cast.fcor(z); N=util.cast.checklen(x);
    ok=(util.cast.ok_scan(x) & util.cast.ok_scan(y) & util.cast.ok_scan(z)); x=x[ok]; y=y[ok]; z=z[ok];
    if(Hmax==0.):
        dx=np.max(x)-np.min(x); dy=np.max(y)-np.min(y); dz=np.max(z)-np.min(z); ddx=np.max([dx,dy,dz]);
        Hmax=5.*ddx*(np.float(N)**(-1./3.)); ## mean inter-particle spacing
    
    ## load the routine we need
    exec_call=util.dir.c_routines_dir()+'/StellarHsml/starhsml.so'
    h_routine=ctypes.cdll[exec_call];

    h_out_cast=ctypes.c_float*N; H_OUT=h_out_cast();
    ## main call to the hsml-finding routine
    h_routine.stellarhsml( ctypes.c_int(N),
                          util.cast.vfloat(x),
                          util.cast.vfloat(y),
                          util.cast.vfloat(z),
                          ctypes.c_int(DesNgb),
                          ctypes.c_float(Hmax),
                          ctypes.byref(H_OUT) )
        
    ## now put the output arrays into a useful format
    h = np.ctypeslib.as_array(H_OUT);
    return h;


def get_gas_density_around_stars( x_gas, y_gas, z_gas, m_gas, x_star, y_star, z_star, DesNgb=32, Hmax=0. , **kwargs ):
    
    x_gas=util.cast.fcor(x_gas); y_gas=util.cast.fcor(y_gas); z_gas=util.cast.fcor(z_gas);
    m_gas=util.cast.fcor(m_gas);
    
    x_star=util.cast.fcor(x_star); y_star=util.cast.fcor(y_star); z_star=util.cast.fcor(z_star);
    N_star=util.cast.checklen(x_star);
    
    ok=(util.cast.ok_scan(x_gas) & util.cast.ok_scan(y_gas) & util.cast.ok_scan(z_gas) & util.cast.ok_scan(m_gas));
    x_gas=x_gas[ok]; y_gas=y_gas[ok]; z_gas=z_gas[ok]; m_gas=m_gas[ok];
    N_gas=util.cast.checklen(x_gas);
    
    
    if(Hmax==0.):
        dx=np.max(x_gas)-np.min(x_gas); dy=np.max(y_gas)-np.min(y_gas); dz=np.max(z_gas)-np.min(z_gas); ddx=np.max([dx,dy,dz]);
        Hmax=5.*ddx*(np.float(N_gas)**(-1./3.)); ## mean inter-particle spacing
    
    ## load the routine we need
    exec_call=util.dir.c_routines_dir()+'/StellarGasDensity/stargasdensity.so'
    h_routine=ctypes.cdll[exec_call];

    h_out_cast=ctypes.c_float*N_star; H_OUT=h_out_cast();
    ## main call to the hsml-finding routine
    h_routine.stellargasdensity( ctypes.c_int(N_gas), ctypes.c_int(N_star),
                          util.cast.vfloat(x_gas),  util.cast.vfloat(y_gas),    util.cast.vfloat(z_gas), util.cast.vfloat(m_gas),
                          util.cast.vfloat(x_star),  util.cast.vfloat(y_star),    util.cast.vfloat(z_star),
                          ctypes.c_int(DesNgb),
                          ctypes.c_float(Hmax),
                          ctypes.byref(H_OUT) )
        
                          ## now put the output arrays into a useful format
    h = np.ctypeslib.as_array(H_OUT);
    return h;


def get_gas_temperature_around_stars( x_gas, y_gas, z_gas, m_gas, rho_gas, u_gas, x_star, y_star, z_star, DesNgb=32, Hmax=0. , **kwargs ):

    x_gas=util.cast.fcor(x_gas); y_gas=util.cast.fcor(y_gas); z_gas=util.cast.fcor(z_gas);
    m_gas=util.cast.fcor(m_gas);
    rho_gas=util.cast.fcor(rho_gas);
    u_gas=util.cast.fcor(u_gas);

    x_star=util.cast.fcor(x_star); y_star=util.cast.fcor(y_star); z_star=util.cast.fcor(z_star);
    N_star=util.cast.checklen(x_star);

    ok=(util.cast.ok_scan(x_gas) & util.cast.ok_scan(y_gas) & util.cast.ok_scan(z_gas) & util.cast.ok_scan(m_gas));
    x_gas=x_gas[ok]; y_gas=y_gas[ok]; z_gas=z_gas[ok]; m_gas=m_gas[ok]; rho_gas=rho_gas[ok]; u_gas=u_gas[ok]
    N_gas=util.cast.checklen(x_gas);


    if(Hmax==0.):
        dx=np.max(x_gas)-np.min(x_gas); dy=np.max(y_gas)-np.min(y_gas); dz=np.max(z_gas)-np.min(z_gas); ddx=np.max([dx,dy,dz]);
        Hmax=5.*ddx*(np.float(N_gas)**(-1./3.)); ## mean inter-particle spacing

    ## load the routine we need
    exec_call=util.dir.c_routines_dir()+'/StellarGasTemperature/stargasavg.so'
    h_routine=ctypes.cdll[exec_call];

    h_out_cast=ctypes.c_float*N_star; H_OUT=h_out_cast();
    ## main call to the hsml-finding routine
    h_routine.stellargastemperature( ctypes.c_int(N_gas), ctypes.c_int(N_star),
                          util.cast.vfloat(x_gas),  util.cast.vfloat(y_gas),    util.cast.vfloat(z_gas), util.cast.vfloat(m_gas), util.cast.vfloat(rho_gas), util.cast.vfloat(u_gas),
                          util.cast.vfloat(x_star),  util.cast.vfloat(y_star),    util.cast.vfloat(z_star),
                          ctypes.c_int(DesNgb),
                          ctypes.c_float(Hmax),
                          ctypes.byref(H_OUT) )

                          ## now put the output arrays into a useful format
    h = np.ctypeslib.as_array(H_OUT);
    return h;

