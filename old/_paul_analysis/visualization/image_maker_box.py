import numpy as np
import matplotlib
#matplotlib.use('Agg') ## this calls matplotlib without a X-windows GUI
from subprocess import call
import matplotlib.colors
import matplotlib.cm
import matplotlib.pyplot as plt
import math
import utilities as util
import gadget
import gadget_lib.load_stellar_hsml as starhsml
import visualization.colors as viscolors
import visualization.get_attenuated_stellar_luminosities as getlum
import visualization.make_threeband_image as makethreepic
import visualization.contour_makepic as cmakepic
import visualization.raytrace_projection as rayproj
import h5py

def checklen(x):
    return len(np.array(x,ndmin=1));

def ok_scan(input,xmax=1.0e10,pos=0):
    if (pos==1):
        return (np.isnan(input)==False) & (np.abs(input)<=xmax) & (input > 0.);
    else:
        return (np.isnan(input)==False) & (np.abs(input)<=xmax);

def snap_ext(snum,four_char=0):
	ext='00'+str(snum);
	if (snum>=10): ext='0'+str(snum)
	if (snum>=100): ext=str(snum)
	if (four_char==1): ext='0'+ext
	if (snum>=1000): ext=str(snum)
	return ext;
	
def get_gas_xray_luminosity(ppp):
    brems = gadget.gas_xray_brems( \
        ppp['m'], ppp['u'], ppp['rho'], ppp['ne'], ppp['nh'] );
    return brems; ## can also add metal line cooling luminosities later

##
## handy function to load everything we need from snapshots, in particular to 
##   concatenate the various star types if we're using non-cosmological snapshots
##
def load_snapshot_brick(ptypes, snapdir, snapnum, h0=1, cosmological=0, \
        skip_bh=0, do_xray=0, four_char=0, use_rundir=0, full_gas_set=0 ):
    have=0; have_h_stars=0;
    ppp_head=gadget.readsnap(snapdir,snapnum,0,h0=h0,cosmological=cosmological,skip_bh=skip_bh,header_only=1)
    time=ppp_head['time']
    if (full_gas_set==1):
        ## here it will just return the entire gas 'brick'
        ppp=gadget.readsnap(snapdir,snapnum,0,h0=h0,cosmological=cosmological,skip_bh=skip_bh);
        m=ppp['m']; p=ppp['p']; x=p[:,0]; y=p[:,1]; z=p[:,2]; zm=ppp['z']; 
        if(checklen(zm.shape)>1): zm=zm[:,0]; lx=get_gas_xray_luminosity(ppp) / (3.9e33);
        sfr0=ppp['SFR']; sfr=1.*sfr0; p0=np.min(ppp['rho'][(sfr0 > 0.)]); 
        lo=(sfr0 <= 0.); sfr[lo]=0.25*(1./(1.+(ppp['rho'][lo]/p0)**(-0.5)))*np.min(sfr0[(sfr0>0.)]);
        return ppp['u'], ppp['rho'], ppp['h'], ppp['nh'], ppp['ne'], sfr, lx, zm, m, x, y, z, time;
        
    for ptype in ptypes:
        ppp=gadget.readsnap(snapdir,snapnum,ptype,h0=h0,cosmological=cosmological,skip_bh=skip_bh);
        if(ppp['k']==1):
            n=checklen(ppp['m']);
            if(n>1):
                m=ppp['m']; p=ppp['p']; x=p[:,0]; y=p[:,1]; z=p[:,2]; 
                if (ptype==0): 
                    cc=gadget.gas_temperature(ppp['u'],ppp['ne']);
                    if (do_xray==1): m=get_gas_xray_luminosity(ppp) / (3.9e33); ## bremstrahhlung (in solar)
                    h=ppp['h']; zm=ppp['z']; 
                    if(checklen(zm.shape)>1): zm=zm[:,0]
                if (ptype==4):
                    cc=gadget.get_stellar_ages(ppp,ppp_head,cosmological=cosmological);
                    zm=ppp['z']; 
                    if(checklen(zm.shape)>1): zm=zm[:,0]
                if (ptype==2): ## need to assign ages and metallicities
                    cc=np.random.rand(n)*(ppp_head['time']+4.0);
                    zm=(np.random.rand(n)*(1.0-0.1)+0.1) * 0.02;
                if (ptype==3): ## need to assign ages and metallicities
                    cc=np.random.rand(n)*(ppp_head['time']+12.0);
                    zm=(np.random.rand(n)*(0.3-0.03)+0.03) * 0.02;
                hstars_should_concat=0;
                if((ptype>0) & (have_h_stars==0)):
                    h=starhsml.load_allstars_hsml(snapdir,snapnum,cosmo=cosmological, \
                        use_rundir=use_rundir,four_char=four_char,use_h0=h0);
                    have_h_stars=1;
                    hstars_should_concat=1;
                if (have==1):
                    m_all=np.concatenate((m_all,m)); x_all=np.concatenate((x_all,x)); 
                    y_all=np.concatenate((y_all,y)); z_all=np.concatenate((z_all,z)); 
                    c_all=np.concatenate((c_all,cc)); zm_all=np.concatenate((zm_all,zm)); 
                    if(hstars_should_concat==1): h_all=np.concatenate((h_all,h)); 
                else:
                    m_all=m; x_all=x; y_all=y; z_all=z; zm_all=zm; c_all=cc; h_all=h; 
                    have=1;
    if (have==1):
        return m_all, x_all, y_all, z_all, c_all, h_all, zm_all, time;
    return 0,0,0,0,0,0,0,0;
    

def overlay_scale_label(xr,yr,figure_axis,c0='w'):
    ddx=0.5*(xr[1]-xr[0]);
    if (ddx>0.00002): dx=0.00001; dxl='0.01 pc';
    if (ddx>0.0002): dx=0.0001; dxl='0.1 pc';
    if (ddx>0.002): dx=0.001; dxl='1 pc';
    if (ddx>0.02): dx=0.01; dxl='10 pc';
    if (ddx>0.2): dx=0.1; dxl='100 pc';
    if (ddx>2.): dx=1.; dxl='1 kpc';
    if (ddx>20.): dx=10.; dxl='10 kpc';
    if (ddx>200.): dx=100.; dxl='100 kpc';
    if (ddx>2000.): dx=1000.; dxl='1 Mpc';
    if (ddx>20000.): dx=10000.; dxl='10 Mpc';
    
    xlen=(xr[1]-xr[0])
    ylen=(yr[1]-yr[0])
    xoff = (0.25+0.02)*ddx / xlen
    yoff = 0.025*(yr[1]-yr[0]) / ylen
    xr_new = np.array([xoff-dx/xlen*0.5,xoff+dx/xlen*0.5])
    yr_new = np.array([yoff,yoff])
    
    plt.text(xoff,1.55*yoff,dxl,color=c0,\
        horizontalalignment='center',verticalalignment='baseline',\
        transform=figure_axis.transAxes,fontsize=32)
    figure_axis.autoscale(False)
    figure_axis.plot(xr_new,yr_new,color=c0,linewidth=4.0,\
        transform=figure_axis.transAxes)

    
def overlay_time_label(time, figure_axis, cosmo=0, c0='w', n_sig_add=0, tunit_suffix='Gyr'):
    if(cosmo==0): 
        time_to_use = time ## absolute time in Gyr
        prefix = ''
        suffix = tunit
    else: ## 
        time_to_use = 1./time - 1. ## redshift
        prefix = 'z='
        suffix = ''
        
    n_sig = 2
    if time_to_use >= 10: n_sig += 0
    if time_to_use < 1.: n_sig += 1
    t_str = round_to_n( time_to_use, n_sig+n_sig_add )
    label_str = prefix+t_str+suffix
    xoff=0.03; yoff=1.-0.025;
    plt.text(xoff,yoff,label_str,color=c0,\
        horizontalalignment='left',verticalalignment='top',\
        transform=figure_axis.transAxes,fontsize=32)


def round_to_n(x, n):
    ''' Utility function used to round labels to significant figures for display purposes  from: http://mail.python.org/pipermail/tutor/2004-July/030324.html'''
    if n < 1:
        raise ValueError("number of significant digits must be >= 1")

    # show everything as floats (preference; can switch using code below to showing eN instead
    format = "%." +str(n-1) +"f"
    as_string=format % x
    return as_string

    # Use %e format to get the n most significant digits, as a string.
    format = "%." + str(n-1) + "e"
    as_string = format % x
    if as_string[-3:] in ['+00', '+01', '+02', '+03','-01', '-02', '-03']:
        #then number is 'small', show this as a float
        format = "%." +str(n-1) +"f"
        as_string=format % x
    return as_string


def coordinates_rotate(x_all, y_all, z_all, theta, phi, coordinates_cylindrical=0):
    ## set viewing angle for image
    ## project into plane defined by vectors perpendicular to the vector of this angle
    x = np.cos(phi)*x_all + np.sin(phi)*y_all + 0.*z_all
    y = -np.cos(theta)*np.sin(phi)*x_all + np.cos(theta)*np.cos(phi)*y_all + np.sin(theta)*z_all
    z =  np.sin(theta)*np.sin(phi)*x_all - np.sin(theta)*np.cos(phi)*y_all + np.cos(theta)*z_all
    if (coordinates_cylindrical==1):
        x=x/abs(x)*sqrt(x*x+z*z); y=y; z=z
    return x,y,z
  
 
#
# given a set of positions, project them to the coordinates they 'would be at'
#   with respect to a given camera (for ray-tracing purposes)
#
def coordinates_project_to_camera(x, y, z, \
        camera_pos=[0.,0.,0.], camera_dir=[0.,0.,1.], \
        screen_distance=1.0 ):
    camera_pos=np.array(camera_pos,dtype='d');
    camera_dir=np.array(camera_dir,dtype='d');
    ## recenter at camera location
    x-=camera_pos[0]; y-=camera_pos[1]; z-=camera_pos[2];
    ## line-of-sight distance along camera viewing angle
    znew = x*camera_dir[0] + y*camera_dir[1] + z*camera_dir[2]; 
    ## get perpendicular axes to determine position in 'x-y' space:
    cperp_x,cperp_y = util.return_perp_vectors(camera_dir);
    xnew = x*cperp_x[0] + y*cperp_x[1] + z*cperp_x[2];
    ynew = x*cperp_y[0] + y*cperp_y[1] + z*cperp_y[2];

    ## project to 'screen' at fixed distance along z-axis :: 
    znew[znew==0.]=-1.0e-10; ## just to prevent nans
    ## need correction factors for other quantities (like hsml, for example)
    r_corr = screen_distance/znew;
    xnew *= r_corr; ## put positions 'at the screen'
    ynew *= r_corr; ## put positions 'at the screen'
    znew = -znew; ## so values in the direction of camera are *negative* (for later sorting)

    return xnew, ynew, znew, r_corr

    

def test():
    smaster='/Users/phopkins/Documents/work/plots/zooms/'
    sdir='hires_bh'
    snum=441
    x00=10.0

    sdir='m12_mr'
    snum=429
    x00=20.0
    x00=50.0

    sdir='m11_hr'
    snum=440
    snum=316
    x00=20.

    sdir='m13_ics'
    snum=56
    x00=50.

    sdir='m10_dw'
    snum=440
    x00=10.0
    x00=100.0
    
    sdir='z5_m10_mr_wionbg_Dec16_2013'
    snum=50
    x00=4.
    
    sdir='m12v_3_Jan2014'
    snum=440
    x00=25.
    x00=50.

    out0,out1 = image_maker(sdir,snum,
        snapdir_master=smaster,outdir_master='./',pixels=720,
        do_with_colors=1,cosmo=1,xrange=np.array([-1.,1.])*x00,yrange=np.array([-1.,1.])*x00,
        #dynrange=3.0e2,maxden=3.8e7/1.0e10, show_gasstarxray='gas',include_lighting=1,
        ###dynrange=1.0e2,maxden=1.1e8/1.0e10, show_gasstarxray='star',include_lighting=0,
        #dynrange=0.6e2,maxden=1.1e8/1.0e10, show_gasstarxray='star',include_lighting=0,
        dynrange=0.5e2,maxden=0.8e8/1.0e10, show_gasstarxray='star',include_lighting=0,
        #theta=10.,phi=0.,threecolor=1, project_to_camera=1) #faceon
        theta=108.,phi=-60.,threecolor=1, project_to_camera=1) #edgeon
        #theta=17.,phi=-60.,threecolor=1, project_to_camera=0) # faceon
        #theta=107.,phi=-60.,threecolor=1, project_to_camera=0) # edgeon
    

def test2():
    smaster='/Users/phopkins/Documents/work/plots/zooms/'
    sdir='m12i_Jul30_2014_wRlimswk'
    snum=441
    x00=22.5

    out0,out1 = image_maker(sdir,snum,
        snapdir_master=smaster,outdir_master='./',pixels=720,
        do_with_colors=1,cosmo=1,xrange=np.array([-1.,1.])*x00,yrange=np.array([-1.,1.])*x00,
        #dynrange=5.e2,maxden=6.e-3, show_gasstarxray='gas',include_lighting=1,
        #theta=17.,phi=-60.,threecolor=1, project_to_camera=0) # faceon
        dynrange=1.e3,maxden=6.e-3, show_gasstarxray='gas',include_lighting=1,
        theta=107.,phi=-30.,threecolor=1, project_to_camera=0) # edgeon

    #map,pic,tmap,tpic = image_maker(sdir,snum,snapdir_master=smaster,outdir_master='./',\
    #    do_with_colors=1,cosmo=1,xrange=np.array([-1.,1.])*x00,yrange=np.array([-1.,1.])*x00,\
    #    dynrange=1.0e4,maxden_rescale=0.1, theta=90., phi=90., \
    #    show_gasstarxray='gas',log_temp_wt=0,set_temp_min=1.0e3,set_temp_max=1.e6,\
    #    project_to_camera=1 );
    #plt.imshow(tpic)
    #return map,pic,tmap,tpic

    if(2==0):
        out0,out1 = image_maker(sdir,snum,
            snapdir_master=smaster,outdir_master='./',
            do_with_colors=1,cosmo=1,xrange=np.array([-1.,1.])*x00,yrange=np.array([-1.,1.])*x00,
            #dynrange=3.0e2*1000.,maxden=3.8e7/1.0e10, show_gasstarxray='gas',
            #dynrange=3.0e2*3.,maxden=3.8e7*3.0/1.0e10, show_gasstarxray='star',include_lighting=0,
            #theta=90., phi=90.,threecolor=1, project_to_camera=1)
            #theta=60., phi=45.,threecolor=1, project_to_camera=1)
            #dynrange=3.0e2*10.,maxden=3.8e7/1.0e10, show_gasstarxray='gas',
            #dynrange=3.0e2*3.,maxden=3.8e7*3.0/1.0e10, show_gasstarxray='star',include_lighting=0,
            dynrange=3.0e2*10.,maxden=3.8e7*3.0/1.0e10, show_gasstarxray='star',include_lighting=0,
            theta=0., phi=0.,threecolor=0, project_to_camera=0)
        #plt.imshow(out2)
        return out0,out1;


    #map,pic = image_maker(sdir,snum,snapdir_master=smaster,outdir_master='./',\
    #    do_with_colors=0,cosmo=1,xrange=np.array([-1.,1.])*x00,yrange=np.array([-1.,1.])*x00, \
    #    show_gasstarxray='gas',dynrange=1.e5,maxden_rescale=100.,theta=75., 
    #     project_to_camera=1);
    #plt.imshow(pic)
    #return map,pic
    
    

def image_maker( sdir, snapnum, \
    snapdir_master='/n/scratch2/hernquist_lab/phopkins/sbw_tests/', 
    outdir_master='/n/scratch2/hernquist_lab/phopkins/images/', 
    theta=0., phi=0., dynrange=1.0e5, maxden_rescale=1., maxden=0.,
    show_gasstarxray = 'gas', #	or 'star' or 'xray'
    add_gas_layer_to_image='', set_added_layer_alpha=0.3, set_added_layer_ctable='heat_purple', 
    add_extension='', show_time_label=1, show_scale_label=1, 
    filename_set_manually='',
    do_with_colors=1, log_temp_wt=1, include_lighting=1, 
    set_percent_maxden=0, set_percent_minden=0, 
    center_on_com=0, center_on_bh=0, use_h0=1, cosmo=0, 
    center=[0., 0., 0.], coordinates_cylindrical=0, 
    pixels=720,xrange=[-1.,1.],yrange=0,zrange=0,set_temp_max=0,set_temp_min=0,
    threecolor=1, nasa_colors=1, sdss_colors=0, use_old_extinction_routine=0, 
    dust_to_gas_ratio_rescale=1.0, 
    project_to_camera=1, camera_opening_angle=45.0,
    center_is_camera_position=0, camera_direction=[0.,0.,-1.], 
    gas_map_temperature_cuts=[1.0e4, 1.0e6], 
    input_data_is_sent_directly=0, \
    m_all=0,x_all=0,y_all=0,z_all=0,c_all=0,h_all=0,zm_all=0,\
    gas_u=0,gas_rho=0,gas_hsml=0,gas_numh=0,gas_nume=0,gas_metallicity=0,\
    gas_mass=0,gas_x=0,gas_y=0,gas_z=0,time=0,
    invert_colors=0,spin_movie=0,scattered_fraction=0.01,
    min_stellar_age=0,h_rescale_factor=1.,h_max=0, h_e_factor=1.e10,
    max_age_for_h_max = 1.e10 ):
    

	## define some of the variables to be used
    ss=snap_ext(snapnum,four_char=1);
    tt=snap_ext(np.around(theta).astype(int));
    theta *= math.pi/180.; phi *= math.pi/180.; # to radians
    nameroot = sdir+'_s'+ss+'_t'+tt;
    outputdir = outdir_master+sdir;
    call(["mkdir",outputdir]);
    outputdir+='/'; snapdir=snapdir_master+sdir;
    suff='_'+show_gasstarxray;
    if (threecolor==1):
        if (sdss_colors==1): suff+='_S3c'
        if (nasa_colors==1): suff+='_N3c'
    suff+=add_extension;
    fname_base=outputdir+nameroot+suff
    do_xray=0; do_stars=0; 
    if((show_gasstarxray=='xr') or (show_gasstarxray=='xray')): do_xray=1; do_with_colors=0;
    if((show_gasstarxray=='star') or (show_gasstarxray=='stars') or (show_gasstarxray=='st')): do_stars=1;
    
    
    ## check whether the data needs to be pulled up, or if its being passed by the calling routine
    if(input_data_is_sent_directly==0):
        ## read in snapshot data and do centering 
        ptypes=[2,3,4]; ## stars in non-cosmological snapshot
        if (cosmo==1): ptypes=[4];
        if (do_stars==0): ptypes=[0];
        m_all, x_all, y_all, z_all, c_all, h_all, zm_all, time = load_snapshot_brick(ptypes, \
            snapdir, snapnum, h0=use_h0, cosmological=cosmo, \
            skip_bh=cosmo, do_xray=do_xray, four_char=0, use_rundir=1, full_gas_set=0);
        h_all *= 1.25;

        if ((do_stars==1) & (threecolor==1)): ## will need gas info to process attenuation
            gas_u, gas_rho, gas_hsml, gas_numh, gas_nume, gas_sfr, gas_lxray, gas_metallicity, gas_mass, gas_x, gas_y, gas_z, time = \
             load_snapshot_brick([0], snapdir, snapnum, h0=use_h0, cosmological=cosmo, full_gas_set=1);
            ## don't allow hot gas to have dust
            gas_temp = gadget.gas_temperature(gas_u,gas_nume); gas_metallicity[gas_temp > 1.0e6] = 0.0; 

    
        if (center[0] == 0.):
            if (cosmo==1): 
                center=gadget.calculate_zoom_center(snapdir,snapnum);
            else:
                if (center_on_bh==1):
                    pbh=gadget.readsnap(snapdir,snapnum,ptype,h0=use_h0,cosmological=0,skip_bh=0);
                    pos=pbh['p']; center=[pos[0,0],pos[0,1],pos[0,2]];
                if (center_on_com==1):
                    center=[np.median(x_all),np.median(y_all),np.median(z_all)];
        x_all-=center[0]; y_all-=center[1]; z_all-=center[2];
        print 'center at ',center


    ## rotate and re-project coordinate frame 
    if (center_is_camera_position==0):
        x,y,z=coordinates_rotate(x_all,y_all,z_all,theta,phi,coordinates_cylindrical=coordinates_cylindrical);
        if ((do_stars==1) & (threecolor==1)):
            gas_x-=center[0]; gas_y-=center[1]; gas_z-=center[2];
            gx,gy,gz=coordinates_rotate(gas_x,gas_y,gas_z,theta,phi,coordinates_cylindrical=coordinates_cylindrical);
    else:
        x=x_all; y=y_all; z=z_all; 
        if ((do_stars==1) & (threecolor==1)): gx=gas_x; gy=gas_y; gz=gas_z;


    ## set dynamic ranges of image
    temp_max=1.0e6; temp_min=1.0e3; ## max/min gas temperature
    if (do_stars==1): temp_max=50.; temp_min=0.01; ## max/min stellar age
    if (set_temp_max != 0): temp_max=set_temp_max
    if (set_temp_min != 0): temp_min=set_temp_min	
    xr=xrange; yr=xr; zr=0;
    if (checklen(yrange)>1): yr=yrange;
    if (checklen(zrange)>1): zr=zrange;
    scale=0.5*np.max([xr[1]-xr[0],yr[1]-yr[0]]); zr=np.array([-1.,1.])*scale;
    if(maxden==0.):
        maxden = 0.1 * 1.e11/1.0e10 * (2.0/scale)**(0.3) * maxden_rescale;
    if ((threecolor==1) & (do_stars)):
        if(dynrange==1.0e5): # reset if it's the default value 
            dynrange=1.0e2;
            if (nasa_colors==1): dynrange=1.0e4;
            if (nasa_colors==1): maxden *= 30.;


    ## now check if we're projecting to a camera (instead of a fixed-plane)
    xlen=0.5*(xr[1]-xr[0]); ylen=0.5*(yr[1]-yr[0]); 
    xr_0=xr; yr_0=yr;
    if (project_to_camera==1):  
        ## determine the angular opening and camera distance given physical x/y range:
        ## note: camera_opening_angle is the half-angle
        xr=np.array([-1.,1.])*np.tan(camera_opening_angle*math.pi/180.); yr=xr*ylen/xlen; ## n- degree opening angle
        ## use this to position camera
        c_dist = xlen/xr[1]; 
        camera_direction /= np.sqrt(camera_direction[0]**2.+camera_direction[1]**2.+camera_direction[2]**2.);
        camera_pos = -c_dist * camera_direction;
        ## and determine z-depth of image
        cpmax=np.sqrt(camera_pos[0]**2.+camera_pos[1]**2.+camera_pos[2]**2.); zrm=np.sqrt(zr[0]**2.+zr[1]**2.);
        ## clip z_range for speed : 
        zrm=5.0*zrm; 
        #zrm=2.5*zrm
        zrm=np.sqrt(zrm*zrm+cpmax*cpmax); zr=[-zrm,0.]
        ## correct if the center given is the camera position:
        if (center_is_camera_position==1): camera_pos=[0.,0.,0.];

        ## now do the actual projection into camera coordinates:        
        x,y,z,rc=coordinates_project_to_camera(x,y,z,camera_pos=camera_pos,camera_dir=camera_direction);
        print "c_dist = ",c_dist
        h_all*=rc; ## correct for size re-scaling
        # XXX temp for VoT
        #if not ((threecolor==0) & (do_stars==0)):
        #    m_all*=1./(z*z + h_all*h_all + (0.25*c_dist)**2.); ## makes into 'fluxes' (dimmer when further away)

	m_all *= 1./(z*z + h_all*h_all + (0.25*c_dist)**2.);
	#m_all *= c_all
	h_all *= 1.15
	#h_all *= 1.0 + (h_all/0.05)**2.0 # good, maybe -too- smoothed?
        h_all *= 1.0 + (h_all/0.2)**1.0 # (h_all/0.1)**1.5 # good

        ## clip particles too close to camera to prevent artifacts
        z_clip = -c_dist/10.0
        print "z_clip = ",z_clip
        m_all[z >= z_clip] = 0.0
        print "Max z: ",np.max(z[m_all > 0]),"; min z: ",np.min(z[m_all > 0]),"; median z: ",\
            np.median(z[m_all > 0])
        print "z = ",z[m_all > 0]
        ## also need to correct gas properties if using it for attenuation
        if ((threecolor==1) & (do_stars==1)):
            gx,gy,gz,grc=coordinates_project_to_camera(gx,gy,gz,camera_pos=camera_pos,camera_dir=camera_direction);
            gas_hsml*=grc; gas_mass*=1./(gz*gz + gas_hsml*gas_hsml + (0.25*c_dist)**2.); ## this conserves surface density, so attenuation units are fine
            gas_mass[gz >= z_clip] = 0.0 ## clipping to prevent artifacts

    #plotting_setup_junk
    plt.close('all')
    format = '.png' ## '.ps','.pdf','.png', '.eps' work well
    axis_ratio = ylen/xlen
    #fig=plt.figure(frameon=False,figsize=(1.,1.*axis_ratio),dpi=pixels)
    fig=plt.figure(frameon=False,dpi=pixels)
    fig.set_size_inches(1,1)
    ax_fig=plt.Axes(fig,[0.,0.,1.,1.*axis_ratio])
    ax_fig.set_axis_off()
    fig.add_axes(ax_fig)

	  
	## trim particles to those inside/near the plotting region
    sidebuf=10.; 
    if (cosmo==1): sidebuf=1.15
    x=np.array(x,dtype='f'); y=np.array(y,dtype='f'); z=np.array(z,dtype='f'); 
    m_all=np.array(m_all,dtype='f'); h_all=np.array(h_all,dtype='f'); 
    c_all=np.array(c_all,dtype='f'); 
    ok =    ok_scan(x,xmax=np.max(np.fabs(xr))*sidebuf) & ok_scan(y,xmax=np.max(np.fabs(yr))*sidebuf) & \
            ok_scan(z,xmax=np.max(np.fabs(zr))*sidebuf) & ok_scan(m_all,pos=1,xmax=1.0e40) & \
            ok_scan(h_all,pos=1) & ok_scan(c_all,pos=1,xmax=1.0e40)
    weights = m_all; color_weights = c_all; ## recall, gadget masses in units of 1.0d10


    TEMPORARY_FOR_THORSTEN_REVIEW_ONLY = 0
    if(TEMPORARY_FOR_THORSTEN_REVIEW_ONLY==1):
        #xxx=1.*x; x=1.*y; y=1.*xxx;
        if ((threecolor==1) & (do_stars==1)):
            #gxxx=1.*gx; gx=1.*gy; gy=1.*gxxx;
            h_all *= 0.8
            #c_all *= (c_all / 1.0);
            gas_metallicity *= 1.e-5;
            c_all *= (c_all / 1.0)**0.5; gas_metallicity *= 0.67;
            #c_all *= (c_all / 1.0)**0.25; gas_metallicity *= 0.67;
            #c_all *= (c_all / 1.0)**0.5; gas_metallicity *= 0.67;
            #c_all *= (c_all / 1.0)**0.25; 
            use_old_extinction_routine = 0;

    TEMPORARY_FOR_SPIN_MOVIE_ONLY = 0.
    if(TEMPORARY_FOR_SPIN_MOVIE_ONLY==1):
        if ((threecolor==1) & (do_stars==1)):
            h_all *= ((scale/30.)**0.25) * (1.0 + (h_all / 0.1)**(0.5)) / (1. + (h_all / 1.0)**(0.5))
            #c_all *= (c_all / 1.0)**0.25; 
            c_all *= (c_all / 1.0)**0.125; 
        else:
            h_all *= (((30.*0.+1.*scale)/30.)**0.25) * (0.8 + (h_all / 0.1)**(0.5)) / (1. + (h_all / 1.0)**(0.5))

    ## adjust softenings if appropriate keywords are set
    h_all *= h_rescale_factor
    if (2==0) & (h_max > 0):
        #h_all = np.minimum(h_all,h_max*np.exp(-(z-z_clip)/(c_dist*h_e_factor)))
        #max_age_for_h_max = 0.5 # prevent foreground particles <500 Myr old from being blurry
        h_all[c_all < max_age_for_h_max] = np.minimum(h_all[c_all < max_age_for_h_max],
            h_max*(1+(-(z[c_all < max_age_for_h_max]-z_clip)/(c_dist*h_e_factor))**5))
            # linear
            #h_max*(-(z[c_all < max_age_for_h_max]-z_clip)/(c_dist*h_e_factor)))
        #h_all = np.minimum(h_all,h_max*(1+(-(z-z_clip)/(c_dist*h_e_factor))**3))
        print "c_dist*h_e_factor = ",c_dist*h_e_factor
        print "-(z-z_clip)/(c_dist*h_e_factor): median = ",\
            np.median(-(z[m_all > 0]-z_clip)/(c_dist*h_e_factor)),"; max = ",\
            np.max(-(z[m_all > 0]-z_clip)/(c_dist*h_e_factor)),"; min = ",\
            np.min(-(z[m_all > 0]-z_clip)/(c_dist*h_e_factor))
        print "1+((z-z_clip)/(c_dist*h_e_factor))**2: median = ",\
            np.median(1+((z[m_all > 0]-z_clip)/(c_dist*h_e_factor))**2),"; max = ",\
            np.max(1+((z[m_all > 0]-z_clip)/(c_dist*h_e_factor))**2),"; min = ",\
            np.min(1+((z[m_all > 0]-z_clip)/(c_dist*h_e_factor))**2)
        print "Mean h_all = ",np.mean(h_all[m_all > 0])
        print "Median h_all = ",np.median(h_all[m_all > 0])
        print "Max h_all = ",np.max(h_all[m_all > 0])
        print "Min h_all = ",np.min(h_all[m_all > 0])

    ## alright, now ready for the main plot construction:
    if (threecolor==1):
        if (do_stars==1):
            ## making a mock 3-color composite image here ::
            
            ## first grab the appropriate, attenuated stellar luminosities
            BAND_IDS=[9,10,11]  ## ugr composite
            ## for now not including BHs, so initialize some appropriate dummies:
            bh_pos=[0,0,0]; bh_luminosity=0.; include_bh=0;
            ## and some limits for the sake of integration:
            SizeMax=20.*np.max([np.fabs(xr[1]-xr[0]),np.fabs(yr[1]-yr[0])]);
            SizeMin=np.min([np.fabs(xr[1]-xr[0]),np.fabs(yr[1]-yr[0])])/250.;
 
            ## now set up the main variables:
            star_pos=np.zeros((3,checklen(x[ok]))); 
            star_pos[0,:]=x[ok]; star_pos[1,:]=y[ok]; star_pos[2,:]=z[ok];
            stellar_age=np.maximum(c_all,min_stellar_age);
            stellar_metallicity=zm_all; stellar_mass=m_all;
            ## gas will get clipped in attenuation routine, don't need to worry about it here
            gas_pos=np.zeros((3,checklen(gx))); 
            gas_pos[0,:]=gx; gas_pos[1,:]=gy; gas_pos[2,:]=gz;

            if ((add_gas_layer_to_image=='')==False):
                gas_wt = 0.*gx
                maxden_for_layer=1.*maxden; dynrange_for_layer=1.*dynrange
                set_percent_maxden_layer=0.; set_percent_minden_layer=0.;
                gas_hsml_for_extra_layer=gas_hsml
                if (add_gas_layer_to_image=='Halpha'): 
                    gas_wt = gas_mass * gas_rho*gas_nume*(1.-gas_numh)*(gadget.gas_temperature(gas_u,gas_nume)**(-0.75))/(gadget.gas_mu(gas_nume)**2.)
                    maxden_for_layer *= 3.e4;
                    #dynrange_for_layer *= 1.5;
                    dynrange_for_layer *= 0.7;
                    dynrange_for_layer *= 2.0;
                    #gas_hsml_for_extra_layer *= 1.1;
                if (add_gas_layer_to_image=='CO'): 
                    gas_tmp = gadget.gas_temperature(gas_u,gas_nume)
                    n_tmp = gas_rho * 176.2; # assuming mean molec weight of 2.3 for dense gas
                    gas_wt = gas_mass * np.exp(-(gas_tmp/8000. + 10./n_tmp));
                    maxden_for_layer *= 3.e4;
                    dynrange_for_layer *= 0.7e4;
                if (add_gas_layer_to_image=='SFR'): 
                    gas_wt = gas_sfr
                    set_percent_maxden_layer=0.9999;
                    set_percent_minden_layer=0.01;
                if (add_gas_layer_to_image=='Xray'): 
                    gas_wt = gas_lxray
                    maxden_for_layer *= 0.6e3;
                    dynrange_for_layer *= 0.1;
                    gas_hsml_for_extra_layer *= 1.2;
                if (add_gas_layer_to_image=='Zmetal'): 
                    gas_wt = gas_mass*gas_metallicity
                    set_percent_maxden_layer=0.9999;
                    set_percent_minden_layer=0.01;
                gas_wt /= np.sum(gas_wt)
                                    
                massmap_gas_extra_layer,image_singledepth_extra_layer = \
                cmakepic.simple_makepic(gx,gy,weights=gas_wt,hsml=gas_hsml_for_extra_layer,\
                    xrange=xr,yrange=yr,
                    set_dynrng=dynrange_for_layer,set_maxden=maxden_for_layer,
                    set_percent_maxden=set_percent_maxden_layer,set_percent_minden=set_percent_minden_layer, 
                    color_temperature=0,pixels=pixels,invert_colorscale=1-invert_colors);


            if (use_old_extinction_routine == 0):
                ##
                ## this is the newer 'single call' attenuation and ray-tracing package: 
                ##  slightly more accurate in diffuse regions, more stable behavior
                ##

                if (x[ok].size <= 3):
                    ## we got a bad return from the previous routine, initialize blank arrays
                    out_gas=out_u=out_g=out_r=np.zeros((pixels,pixels))
                    image24=massmap=np.zeros((pixels,pixels,3))
                else:
                    ## actually call the ray-trace:
                    out_gas,out_u,out_g,out_r = rayproj.stellar_raytrace( BAND_IDS, \
                        x[ok], y[ok], z[ok], \
                        stellar_mass[ok], stellar_age[ok], stellar_metallicity[ok], h_all[ok], \
                        gx, gy, gz, gas_mass, \
                        gas_metallicity*dust_to_gas_ratio_rescale, gas_hsml, \
                        xrange=xr, yrange=yr, zrange=zr, pixels=pixels, \
                        ADD_BASE_METALLICITY=0.1*0.02, ADD_BASE_AGE=0.0003, 
                        #ADD_BASE_METALLICITY=0.001*0.02, ADD_BASE_AGE=0.0003, 
                        IMF_SALPETER=0, IMF_CHABRIER=1 );
                
                    if(np.array(out_gas).size<=1):
                        ## we got a bad return from the previous routine, initialize blank arrays
                        out_gas=out_u=out_g=out_r=np.zeros((pixels,pixels))
                        image24=massmap=np.zeros((pixels,pixels,3))
                    else:
                        ## make the resulting maps into an image
                        image24, massmap = \
                            makethreepic.make_threeband_image_process_bandmaps( out_r,out_g,out_u, \
                            maxden=maxden,dynrange=dynrange,pixels=pixels, \
                            color_scheme_nasa=nasa_colors,color_scheme_sdss=sdss_colors );  
                
            else: ## use_old_extinction_routine==1; can be useful for certain types of images
                ## (for example, this can allow for a scattered fraction, the other does not currently)

                ## call routine to get the post-attenuation band-specific luminosities
                lum_noatten, lum_losNH = \
                getlum.get_attenuated_stellar_luminosities( BAND_IDS, star_pos, gas_pos, bh_pos, \
                    stellar_age[ok], stellar_metallicity[ok], stellar_mass[ok], \
                    gas_u, gas_rho, gas_hsml, gas_numh, gas_nume, \
                    gas_metallicity*dust_to_gas_ratio_rescale, gas_mass, \
                    bh_luminosity, \
                    xrange=xr, yrange=yr, zrange=zr, \
                    INCLUDE_BH=include_bh, SKIP_ATTENUATION=0, 
                    #ADD_BASE_METALLICITY=1.0, ADD_BASE_AGE=0., 
                    ADD_BASE_METALLICITY=1.0e-2, ADD_BASE_AGE=0., 
                    IMF_SALPETER=0, IMF_CHABRIER=1, \
                    MIN_CELL_SIZE=SizeMin, OUTER_RANGE_OF_INT=SizeMax, \
                    SCATTERED_FRACTION=scattered_fraction, \
                    REDDENING_SMC=1, REDDENING_LMC=0, REDDENING_MW=0, \
                    AGN_MARCONI=0, AGN_HRH=1, AGN_RICHARDS=0, AGN_SDSS=0 );
                
                ## now pass these into the 3-band image routine
                image24, massmap = \
                makethreepic.make_threeband_image(x[ok],y[ok],lum_losNH,hsml=h_all[ok],\
                    xrange=xr,yrange=yr,maxden=maxden,dynrange=dynrange,pixels=pixels, \
                    color_scheme_nasa=nasa_colors,color_scheme_sdss=sdss_colors );           

        else: 
            ## 
            ## threecolor==1, but do_stars==0, so doing a multi-pass gas image
            ## 
            ##  -- many experiments here: doing gas isosurfaces with broad kernels
            ##       and overlaying a custom set of color tables after the fact seems
            ##       best. adjust opacity (kappa_units) and kernel width as needed. 
            ##       also, can set 'dynrange=0' to automatically produce a given 
            ##       fraction of saturated/unsaturated pixels
            ##
            #gas_map_temperature_cuts=np.array([300., 2.0e4, 3.0e5 ])
            #kernel_widths=np.array([0.8,0.3,0.6])
            #kernel_widths=np.array([0.7,0.3,0.7])

            # layer below is for zoom_dw large-scale GAS sims ONLY
            gas_map_temperature_cuts=np.array([300., 1.0e4, 1.0e5 ])
            kernel_widths=np.array([0.5,0.25,0.6])


            gas_map_temperature_cuts=np.array([1.0e4, 1.0e5, 1.0e6 ])
            kernel_widths=np.array([0.3,0.3,0.5])
	    #kappas = 2.0885*np.array([1.1,2.0,5.0]) 
            kappas = 2.0885*np.array([0.5,2.0,10.0])
	    ##kappas = 2.0885*np.array([0.5,2.0,1.0])

            out_gas,out_u,out_g,out_r = rayproj.gas_raytrace_temperature( \
                gas_map_temperature_cuts, \
                x[ok], y[ok], z[ok], color_weights[ok], weights[ok], h_all[ok], \
                xrange=xr, yrange=yr, zrange=zr, pixels=pixels, \
                isosurfaces = 1, kernel_width=kernel_widths, \
                add_temperature_weights = 0 , KAPPA_UNITS = kappas);
                #add_temperature_weights = 0 , KAPPA_UNITS = 2.0885*np.array([4.1,2.0,2.0]));
            
	    #out_u *= 1.; out_g *= 6.; out_r *= 20.; # weight correction to balance current B1 color scheme
            out_u *= 1.; out_g *= 6.; out_r *= 50.; # weight correction to balance current B1 color scheme
	    dr0 = 30.;

            image24, massmap = \
                makethreepic.make_threeband_image_process_bandmaps( out_r,out_g,out_u,pixels=pixels, \
		maxden=0,dynrange=dr0,color_scheme_nasa=1,color_scheme_sdss=0);
	
	    #massmap[:,:,0] *= 20.0;
	    #massmap[:,:,1] *=  6.0;
            #massmap[:,:,2] *=  1.0;

            #image24 = makethreepic.layer_band_images(image24, massmap);                

    else:   ## threecolor==0    ( not a three-color image )
    
        if (do_with_colors==1):
            ##
            ## here we have a simple two-pass image where surface brightness is 
            ##   luminosity and a 'temperature' is color-encoded
            ##
            if (log_temp_wt==1): 
                color_weights=np.log10(color_weights); temp_min=np.log10(temp_min); temp_max=np.log10(temp_max);

            massmap_singlelayer,pic_singlelayer,massmap,image24 = \
            cmakepic.simple_makepic(x[ok],y[ok],weights=weights[ok],hsml=h_all[ok],\
                xrange=xr,yrange=yr,
                set_dynrng=dynrange,set_maxden=maxden,
                set_percent_maxden=set_percent_maxden,set_percent_minden=set_percent_minden, 
                color_temperature=1,temp_weights=color_weights[ok],
                pixels=pixels,invert_colorscale=invert_colors,
                set_temp_max=temp_max,set_temp_min=temp_min);

        else: ## single color-scale image
            ##
            ## and here, a single-color scale image for projected density of the quantity
            ##
            massmap,image_singledepth = \
            cmakepic.simple_makepic(x[ok],y[ok],weights=weights[ok],hsml=h_all[ok],\
                xrange=xr,yrange=yr,
                set_dynrng=dynrange,set_maxden=maxden,
                set_percent_maxden=set_percent_maxden,set_percent_minden=set_percent_minden, 
                color_temperature=0,pixels=pixels,invert_colorscale=1-invert_colors);

            ## XXX temporary for VoT; make mass-weighted temperature map
            print "weights -- should be mass: max = ",np.max(weights),"   min = ",\
                np.min(weights),"    sum = ",weights.sum()
            print "color_weights -- should be T: max = ",np.max(color_weights),\
                "   min = ",np.min(color_weights)
            #if (log_temp_wt==1):
            #    print "log_temp_wt set"
            #    color_weights=np.log10(color_weights)

            tempmap_temp,tempimage_singledepth = \
            cmakepic.simple_makepic(x[ok],y[ok],weights=weights[ok]*color_weights[ok],
                hsml=h_all[ok],\
                xrange=xr,yrange=yr,
                set_dynrng=dynrange,set_maxden=maxden,
                set_percent_maxden=set_percent_maxden,set_percent_minden=set_percent_minden,
                color_temperature=0,pixels=pixels,invert_colorscale=1-invert_colors);

            tempmap = tempmap_temp/massmap
            print "TempMap: max = ", np.max(tempmap),"   min = ",np.min(tempmap),\
                "    mean = ", np.mean(tempmap), "    median = ",np.median(tempmap)
            #Tsort=np.sort(np.ndarray.flatten(tempmap))
            #ta=Tsort[0.99*float(checklen(tempmap)-1)]
            #ti=Tsort[0.01*float(checklen(tempmap)-1)]
            ta=2.0e6
            ti=8.0e3
            print "Clipping at   ta= ", ta, " ti= ", ti
            tempmap[tempmap < ti]=ti; tempmap[tempmap > ta]=ta

            cols=255. # number of colors
            # XXX temp for VoT; IDL script has cols-2
            #pic=(np.log(tempmap/ti)/np.log(ta/ti)) * (cols-3.) + 2.
            pic=(np.log(tempmap/ti)/np.log(ta/ti)) * (cols-2.) + 2.
            #pic[pic > cols]=cols; pic[pic < 1]=1.
            #backgrd = np.where((pic<=2) | (np.isnan(pic)))
            #invert_colorscale=1-invert_colors
            #if (invert_colorscale==0):
            #pic=256-pic;
            #pic[backgrd] = 1; # white
            #else:
            #    pic[backgrd] = 0; # black

            ## convert to actual image using a colortable:
            my_cmap=matplotlib.cm.get_cmap('hot');
            # XXX temp for VoT
            #image24 = my_cmap(image_singledepth/255.);
            image24 = my_cmap(pic/255.)

    ##
    ## whichever sub-routine you went to, you should now have a massmap (or set of them) 
    ##   and processed rgb 'image24' (which should always be re-makable from the massmaps)
    ##
    ## optionally can have great fun with some lighting effects to give depth:
    ##   (looks great, with some tuning)
    ##
    if (include_lighting==1):
        #light = matplotlib.colors.LightSource(azdeg=0,altdeg=65)
        light = viscolors.CustomLightSource(azdeg=0,altdeg=65)
        if (len(massmap.shape)>2):
            ## do some clipping to regulate the lighting:
            elevation = massmap.sum(axis=2)
            minden = maxden / dynrange
            elevation = (elevation - minden) / (maxden - minden)
            elevation[elevation < 0.] = 0.
            elevation[elevation > 1.] = 1.
            elevation *= maxden
            grad_max = maxden / 5.
            grad_max = maxden / 6.
            #image24_lit = light.shade_rgb(image24, massmap.sum(axis=2))
            image24_lit = light.shade_rgb(image24, elevation, vmin=-grad_max, vmax=grad_max)
        else:
            image24_lit = light.shade(image24, massmap)
        image24 = image24_lit
        
    plt.imshow(image24,origin='lower',interpolation='bicubic',aspect='auto');		# ptorrey aspect normal->auto

    if ((add_gas_layer_to_image=='')==False):
        viscolors.load_my_custom_color_tables();
        plt.imshow(image_singledepth_extra_layer,origin='lower',interpolation='bicubic',aspect='auto',	# ptorrey aspect normal->auto
            cmap=set_added_layer_ctable,alpha=set_added_layer_alpha)

    ## slap on the figure labels/scale bars
    if(show_scale_label==1): overlay_scale_label(xr_0,yr_0,ax_fig,c0='w')
    if(show_time_label==1): overlay_time_label(time,\
        ax_fig,cosmo=cosmo,c0='w',n_sig_add=0,tunit_suffix='Gyr')    
    
    filename=fname_base
    if(filename_set_manually!=''): filename=filename_set_manually
    plt.savefig(filename+format,dpi=pixels)#,bbox_inches='tight',pad_inches=0)
    plt.close('all')

    ## save them:
    outfi = h5py.File(filename+'.dat','w')
    dset_im = outfi.create_dataset('image24',data=image24)
    dset_mm = outfi.create_dataset('massmap',data=massmap)
    dset_xr = outfi.create_dataset('xrange',data=xr_0)
    dset_yr = outfi.create_dataset('yrange',data=yr_0)
    dset_time = outfi.create_dataset('time',data=time)
    dset_cosmo = outfi.create_dataset('cosmological',data=cosmo)
    if ((add_gas_layer_to_image=='')==False):
        dset_im2 = outfi.create_dataset('image24_extralayer',data=image_singledepth_extra_layer)
    else:
        dset_im2 = outfi.create_dataset('image24_extralayer',data=np.array([-1]))
    outfi.close()
    ## (to read) :::
    #infiname = fname_base+'.dat'
    #infi=h5py.File(infiname,'r')
    #image24 = np.array(infi["image24"])
    #massmap = np.array(infi["massmap"])
    #infi.close()
    
    ## and return the arrays to the parent routine!
    return image24, massmap;


