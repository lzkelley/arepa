import visualization.image_maker as im_make
import numpy as np
import matplotlib
#matplotlib.use('Agg') ## this calls matplotlib without a X-windows GUI
import matplotlib.pylab as pylab
import matplotlib.pyplot as plot
import h5py
import math
import visualization.contour_makepic as cmakepic
import visualization.colors as viscolors


def output():
    construct_spin_movie_single_sim(nf_base=10, key='star', add_key='nodust', pixels=800,frame_min=27.,frame_max=27.)
    #construct_spin_movie_single_sim(nf_base=10, key='star', add_key='nodust', pixels=800,frame_min=27.,frame_max=37.)
    #construct_spin_movie_single_sim(nf_base=100,key='star',add_key='Zmetal',pixels=800,
    #construct_spin_movie_single_sim(nf_base=100,key='star',add_key='Halpha',pixels=800,frame_min=250,frame_max=250)
    #construct_spin_movie_single_sim(nf_base=100,key='star',
    #    add_key='Halpha',pixels=800,frame_min=60,frame_max=60)


def construct_spin_movie_single_sim( \
    smaster = '/Users/phopkins/Documents/work/plots/zooms/', \
    sdir = 'm12v_3_Jan2014', \
    snum = 440, \
    omaster = '/Users/phopkins/Documents/work/plots/zooms/m12v_3_Jan2014/spin_movie_frames/', \
    nf_base = 100, \
    key = 'star', \
    add_key = '', \
    x00_in = 2.0, \
    x00_med = 30.0, \
    frac_to_rotate_initial_spin = 0.76, \
    pixels = 800., \
    theta_initial = 11., \
    theta_median = 109., \
    phi_initial = 0., \
    phi_median = -60., \
    frame_min = 0, \
    frame_max = -1,\
    scattered_fraction = 0.05,\
    h_rescale_factor=1.,\
    h_max=0):
    
    nf1 = nf_base
    nf2 = nf_base
    nf2x = np.int(np.round(nf_base/1.5))
    nf3 = nf2
    frac1 = frac_to_rotate_initial_spin
    t0=theta_initial; t1=theta_median; p0=phi_initial; p1=phi_median;
    nf2 = np.int(np.round(nf_base*1.5))
    nf2x = np.int(np.round(nf_base*1.0))
    nf3 = nf_base
    
    n = nf1
    dt = frac1*(t1-t0)/n; dt0=dt;
    dp = frac1*(p1-p0)/n; dp0=dp;
    dlnx = 0.0; dlnx0=dlnx;
    theta_grid = t0 + dt*np.arange(0,n)
    phi_grid = p0 + dp*np.arange(0,n)
    x00_grid = x00_med * np.exp(dlnx*np.arange(0,n))

    n = nf2
    dt = 1.*(1.-frac1)*(t1-t0)/n;
    dp = (1.-frac1)*(p1-p0)/n;
    theta_grid,dt0 = tcompile(theta_grid,dt,dt0,n)
    phi_grid,dp0 = tcompile(phi_grid,dp,dp0,n)
    x00_grid,dlnx0 = xcompile(x00_grid,x00_in,dlnx0,n)

    n = nf2x
    dt = 1.*(1.-frac1)*(t1-t0)/n
    dp = 2.*(1.-frac1)*(p1-p0)/n
    theta_grid,dt0 = tcompile(theta_grid,dt,dt0,n)
    phi_grid,dp0 = tcompile(phi_grid,dp,dp0,n)
    x00_grid,dlnx0 = xcompile(x00_grid,0.5*x00_med,dlnx0,n)

    n = nf1
    dt = (frac1-1.*(1.-frac1))*frac1*(t1-t0)/n
    dp = (frac1-2.*(1.-frac1))*(p1-p0)/n
    theta_grid,dt0 = tcompile(theta_grid,dt,dt0,n)
    phi_grid,dp0 = tcompile(phi_grid,dp,dp0,n)
    x00_grid,dlnx0 = xcompile(x00_grid,0.8*x00_med,dlnx0,n)

    n = nf3
    dt = 0.33*frac1*(t1-t0)/n
    dp = 0.5*180.0/n
    theta_grid,dt0 = tcompile(theta_grid,dt,dt0,n)
    phi_grid,dp0 = tcompile(phi_grid,dp,dp0,n)
    x00_grid,dlnx0 = xcompile(x00_grid,1.2*x00_med,dlnx0,n)

    if(frame_max<=0): frame_max=x00_grid.size-1;
    if(frame_min<0): frame_min=0;
    i_loop_list = np.arange(frame_min,frame_max+1,1)
    for i in i_loop_list:
        x00=x00_grid[i]; theta_0=theta_grid[i]; phi_0=phi_grid[i]; 
        addlayer='';set_added_layer_alpha=0.0;lightk=0;dustgas=1;ctable='heat_purple'
        threecolor=1; invert_colors=0;
        if(key=='star'):
            dynr = 2.e3; maxd = 1.1e-2 * (x00 / 30.)**(-1.5); lightk=0; 
            dustgas=0.67;
            
            qq=1.e10; maxd*=10.; dynr*=10.;
            maxd=1.e20; dynr=1.e6;
            
            if (add_key=='nodust'): 
                dustgas = 1.e-5
            if (add_key=='CO'): 
                addlayer='CO'; ctable='heat_blue'; set_added_layer_alpha=0.4;
            if (add_key=='Halpha'): 
                addlayer='Halpha'; ctable='heat_green'; set_added_layer_alpha=0.3;
            if (add_key=='SFR'):
                addlayer='SFR'; ctable='heat_purple'; set_added_layer_alpha=0.3;
            if (add_key=='Xray'):
                addlayer='Xray'; ctable='heat_redyellow'; set_added_layer_alpha=0.4;
            if (add_key=='Zmetal'):
                addlayer='Zmetal'; ctable='heat_orange'; set_added_layer_alpha=0.8;
        if(key=='gas'):
            dynr = 1.e3; maxd = 6.0e-3 * (x00 / 30.)**(-1.5);
            threecolor=0; maxd*=10.0; dynr*=3.; lightk=0; invert_colors=0;
            #threecolor=1; lightk=1; ## not good for edge-on projections inside disk
        
        added_label='_full'
        if(add_key!=''): added_label='_'+add_key
        fname=omaster+'/'+sdir+'_'+key+added_label+'_f'+frame_ext(i)

        out0,out1 = im_make.image_maker(sdir,snum,
            snapdir_master=smaster,outdir_master=omaster,pixels=pixels,
            do_with_colors=1,threecolor=threecolor, project_to_camera=1,cosmo=1,
            xrange=np.array([-1.,1.])*x00,yrange=np.array([-1.,1.])*x00,
            dynrange=dynr,maxden=maxd, show_gasstarxray=key,include_lighting=lightk,
            theta=theta_0,phi=phi_0,filename_set_manually=fname,set_added_layer_ctable=ctable,
            add_gas_layer_to_image=addlayer,dust_to_gas_ratio_rescale=dustgas,set_added_layer_alpha=set_added_layer_alpha,
            invert_colors=invert_colors,
            show_scale_label=0,show_time_label=0,spin_movie=1,
            scattered_fraction=scattered_fraction,
            h_rescale_factor=h_rescale_factor,h_max=h_max);
                
                
def check_input():
    infiname='m12v_3_Jan2014/m12v_3_Jan2014_s0440_t011_star_N3c.dat'
    infi=h5py.File(infiname,'r')
    image24 = np.array(infi["image24"])
    image24x = np.array(infi["image24_extralayer"])
    massmap = np.array(infi["massmap"])
    infi.close()
    
    print np.shape(image24)
    print np.shape(image24x)
    print np.max(image24x)
    
    viscolors.load_my_custom_color_tables();
    
    pylab.imshow(image24,origin='lower',interpolation='bicubic',aspect='normal',
        rasterized=True,zorder=0);
    pylab.imshow(image24x,origin='lower',interpolation='bicubic',aspect='normal',
        cmap='heat_purple',alpha=0.9,rasterized=True,zorder=0)

    pylab.savefig('m12v_3_Jan2014/tst.pdf',rasterized=True)#,bbox_inches='tight',pad_inches=0)


def xcompile(x_grid_previous,x_final,dlnx_previous,n):
    x_0 = x_grid_previous[-1]
    x_f = x_final
    dx_desired = np.log(x_f/x_0) / n
    lnx_grid_0 = np.log(x_grid_previous)
    lnx_new, dlnx_new = tcompile(lnx_grid_0,dx_desired,dlnx_previous,n)
    return np.exp(lnx_new), dlnx_new
    #return np.concatenate([x00_grid, x00_i + (x00_f-x00_i)/n*np.arange(0,n)])
    #return np.concatenate([x00_grid, x00_i * np.exp( np.log(x00_f/x00_i)/n*np.arange(0,n) )])

def tcompile(t_grid,dt_desired,dt_previous,n):
    f_transition = 0.1
    dt_0 = dt_previous
    dt_1 = (dt_desired - 0.5*f_transition*dt_previous)/(1.0 - 0.5*f_transition)
    dt = np.zeros(n)
    ng = 1.*np.arange(1,n+1)
    ok = (ng <= f_transition*n)
    dt[ok] = dt_0 + (dt_1-dt_0) * (1.*ng[ok])/(f_transition*n)
    ok = (ng > f_transition*n)
    dt[ok] = dt_1
    tg = np.cumsum(dt)
    t_grid_new = np.concatenate([t_grid, t_grid[-1] + tg])
    return t_grid_new, dt_1


def frame_ext(snum,four_char=1):
	ext='00'+str(snum);
	if (snum>=10): ext='0'+str(snum)
	if (snum>=100): ext=str(snum)
	if (four_char==1): ext='0'+ext
	if (snum>=1000): ext=str(snum)
	return ext;
