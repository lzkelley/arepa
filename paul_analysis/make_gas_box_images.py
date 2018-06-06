import plotting.images.gas_images as gas_images
import matplotlib.pyplot as plt
import simread.readsubfHDF5 as subf
import numpy as np
import glob
import simread.readsnapHDF5 as ws


npixels=720
boxsize=10.0

for run_base in [#'explicit_feedback_256', \
                 'explicit_feedback_256_soft', \
#                 'explicit_feedback_512', \
                 'explicit_feedback_512_soft_amd']:         #gas_prod1_runs', 'gas_cool_runs', 'gas_adiab_runs', 'dm_runs']:
    print "\n\n\n"
    print '../'+run_base+'/output/snap*'
    snaplist = glob.glob( '../'+run_base+'/output/snap*')
    n_snaps = len(snaplist)
    print snaplist
    for snapnum in [10]:        #range(n_snaps):        #range(n_snaps):        #range(13):
        cat = subf.subfind_catalog( '../'+run_base+'/', snapnum, keysel=['GroupPos', 'GroupLenType'] )
        for groupnr in range(10):        #range(np.min([10, len(cat.GroupLenType[:,0])] ) ):
            if cat.GroupLenType[groupnr,4] > 500:
                center = cat.GroupPos[groupnr,:]

                if False:               # full box gas images.  Don't delete, but not that useful
                    fig,(ax1,ax2,ax3)=plt.subplots( 1, 3, figsize=(3,1) )

                    image,massmap = gas_images.gas_image( '../'+run_base+'/output/', snapnum, center=center, xrange=[-5.0,5.0], yrange=[-5.0,5.0], pixels=npixels , cosmo_wrap=True, massmap=True,  projaxis=0)
                    print np.min(image), np.min(image[np.isfinite(image)]), np.max(image), np.median(image) 
                    ax1.imshow(image, vmin=0.25, vmax=3.5 )


                    image,massmap = gas_images.gas_image( '../'+run_base+'/output/', snapnum, center=center, xrange=[-5.0,5.0], yrange=[-5.0,5.0], pixels=npixels , cosmo_wrap=True, massmap=True,  projaxis=1)
                    print np.min(image), np.min(image[np.isfinite(image)]), np.max(image), np.median(image) 
                    ax2.imshow(image, vmin=0.25, vmax=3.5 )

                    image,massmap = gas_images.gas_image( '../'+run_base+'/output/', snapnum, center=center, xrange=[-5.0,5.0], yrange=[-5.0,5.0], pixels=npixels , cosmo_wrap=True, massmap=True,  projaxis=2)
                    print np.min(image), np.min(image[np.isfinite(image)]), np.max(image), np.median(image) 
                    ax3.imshow(image, vmin=0.25, vmax=3.5 )

                    fig.subplots_adjust( left=0.0, bottom=0.0, top=1.0, right=1.0)
                    ax1.axis('off');    ax2.axis('off');        ax3.axis('off')

                    filename='{:s}_gas_box_image_snap_{:.0f}_group_{:.0f}_res_{:.0f}.png'.format( run_base, snapnum, groupnr, npixels)
                    fig.savefig( './plots/'+filename, dpi=npixels )



                if True:
                    fig_,(ax1_,ax2_,ax3_)=plt.subplots( 1, 3, figsize=(3,1) )

 #                   image,massmap = gas_images.gas_image( '../'+run_base+'/output/', snapnum, center=center, xrange=[-0.01,0.01], yrange=[-0.01,0.01], pixels=npixels , cosmo_wrap=True, massmap=True, projaxis=0)
 #                   image = np.transpose( image )
  #                  print np.min(image), np.min(image[np.isfinite(image)]), np.max(image), np.median(image) 
  #                  ax1_.imshow(image, vmin=0.25, vmax=5.5, extent=(-0.01, 0.01, -0.01, 0.01), origin='lower' )

                    #image,massmap = gas_images.gas_image( '../'+run_base+'/output/', snapnum, center=center, xrange=[-0.01,0.01], yrange=[-0.01,0.01], pixels=npixels , cosmo_wrap=True, massmap=True, projaxis=1)
                    #image = np.transpose( image )
                    #ax2_.imshow(image, vmin=0.25, vmax=5.5, extent=(-0.01, 0.01, -0.01, 0.01), origin='lower'  )

                    #image,massmap = gas_images.gas_image( '../'+run_base+'/output/', snapnum, center=center, xrange=[-0.01,0.01], yrange=[-0.01,0.01], pixels=npixels , cosmo_wrap=True, massmap=True, projaxis=2)
                    #image = np.transpose( image )
                    #ax3_.imshow(image, vmin=0.25, vmax=5.5, extent=(-0.01, 0.01, -0.01, 0.01), origin='lower'  )

                    fig_.subplots_adjust( left=0.0, bottom=0.0, top=1.0, right=1.0, hspace=0, wspace=0)
                    ax1_.axis('off')
                    ax2_.axis('off')
                    ax3_.axis('off')


                    import units.springel_units as units

                    header   = ws.snapshot_header(  '../'+run_base+'/output/snapdir_'+str(snapnum).zfill(3)+'/snapshot_'+str(snapnum).zfill(3) )
                    star_pos = ws.read_block( '../'+run_base+'/output/snapdir_'+str(snapnum).zfill(3)+'/snapshot_'+str(snapnum).zfill(3), "POS ", 4 )
                    star_bt  = ws.read_block( '../'+run_base+'/output/snapdir_'+str(snapnum).zfill(3)+'/snapshot_'+str(snapnum).zfill(3), "GAGE", 4 )
                    star_hsml  = ws.read_block( '../'+run_base+'/output/snapdir_'+str(snapnum).zfill(3)+'/snapshot_'+str(snapnum).zfill(3), "STH ", 4 )
                    star_cp1  = ws.read_block( '../'+run_base+'/output/snapdir_'+str(snapnum).zfill(3)+'/snapshot_'+str(snapnum).zfill(3), "CP1 ", 4 )
                    star_cp2  = ws.read_block( '../'+run_base+'/output/snapdir_'+str(snapnum).zfill(3)+'/snapshot_'+str(snapnum).zfill(3), "CP2 ", 4 )


                    star_cpr1  = ws.read_block( '../'+run_base+'/output/snapdir_'+str(snapnum).zfill(3)+'/snapshot_'+str(snapnum).zfill(3), "CPR1", 4 )
                    star_cpr2  = ws.read_block( '../'+run_base+'/output/snapdir_'+str(snapnum).zfill(3)+'/snapshot_'+str(snapnum).zfill(3), "CPR2", 4 )

                    ages  = units.age_from_a( star_bt )
                    current_age = units.age_from_a( header.time )

                    ages -= current_age
                    if np.min(ages) < 0:  ages -= np.min(ages)

                    for ijk in range(3):
                        star_pos[:,ijk] -= center[ijk]
                        star_pos[ star_pos[:,ijk] >  boxsize/2.0 ,ijk] -= boxsize
                        star_pos[ star_pos[:,ijk] < -boxsize/2.0 ,ijk] += boxsize

                    in_plot_index = (star_pos[:,0] > -0.01) & (star_pos[:,0] < 0.01) & \
                                    (star_pos[:,1] > -0.01) & (star_pos[:,1] < 0.01) & \
                                    (star_pos[:,2] > -0.01) & (star_pos[:,2] < 0.01) & \
                                    (ages < 0.030 )

                    print np.min( star_hsml[ in_plot_index] ) * 1e3, np.max( star_hsml[ in_plot_index ] ) * 1e3, np.median( star_hsml[ in_plot_index ] ) * 1e3 
                    print star_cp2/star_cp1     # expect larger than unity
                    print np.min( star_cp2 / star_cp1), np.max(star_cp2/star_cp1), np.median( star_cp2 / star_cp1 )



                    print star_cpr2/star_cpr1     # expect larger than unity
                    print np.min( star_cpr2 / star_cpr1), np.max(star_cpr2/star_cpr1), np.median( star_cpr2 / star_cpr1 )

                    sys.exit()

                    ax1_.scatter( star_pos[in_plot_index,0], star_pos[in_plot_index,1], s=0.1, linewidths=0.0, c='k',  alpha=0.5 )
                    ax2_.scatter( star_pos[in_plot_index,2], star_pos[in_plot_index,0], s=0.1, linewidths=0.0, c='k',  alpha=0.4 )
                    ax3_.scatter( star_pos[in_plot_index,1], star_pos[in_plot_index,2], s=0.1, linewidths=0.0, c='k',  alpha=0.3 )

                    for ijk in range(3):
                        print star_pos[in_plot_index,ijk]


                    filename='{:s}_gas_central_image_with_stars_snap_{:.0f}_group_{:.0f}_res_{:.0f}.png'.format( run_base, snapnum, groupnr, npixels)
                    print filename

                    fig_.savefig( './plots/'+filename, dpi=npixels )





                if False:       # this has moved to gas_temperature_...
                    fig_,(ax1_,ax2_,ax3_)=plt.subplots( 1, 3, figsize=(3,1) )
        
                    image,massmap = gas_images.gas_image( '../'+run_base+'/output/', snapnum, center=center, xrange=[-0.01,0.01], yrange=[-0.01,0.01], pixels=npixels , cosmo_wrap=True, massmap=False, projaxis=0, \
                                                    #maxden = 50, dynrange = 1e3, 
                                                    unit_length_mpc=True)
                    print "\n\n\n"
                    print image.shape
                    print np.min(image), np.min(image[np.isfinite(image)]), np.max(image), np.median(image) 
                    print "\n\n\n"
                    ax1_.imshow(image ) 
    
    
    
    #                image,massmap = gas_images.gas_image( '../'+run_base+'/output/', snapnum, center=center, xrange=[-0.01,0.01], yrange=[-0.01,0.01], pixels=npixels , cosmo_wrap=True, massmap=False, projaxis=1)
    #                print "\n\n\n"
    #                print image.shape
    #                print np.min(image), np.min(image[np.isfinite(image)]), np.max(image), np.median(image) 
    #                print "\n\n\n"
    #                ax2_.imshow(image, vmin=0.25, vmax=5.5 )
    
    
    #                image,massmap = gas_images.gas_image( '../'+run_base+'/output/', snapnum, center=center, xrange=[-0.01,0.01], yrange=[-0.01,0.01], pixels=npixels , cosmo_wrap=True, massmap=False, projaxis=2)
    #                print "\n\n\n"
    #                print image.shape
    #                print np.min(image), np.min(image[np.isfinite(image)]), np.max(image), np.median(image) 
    #                print "\n\n\n"
    #                ax3_.imshow(image, vmin=0.25, vmax=5.5 )
    
                    fig_.subplots_adjust( left=0.0, bottom=0.0, top=1.0, right=1.0, hspace=0, wspace=0)
                    ax1_.axis('off')
                    ax2_.axis('off')
                    ax3_.axis('off')
                    
    
                    filename='{:s}_gas_central_temperature_image_snap_{:.0f}_group_{:.0f}_res_{:.0f}.png'.format( run_base, snapnum, groupnr, npixels)
                    print filename
                    
                    fig_.savefig( './plots/'+filename, dpi=npixels )

