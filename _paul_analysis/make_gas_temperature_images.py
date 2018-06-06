

import plotting.images.gas_images as gas_images
import matplotlib.pyplot as plt
import numpy as np
import glob

from . simread import readsubfHDF5 as subf
# import simread.readsnapHDF5 as ws

npixels = 720
boxsize = 10.0

run_names = ['explicit_feedback_256', 'explicit_feedback_256_soft',
             'explicit_feedback_512', 'explicit_feedback_512_soft_amd']

for run_base in run_names:
    print("\n\n\n")
    print('../' + run_base + "/output/snap*")
    snaplist = glob.glob("../" + run_base + "/output/snap*")
    n_snaps = len(snaplist)
    print(snaplist)
    for snapnum in [10]:
        cat = subf.subfind_catalog(
            '../' + run_base + '/', snapnum, keysel=['GroupPos', 'GroupLenType'])
        # range(np.min([10, len(cat.GroupLenType[:,0])] ) ):
        for groupnr in range(1):
            if cat.GroupLenType[groupnr, 4] > 500:
                center = cat.GroupPos[groupnr, :]

                fig_, (ax1_, ax2_, ax3_) = plt.subplots(1, 3, figsize=(3, 1))
                image, massmap = gas_images.gas_image(
                    '../' + run_base + '/output/', snapnum, center=center,
                    xrange=[-0.01, 0.01], yrange=[-0.01, 0.01], zrange=[-0.01, 0.01],
                    pixels=npixels, cosmo_wrap=True, massmap=False, projaxis=0,
                    maxden=2e3, dynrange=1e3,
                    unit_length_mpc=True)

                for ijk in range(3):
                    image[:, :, ijk] = np.transpose(image[:, :, ijk])
                ax1_.imshow(image, origin='lower', extent=(-0.01, 0.01, -0.01, 0.01))

                # image,massmap = gas_images.gas_image('../'+run_base+'/output/', snapnum,
                #     center=center, xrange=[-0.01,0.01], yrange=[-0.01,0.01], zrange=[-0.01, 0.01],
                #     pixels=npixels, cosmo_wrap=True, massmap=False, projaxis=1, maxden = 2e3,
                #     dynrange=1e3, unit_length_mpc=True)

                for ijk in range(3):
                    image[:, :, ijk] = np.transpose(image[:, :, ijk])
                ax2_.imshow(image, origin='lower', extent=(-0.01, 0.01, -0.01, 0.01))

                # image, massmap = gas_images.gas_image('../' + run_base + '/output/', snapnum,
                #     center=center, xrange=[-0.01,0.01], yrange=[-0.01,0.01], zrange=[-0.01, 0.01],
                #     pixels=npixels , cosmo_wrap=True, massmap=False, projaxis=2, maxden=2e3,
                #     dynrange=1e3, unit_length_mpc=True)

                for ijk in range(3):
                    image[:, :, ijk] = np.transpose(image[:, :, ijk])
                ax3_.imshow(image, origin='lower', extent=(-0.01, 0.01, -0.01, 0.01))

                # image, massmap = gas_images.gas_image(
                #     '../' + run_base + '/output/', snapnum,
                #     center=center, xrange=[-0.01, 0.01], yrange=[-0.01, 0.01], pixels=npixels,
                #     cosmo_wrap=True, massmap=False, projaxis=1)
                # print(np.min(image), np.min(image[np.isfinite(image)]), np.max(image), np.median(image))
                # ax2_.imshow(image, vmin=0.25, vmax=5.5)
                #
                # image, massmap = gas_images.gas_image(
                #     '../' + run_base + '/output/', snapnum,
                #     center=center, xrange=[-0.01, 0.01], yrange=[-0.01, 0.01], pixels=npixels,
                #     cosmo_wrap=True, massmap=False, projaxis=2)
                # print(np.min(image), np.min(image[np.isfinite(image)]), np.max(image), np.median(image))
                # ax3_.imshow(image, vmin=0.25, vmax=5.5)

                fig_.subplots_adjust(left=0.0, bottom=0.0, top=1.0, right=1.0, hspace=0, wspace=0)
                ax1_.axis('off')
                ax2_.axis('off')
                ax3_.axis('off')

                filename = \
                    '{:s}_gas_central_temperature_image_snap_{:.0f}_group_{:.0f}_res_{:.0f}.png'.format(run_base, snapnum, groupnr, npixels)
                print(filename)
                fig_.savefig('./plots/'+filename, dpi=npixels)
