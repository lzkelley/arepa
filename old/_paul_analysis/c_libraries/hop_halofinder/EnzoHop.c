/************************************************************************
* Copyright (C) 2008-2011 Matthew Turk.  All Rights Reserved.
*
* This file is part of yt.
*
* yt is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*
************************************************************************/

//
// EnzoHop
//   A module for running HOP halo finding on a set of particles
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#include <ctype.h>
#include <assert.h>
#include "kd.h"
#include "hop.h"
#include "slice.h"
#include "smooth.h"


void initgrouplist(Grouplist *g);
void hop_main(KD kd, HC *my_comm, float densthres);
void regroup_main(float dens_outer, HC *my_comm);


int hop_shared(
    int num_particles_in, // number of particles in vectors
    float* xpos, 
    float* ypos, 
    float* zpos, 
    float* mass, // vectors with appropriate quantities
    int* particle_group_id, // vector for output (group ID of membership of each input particle)
    float thresh_in) // thres=160.0, recommended
{
    int i;
    int num_particles;
    float thresh;
    num_particles = (int)num_particles_in;
    num_particles += 0;
    thresh = (float)thresh_in;
    thresh *= 1.00000000000000001;
    
    if (num_particles < 0) goto _fail;
    float normalize_to = 1.0;
    printf("np %d thres %g \n",num_particles,thresh);

    float totalmass = 0.0;
    for(i = 0; i < num_particles; i++) totalmass += mass[i];
    totalmass /= normalize_to;

    /* initialize the kd hop structure */
    KD kd;
    int nBucket = 16, kdcount = 0;
    kdInit(&kd, nBucket);
    kd->nActive = num_particles;
    kd->p = malloc(sizeof(PARTICLE)*num_particles);
    kd->np_densities = malloc(sizeof(npy_float64)*num_particles);
    kd->np_masses = malloc(sizeof(npy_float64)*num_particles);
    kd->np_pos[0] = malloc(sizeof(npy_float64)*num_particles);
    kd->np_pos[1] = malloc(sizeof(npy_float64)*num_particles);
    kd->np_pos[2] = malloc(sizeof(npy_float64)*num_particles);
    if (kd->p == NULL) {
        fprintf(stderr, "failed allocating particles.\n");
        goto _fail;
    }
  
 	/* Copy positions into kd structure. */
    fprintf(stdout, "Copying arrays for %d particles\n", num_particles);
	for (i = 0; i < num_particles; i++) 
	{
    kd->np_masses[i] = mass[i];
    kd->np_pos[0][i] = xpos[i];
    kd->np_pos[1][i] = ypos[i];
    kd->np_pos[2][i] = zpos[i];
    kd->p[i].np_index = i;
    kd->np_densities[i] = 1.0;
    }
    kd->totalmass = totalmass;

    /*  // use these lines for debugging only //
	for (i = 0; i < num_particles; i++) 
	{
    printf("m %g p0 %g p1 %g p2 %g i %d dd %g \n", 
        NP_MASS(kd,i), NP_POS(kd,i,0), NP_POS(kd,i,1), NP_POS(kd,i,2), i, NP_DENS(kd,i) );
    fflush(stdout);
    }
    printf("mtot = %g \n",kd->totalmass);
    fflush(stdout);
    */ 

    HC my_comm;
    my_comm.s = newslice();
    my_comm.gl = (Grouplist*)malloc(sizeof(Grouplist));
    if(my_comm.gl == NULL) {
        fprintf(stderr, "failed allocating Grouplist\n");
        goto _fail;
    }
    initgrouplist(my_comm.gl);

    fprintf(stderr, "Calling hop... %d %0.3e\n",num_particles,thresh);
    hop_main(kd, &my_comm, thresh);

    fprintf(stderr, "Calling regroup... %d %0.3e\n",num_particles,thresh);
    printf("thresh %g \n",thresh);
    regroup_main(thresh, &my_comm);

    // Now we need to get the groupID, realID and the density.
    // This will give us the index into the original array.
    // Additionally, note that we don't really need to tie the index
    // back to the ID in this code, as we can do that back in the python code.
    // All we need to do is provide density and group information.
    
    // Tags (as per writetagsf77) are in gl.s->ntag+1 and there are gl.s->numlist of them.
//    PyArrayObject *particle_group_id = (PyArrayObject *)
//            PyArray_SimpleNewFromDescr(1, PyArray_DIMS(xpos),
//                    PyArray_DescrFromType(NPY_INT32));
    
    for (i = 0; i < num_particles; i++) {
      // tag is in gl.s->ntag[i+1]
      particle_group_id[i] = my_comm.s->ntag[i+1];
//      if(particle_group_id[i] != 0) {printf("pid = %d %d \n",i,particle_group_id[i]); fflush(stdout);}
//      *(npy_int32*)(PyArray_GETPTR1(particle_group_id, i)) =
//            (npy_int32) my_comm.s->ntag[i+1];
    }

	kdFinish(kd);
    free(my_comm.gl);
    free_slice(my_comm.s);
    
    return 1;


_fail:
    if(kd->p!=NULL)free(kd->p);
    return -1;
}
