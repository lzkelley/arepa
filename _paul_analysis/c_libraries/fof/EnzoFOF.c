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
// EnzoFOF
//   A module for running friends-of-friends halo finding on a set of particles
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <signal.h>
#include <ctype.h>
#include "kd.h"
#include "tipsydefs.h"


int fof_shared(
    int num_particles_in, // number of particles in vectors
    float* xpos, 
    float* ypos, 
    float* zpos, 
    float* mass, // vectors with appropriate quantities
    int* particle_group_id, // vector for output (group ID of membership of each input particle)
    float link_in) // link=0.2, recommended (linking length relative to mean particle separation)
{
	KDFOF kd;
	int nBucket,i,j;
	float fPeriod[3],fEps;
	int nMembers,nGroup,bVerbose=1;
	int sec,usec;
    int num_particles = (int)num_particles_in;
    num_particles += 0;
    float link = (float)link_in;
    link *= 1.00000000000000001;
	
	/* linking length */
	fprintf(stdout, "Link length is %f\n", link);
	fEps = link;
	/* hard-wire a couple of choices for particle decomposition in the kd-tree */
	nBucket = 16;
	nMembers = 8;
	for (j=0;j<3;++j) fPeriod[j] = 1.0; /* box size assumed to be 1.0 here, hard-wired! */

    /* initialize the kd FOF structure */
	kdInitFoF(&kd,nBucket,fPeriod);
    kd->nActive = num_particles;
    kd->p = (PARTICLEFOF *)malloc(sizeof(PARTICLEFOF)*num_particles);
	assert(kd->p != NULL);
 	/* Copy positions into kd structure. */
    fprintf(stdout, "Filling in %d particles\n", num_particles);
	for (i = 0; i < num_particles; i++) 
	{
	  kd->p[i].iOrder = i;
	  kd->p[i].r[0] = xpos[i];
	  kd->p[i].r[1] = ypos[i];
	  kd->p[i].r[2] = zpos[i];
    }

	
	kdBuildTreeFoF(kd);
	kdTimeFoF(kd,&sec,&usec);
	nGroup = kdFoF(kd,fEps);
	kdTimeFoF(kd,&sec,&usec);
	if (bVerbose) printf("Number of initial groups:%d\n",nGroup);
	nGroup = kdTooSmallFoF(kd,nMembers);
	if (bVerbose) {
		printf("Number of groups:%d\n",nGroup);
		printf("FOF CPU TIME: %d.%06d secs\n",sec,usec);
		}
	kdOrderFoF(kd);

	/* kdOutGroupFoF(kd,ach); */
    // Now we need to get the groupID, realID.
    // This will give us the index into the original array.
    // Additionally, note that we don't really need to tie the index
    // back to the ID in this code, as we can do that back in the python code.
    // All we need to do is group information.
    for (i = 0; i < num_particles; i++) {
      // group tag is in kd->p[i].iGroup
      particle_group_id[i] = kd->p[i].iGroup;
    }

	kdFinishFoF(kd);
    return 1;

_fail:
    if(kd->p!=NULL)free(kd->p);
    return -1;
}
