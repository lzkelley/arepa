#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "overhead.h"
#include "proto.h"
#include "ngbtree3d.h"

#define DEBUG NO
#define NAMELEN 60 

long NUMBER_OF_RAYS;
long Ngas,Nstar,Nbh;

struct particle_3d
{
  float Pos[3];
  float u;
  float rho;
  float hsml;
  float numh;
  float nume;
  float z;
  float mass;
} **P3d;


// re-written call for python, rather than IDL //
/* long getnh(long argc, char **argv) */
// long getnh(long argc, void *argv[]) //
int getnh(int Ngast, int Nstart, int Nbht, float theta, float phi, 
    float* Pos, float* OUT_NH_HOT, float* OUT_NH, float* OUT_Z, 
    float MAX_DISTANCE_FROM0_in, float TREE_MIN_SIZE_in)
{
  long i, j, n_gas, ray_num, N_rays, pid, tmpval; 
  long first_bh_ray, first_star_ray, first_disk_ray, first_bulge_ray; 
  long RAY_ORIGIN_CELL; 
  float dummy[3];
  float nh_prefactor;
  float *on,*onh,*oz; 
  Ray_struct *ray;

  void allocate_3d(void);
  void set_particle_pointer(float* Pos);
  void free_memory_3d(void);
  void FreeCoolMemory(void);
  void free_the_tree(void);

  /*  // brick deprecated with python-style call
  if((argc < 9) || (argc > 11)) {
    fprintf(stderr, "Expected 9-11 arguments (found %ld)\n",argc); 
    exit(0); 
  }
  float theta, phi; 
  int Ngast,Nstart,Nbht;
  float *Pos;
  float *OUT_NH,*OUT_NH_HOT,*OUT_Z; 	// output NH //

  Ngast = *(int *)argv[0];
  Nstart = *(int *)argv[1];
  Nbht = *(int *)argv[2];
  theta = *(float *)argv[3]; 
  phi = *(float *)argv[4]; 
  Pos = (float *)argv[5];
  // returned quantities //
  OUT_NH_HOT = (float *)argv[6]; 
  OUT_NH = (float *)argv[7];
  OUT_Z = (float *)argv[8]; 
  if(argc > 9) {
  	printf("Detected variable input 8 :: \n");
  	MAX_DISTANCE_FROM0=*(float *)argv[9]; 
  	printf("MAXDIST RESET = %g \n",MAX_DISTANCE_FROM0);
  } else {
  	MAX_DISTANCE_FROM0=1.0;
  }
  if(argc > 10) {
  	printf("Detected variable input 9 :: \n");
  	TREE_MIN_SIZE=*(float *)argv[10]; 
  	printf("MINCELL RESET = %g \n",TREE_MIN_SIZE);
  } else {
  	TREE_MIN_SIZE=0.001;
  }
  */
  MAX_DISTANCE_FROM0 = MAX_DISTANCE_FROM0_in;
  TREE_MIN_SIZE = TREE_MIN_SIZE_in;
  All.N_gas=(long)Ngast;
  All.N_star=(long)Nstart;
  All.N_bh=(long)Nbht;

  All.N_halo = 0; 
  All.N_disk = 0; 
  All.N_bulge = 0; 
  All.N_total = All.N_gas + All.N_star + All.N_bh; 
  Ngas = All.N_gas; 
  Nstar = All.N_star; 
  Nbh = All.N_bh; 

  fprintf(stdout,"theta, phi, maxdist, treemin = %.2f %.2f %g %g\n", 
    theta, phi, MAX_DISTANCE_FROM0, TREE_MIN_SIZE); 
  fprintf(stdout, "Ngas = %ld Nstar = %ld Nbh = %ld Ntotal = %ld\n", 
  	(long)Ngas, (long)Nstar, (long)Nbh, (long)All.N_total);

  if (All.N_bh < 1) { 
    fprintf(stdout,"Warning: No black holes! \n"); 
  }

  /* check that we don't go try to integrate exactly along cell boundaries */

  if (sin(theta) == 0.) theta += 1.e-30; 
  if (sin(theta) == 1.) theta += 1.e-30;
  if (sin(theta) == -1.) theta += 1.e-30;
  if (sin(phi) == 0.) phi += 1.e-30; 
  if (sin(phi) == 1.) phi += 1.e-30;
  if (sin(phi) == -1.) phi += 1.e-30;

  allocate_3d();
  set_particle_pointer(Pos);
  allocate_gas(); 
  allocate_star(); 
  allocate_bh(); 

  fprintf(stdout, "star, gas, and bh allocation done too\n"); 

  /* --------------------------------------------------

     OK, now we're ready to use whatever variables we
     need.

     We can access the various fields via:

     P3d[i+1]->Pos[0], or P3d[30]->rho= 100.0


  here's an example:

  printf("\n\n");
  printf("Here are some of the P3d variables\n");
  printf("P3d[1]->Pos[0,1,2]= %g|%g|%g\n",P3d[1]->Pos[0],P3d[1]->Pos[1],P3d[1]->Pos[2]);
  printf("P3d[3003]->u= %g\n",P3d[3003]->u);
  printf("P3d[%ld]->Pos[0,1,2]= %g|%g|%g\n",Ngas+Nstar,P3d[Ngas+Nstar]->Pos[0],P3d[Ngas+Nstar]->Pos[1],P3d[Ngas+Nstar]->Pos[2]);
  printf("P3d[%ld]->Pos[0,1,2]= %g|%g|%g\n",Ngas+Nstar+1,P3d[Ngas+Nstar+1]->Pos[0],P3d[Ngas+Nstar+1]->Pos[1],P3d[Ngas+Nstar+1]->Pos[2]);
  printf("\n\n");
   */


  /* Fill in the structures that are normally done in readsnap */ 
  
  for (i=0; i<All.N_gas; i++) {
    PG[i].mass = P3d[i+1]->mass;
    PG[i].pos[0] = P3d[i+1]->Pos[0]; 
    PG[i].pos[1] = P3d[i+1]->Pos[1];
    PG[i].pos[2] = P3d[i+1]->Pos[2];
    PG[i].u = P3d[i+1]->u; 
    PG[i].rho = P3d[i+1]->rho; 
    PG[i].hsml = P3d[i+1]->hsml; 
    PG[i].nh = P3d[i+1]->numh; 
    PG[i].ne = P3d[i+1]->nume; 
    PG[i].z = P3d[i+1]->z; 
	//printf(" check gas masses = %e \n",PG[i].mass);
	//printf(" check gas metals = %e \n",PG[i].z);
        //printf(" check gas energy = %e \n",PG[i].u);
        //printf(" check gas ioniza = %e \n",PG[i].ne);
  } 
  tmpval = All.N_gas+All.N_star; 
  for (i=All.N_gas,pid=0; i<tmpval; i++,pid++) {
    PS[pid].pos[0] = P3d[i+1]->Pos[0]; 
    PS[pid].pos[1] = P3d[i+1]->Pos[1];
    PS[pid].pos[2] = P3d[i+1]->Pos[2];
  }
  for (i=tmpval,pid=0; i<All.N_total; i++,pid++) {
    PBH[pid].pos[0] = P3d[i+1]->Pos[0]; 
    PBH[pid].pos[1] = P3d[i+1]->Pos[1];
    PBH[pid].pos[2] = P3d[i+1]->Pos[2];
  }

if(DEBUG)
{
  fprintf(stderr, "\n"); 
  fprintf(stderr, "Quick test: \n"); 
  fprintf(stderr, "PG[0].pos = %g %g %g\n", PG[0].pos[0], PG[0].pos[1], PG[0].pos[2]); 
  fprintf(stderr, "PG[1].pos = %g %g %g\n", PG[1].pos[0], PG[1].pos[1], PG[1].pos[2]); 
  fprintf(stderr, "PG[3002].pos = %g %g %g PG[3002].u = %g\n", PG[3002].pos[0], PG[3002].pos[1], PG[3002].pos[2], PG[3002].u); 
  fprintf(stderr, "PS[0].pos = %g %g %g\n", PS[0].pos[0], PS[0].pos[1], PS[0].pos[2]); 
  fprintf(stderr, "PS[1].pos = %g %g %g\n", PS[1].pos[0], PS[1].pos[1], PS[1].pos[2]); 
  fprintf(stderr, "PS[163438].pos = %g %g %g\n", PS[163438].pos[0], PS[163438].pos[1], PS[163438].pos[2]); 
  fprintf(stderr, "\n\n");
}

  /* free_memory_3d();  */

  /* 
   *  initialize the tree 
   */

  setup_lineofsight_overhead();
  if(DEBUG) for (i=0; i<10; i++) printf(" mass = %e %e \n",PG[i].mass,P3d[i+1]->mass);



#ifndef DIRECT_RAY_INTEGRATION
  printf("Using: Faster TREE method\n");
  //float dummy[3];
  ngb3d_treeallocate(All.N_gas, 10*All.N_gas);
  ngb3d_treebuild((float **)&P3d[1], All.N_gas, 0, dummy, dummy);
  allocate_ngblists();
  if (USE_FULL_NEIGHBOR_CALC == 0) initialize_the_tree();
  fprintf(stdout, "tree built... \n"); 
#else
  printf("Using: DIRECT_RAY_INTEGRATION\n");
  ngb3d_treeallocate(All.N_gas, 2*All.N_gas);
  ngb3d_treebuild((float **)&P3d[1], All.N_gas, 0,dummy,dummy);
  printf("Using New 3D tree.\n");
  allocate_ngblists();
#endif

  /* 
   *  initialize the rays 
   */ 

  NUMBER_OF_RAYS  = (All.N_bh + All.N_star + All.N_disk + All.N_bulge);
  N_rays = NUMBER_OF_RAYS;
  printf("NUMBER_OF_RAYS = %ld, N_rays = %ld, sizeof(Ray_struct) = %ld\n", NUMBER_OF_RAYS, N_rays, sizeof(Ray_struct) ); 
  printf("Initializing (%ld) Rays... \n",N_rays);
  /* ray = calloc(N_rays,sizeof(Ray_struct)); */
  /* ray = calloc(N_rays,sizeof(float)); */
  if(!(ray = malloc(N_rays*sizeof(struct Ray_s))))
     {
        printf("failed to allocate memory: N_rays= %ld  size= %ld\n",N_rays,N_rays*sizeof(struct Ray_s));
        exit(0);
     }
  if (DEBUG) fprintf(stdout, "Done allocating memory for rays\n"); 
  ALL_CELL_COUNTER = 0;
  nh_prefactor = (All.UnitMass_in_g/PROTONMASS)/(SQR(All.UnitLength_in_cm));
  printf("nh_prefactor= %g\n",nh_prefactor);

  /* 
   *  set up the rays 
   */ 

  first_bh_ray = 0; 
  first_star_ray = first_bh_ray + All.N_bh; 
  first_disk_ray = first_star_ray + All.N_star; 
  first_bulge_ray = first_disk_ray + All.N_disk; 

/*
 * allocate memory for the IDL arrays 
 */


  /* rays from BH particles */ 
  for (ray_num=first_bh_ray,pid=0; ray_num<first_star_ray; ray_num++,pid++) { 
    for (j=0; j<3; j++) ray[ray_num].pos[j] = PBH[pid].pos[j]; 
  } 
  /* rays from star particles */ 
  for (ray_num=first_star_ray,pid=0; ray_num<first_disk_ray; ray_num++,pid++) { 
    for (j=0; j<3; j++) ray[ray_num].pos[j] = PS[pid].pos[j]; 
  } 
  /* rays from disk particles */ 
  for (ray_num=first_disk_ray,pid=0; ray_num<first_bulge_ray; ray_num++,pid++) { 
    for (j=0; j<3; j++) ray[ray_num].pos[j] = PD[pid].pos[j]; 
  } 
  /* rays from bulge particles */ 
  for (ray_num=first_bulge_ray,pid=0; ray_num<NUMBER_OF_RAYS; ray_num++,pid++) { 
    for (j=0; j<3; j++) ray[ray_num].pos[j] = PB[pid].pos[j]; 
  } 
  /* initialize the remaining rays */ 
  for (ray_num=NUMBER_OF_RAYS; ray_num<N_rays; ray_num++) { 
    for (j=0; j<3; j++) ray[ray_num].pos[j] = 0.; 
  } 
  if (DEBUG) fprintf(stdout, "Done setting up the ray IDs and positions\n"); 

  /* 
   *  now perform the calculations 
   */ 

  for (ray_num=0; ray_num<NUMBER_OF_RAYS; ray_num++) { 
    if(!(ray_num%10000)) {printf("%ld ....",ray_num); fflush(stdout);}

    ray[ray_num].theta = theta; 
    ray[ray_num].phi = phi; 
    ray[ray_num].n_hat[0] = sin(ray[ray_num].theta) * cos(ray[ray_num].phi);
    ray[ray_num].n_hat[1] = sin(ray[ray_num].theta) * sin(ray[ray_num].phi);
    ray[ray_num].n_hat[2] = cos(ray[ray_num].theta);
    ray[ray_num].nh=1.e-40;
    ray[ray_num].nh_hot=1.e-40; 
    ray[ray_num].Z=0.; 
    ray[ray_num].neutral_frac=0.;

    if(DEBUG) 
        if(!(ray_num%100) || ray_num<10) 
        {
        printf("..%ld..",ray_num);
        printf("to find: num %ld px py pz %g %g %g \n",ray_num,ray[ray_num].pos[0],ray[ray_num].pos[1],ray[ray_num].pos[2]);
        fflush(stdout);
        }
        
#ifndef DIRECT_RAY_INTEGRATION
    if (USE_FULL_NEIGHBOR_CALC == 0) {
	RAY_ORIGIN_CELL = find_cell_from_scratch(ray[ray_num].pos); 
    } else {
	RAY_ORIGIN_CELL = 0;
    }
#else
    RAY_ORIGIN_CELL = 0;
#endif
    
    integrate_ray_to_escape(&ray[ray_num],RAY_ORIGIN_CELL);

    ray[ray_num].Z /= ray[ray_num].nh;
    ray[ray_num].neutral_frac /= ray[ray_num].nh;
    ray[ray_num].nh   *= nh_prefactor;
    ray[ray_num].nh_hot *= nh_prefactor;

    if (!(ray_num%10000)) {
        printf("ray_num %ld origin_cell %ld theta %.1f phi %.1f nh %.3f Z %f\n",ray_num, RAY_ORIGIN_CELL,ray[ray_num].theta,ray[ray_num].phi,log10(ray[ray_num].nh_hot),ray[ray_num].Z);
        fflush(stdout);
    }

    /* write the results to the input arrays */ 
    on = OUT_NH + ray_num; 
    *on = ray[ray_num].nh; 

    onh = OUT_NH_HOT + ray_num; 
    //*onh = ray[ray_num].nh_hot; 
    *onh = ray[ray_num].neutral_frac; 

    oz = OUT_Z + ray_num;
    *oz = ray[ray_num].Z;
/*
    OUT_NH[ray_num] = ray[ray_num].nh_hot; 
*/


if (DEBUG) {
if (!(ray_num%100)) {
fprintf(stderr,"ray_num %ld origin_cell %ld theta %.1f phi %.1f nh %.3f Z %f\n",ray_num, RAY_ORIGIN_CELL,ray[ray_num].theta,ray[ray_num].phi,log10(ray[ray_num].nh_hot),ray[ray_num].Z);
fprintf(stderr, "(%ld) %f -- %f -- %f\n", ray_num, OUT_NH[ray_num], *on, ray[ray_num].nh_hot); 
}
if(ray_num<10)
printf("ray_num= %ld  OUT_NH= %g  OUT_Z= %g\n",ray_num,ray[ray_num].nh_hot,ray[ray_num].Z);
}


  }
    printf("done with main loop, freeing memory \n");fflush(stdout);

  free_the_tree();
    if(DEBUG){printf(".1.");fflush(stdout);}
  ngb3d_treefree();
    if(DEBUG){printf(".2.");fflush(stdout);}
  free_ngblists();
    if(DEBUG){printf(".3.");fflush(stdout);}
  free(PG);
    if(DEBUG){printf(".4.");fflush(stdout);}
  free(PS);
    if(DEBUG){printf(".5.");fflush(stdout);}
  free(PBH);
    if(DEBUG){printf(".6.");fflush(stdout);}
  FreeCoolMemory();
    if(DEBUG){printf(".7.");fflush(stdout);}
  free(ray);
    if(DEBUG){printf(".8.");fflush(stdout);}
  free_memory_3d();
    if(DEBUG){printf(".9.\n");fflush(stdout);}

  printf("\n Finished call to getnh.so ! \n");
  return 0; 
}


/* grab-bag for random initializations, etc (as 'begrun' from main code) */
void setup_lineofsight_overhead(void)
{
	set_sph_kernel();
	InitCool();
	set_All_struct_terms();
}


/* Generate the ray angles (generated such that they cover solid angle 
 *  equally, trivial to switch around to focus on particular direction)
 */
void generate_ray_angles(Ray_struct *R, long N_rays)
{
  long i0, theta_i, phi_i, n_ang_max = (long)(sqrt(N_rays));
  float theta, costheta, phi;
  float d_costheta = 2.0 / ((float)n_ang_max);
  float d_phi = 2.0 * PI / ((float)n_ang_max);
  
  /* Move in increments of cos(theta) and phi */
  i0 = 0;
  for (theta_i = 0; theta_i < n_ang_max; theta_i++) {
    costheta = 1.0 - (0.5 + (float)theta_i) * d_costheta;
      if (costheta < -1.0) costheta = -1.0;
      if (costheta >  1.0) costheta =  1.0;
    theta = acos(costheta);
    for (phi_i = 0; phi_i < n_ang_max; phi_i++) {
      phi = (0.5 + (float)phi_i) * d_phi;
      
      R[i0].theta = theta;
      R[i0].phi = phi;  
	  i0++;
  }}
}






void set_particle_pointer(float* Pos)
{
  long i;
  float *pos;

  for(i=1,pos=Pos;i<=Ngas;i++) 
    {
      
      P3d[i]->Pos[0] = pos[0];
      P3d[i]->Pos[1] = pos[1];
      P3d[i]->Pos[2] = pos[2];
  
      P3d[i]->u = pos[3];
      P3d[i]->rho = pos[4];
      P3d[i]->hsml = pos[5];
      P3d[i]->numh = pos[6];
      P3d[i]->nume = pos[7];
      P3d[i]->z = pos[8];
      P3d[i]->mass = pos[9];

      pos+=10;
    }

  for(i=Ngas+1;i<=(Ngas+Nstar+Nbh);i++)
    {

      P3d[i]->Pos[0] = pos[0];
      P3d[i]->Pos[1] = pos[1];
      P3d[i]->Pos[2] = pos[2];

      P3d[i]->u = 0.0;
      P3d[i]->rho = 0.0;
      P3d[i]->hsml = 0.0;
      P3d[i]->numh = 0.0;
      P3d[i]->nume = 0.0;
      P3d[i]->z = 0.0;
      P3d[i]->mass = 0.0;

      pos+=10;
    }
}

void allocate_3d(void)
{
  long i;
  long Nsize;

  Nsize=Ngas+Nstar+Nbh;

  printf("allocating memory...\n");

  if(Nsize>0)
    {
      if(!(P3d=malloc(Nsize*sizeof(struct particle_3d *))))
        {
          printf("failed to allocate memory. (A)\n");
          exit(0);
        }

      P3d--;   /* start with offset 1 */

      if(!(P3d[1]=malloc(Nsize*sizeof(struct particle_3d)))) 
        {
          printf("failed to allocate memory. (B)\n");
          exit(0);
        }

      for(i=2;i<=Nsize;i++)   /* initiliaze pointer table */
        P3d[i]=P3d[i-1]+1;

    }

  printf("allocating memory...done\n");
}

void free_memory_3d(void)
{ 
  long Nsize;
  Nsize=Ngas+Nstar+Nbh;
  if(Nsize>0)
    {
      free(P3d[1]);
      P3d++;
      free(P3d);
    }
}
