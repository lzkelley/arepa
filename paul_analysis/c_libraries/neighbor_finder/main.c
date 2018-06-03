#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "proto.h"

void allocate_3d(int Nsize);
void set_particle_pointer(int Nsize, float* x, float* y, float* z, float* target_weights);
void free_memory_3d(int Nsize);

struct particle_3d
{
    float Pos[3];
    float target_weights;
} **P3d;

/* extremely fast approximation function for the exponential,
 useful here since fractional accuracy errors are smaller than the kernel sources anyways */
inline float fast_exp(float y) {
    float d;
    *((int*)(&d) + 0) = 0;
    *((int*)(&d) + 1) = (int)(1512775 * y + 1072632447);
    return d;
}


// revised call for python calling //
//int stellarhsml(int argc,void *argv[])
int kde(int N_in_search, float* x_search, float* y_search, float* z_search,
                 int N_in_find, float* x, float* y, float* z, float* target_weights,
                 float box_size_x, float box_size_y, float box_size_z,
                 int NumDims, int DesNgb, float Hmin, float Hmax, float* H_OUT, float* NGB_OUT)
{
    int *ngblist; float xyz[3], *r2list, h2, dummy[3];
    float h_guess, h_guess_0, r2, dx, dy, dx_i, dy_i, dx_n, dy_n, i_x_flt, i_y_flt, d_ij, h, hmax, hmin;
    float x2_n, y2_n, r2_n, h2_i, wk, wt_sum, hkernel_over_hsml_to_use, *Kernel;
    float dpi=3.1415926535897932384626433832795;
    long i, j, n, k,imin,imax,jmin,jmax,N_KERNEL_TABLE, N_in = N_in_find;
    int ngbfound;
    allocate_3d(N_in);

    if(box_size_x<=0) box_size_x=1.0e37;
    if(box_size_y<=0) box_size_y=1.0e37;
    if(box_size_z<=0) box_size_z=1.0e37;
    BoxSizeX = box_size_x;
    BoxSizeY = box_size_y;
    BoxSizeZ = box_size_z;

    printf("Entering neighbor search routine: \n");
    printf("  Number of dimensions in which to do search (1=x only, 2=xy, 3=xyz) = %d \n",NumDims);
    printf("  Number of searching points = %d \n",N_in_search);
    printf("  Number of potential neighbor points = %d \n",N_in);
    printf("  Desired number of neighbors (must be >0) = %d \n",DesNgb);
    printf("  Minimum Search Radius (>=0) = %g \n",Hmin);
    printf("  Maximum Search Radius (>0) = %g \n",Hmax);
    printf("  Size of periodic box to search (set to arbitrarily large values, or <0, for non-periodic) in x/y/z = %g/%g/%g \n",box_size_x,box_size_y,box_size_z);
    fflush(stdout);
    
    
    // build a kernel lookup table to save on the calculation below //
    hkernel_over_hsml_to_use = 1.4;
    //hkernel_over_hsml_to_use = 2.0;
    // default = 1 (integrate out to kernel), but low-density regions are represented
    //   more accurately if this is larger (~2); code expense increases though!
    N_KERNEL_TABLE = 10000;
    Kernel = calloc(N_KERNEL_TABLE+1, sizeof(double));
    dx_n=(hkernel_over_hsml_to_use*hkernel_over_hsml_to_use)/((double)N_KERNEL_TABLE); r2_n=0.;
    double kernel_spacing_inv = 1./dx_n;
    for(n=0;n<N_KERNEL_TABLE;n++)
    {
        // radius (to save sqrt operation we're interpolating in r^2 //
        // approximate gaussian for projected, integrate kernel: //
        /* quintic spline (approximate) */
        if(NumDims==1) {wk = sqrt(270./(7.*dpi)) * exp(-(135./14.)*r2_n);}
        if(NumDims==2) {wk = (135./(14.*dpi)) * exp(-(135./14.)*r2_n);}
        if(NumDims==3) {wk = (405./(14.*dpi))*sqrt(15./(14.*dpi)) * exp(-(135./14.)*r2_n);}
        //wk = (135./(14.*dpi)) * exp(-(135./14.)*r2_n/(1.+sqrt(r2_n)));
        // this has more extended tails, designed to reduce artifacts in low-density regions
        //  (where we would really want to use a proper volume-render)
        // wk = 1.91 * exp(-5.56*r2_n); // cubic spline
        // cubic spline kernel
        //h2=(1.-h); if(h<=0.5) {wk=(1.-6.*r2_n*h2);} else {wk=2.*h2*h2*h2;} wk*=40./(7.0*dpi);
        Kernel[n] = wk;
        r2_n += dx_n; // radius at this point in the table
    }
    Kernel[N_KERNEL_TABLE]=0;
    
    
    ngb3d_treeallocate(N_in, 2*N_in);
    set_particle_pointer(N_in, x,y,z, target_weights);
    ngb3d_treebuild((float **)&P3d[1], N_in, 0, dummy, dummy);
    h_guess = Hmax/150.0e0; h_guess_0=h_guess;
    int Desngb_tolerance = DesNgb/4;
    for(i=0;i<N_in_search;i++)
    {
        xyz[0]=x_search[i];
        xyz[1]=y_search[i];
        xyz[2]=z_search[i];
        h2=ngb3d_treefind( xyz, DesNgb, Desngb_tolerance, 1.04*h_guess, &ngblist, &r2list, Hmin, Hmax, &ngbfound);
        
        /* loop over the returned particle list */
        float wk_sum = 0;
        //printf("h2=%g ngb=%d \n",h2,ngbfound); fflush(stdout);
        for(n=0;n<ngbfound;n++)
        {
            r2 = (float)r2list[n];
            h2_i = 1. / ((float)h2);
            //printf("n=%d r2=%g h2=%g kinv=%g \n",n,r2,h2,kernel_spacing_inv);
            if(r2 < h2)
            {
                // use the kernel lookup table compiled above for the weighting //
                double h_i = sqrt(h2_i);
                r2_n = ((float)r2) / ((float)h2) * ((float)kernel_spacing_inv); // ok now have the separation in units of hsml, then kernel_table //
                k = (long)r2_n;
                double prefac;
                if(NumDims==1) prefac = h_i;
                if(NumDims==2) prefac = h2_i;
                if(NumDims==3) prefac = h2_i * h_i;
                prefac *= target_weights[ngblist[n]];
                //printf("n=%d r2=%g h2=%g r2_n=%g k=%ld \n",n,r2,h2,r2_n,k);
                wk = prefac * (Kernel[k] + (Kernel[k+1]-Kernel[k])*(r2_n-k)); // ok that's the weighted result
                wk_sum += wk;
            }
        }
        
        if(!(i%10000))
        {
            printf("i=%d h=%g xyz=%g|%g|%g ngb=%d wk_sum=%g\n",
                   i,sqrt(h2),xyz[0],xyz[1],xyz[2],ngbfound,wk_sum); fflush(stdout);
        }
        H_OUT[i] = sqrt(h2);
        NGB_OUT[i] = wk_sum;
    }
    
    ngb3d_treefree();
    free_memory_3d(N_in);
    printf("done\n");
    return 0;
}



void set_particle_pointer(int Nsize, float* x, float* y, float* z, float* target_weights)
{
    int i;
    float *pos;
    for(i=1;i<=Nsize;i++)
    {
        P3d[i]->Pos[0] = x[i-1];
        P3d[i]->Pos[1] = y[i-1];
        P3d[i]->Pos[2] = z[i-1];
        P3d[i]->target_weights = target_weights[i-1];
    }
}


void allocate_3d(int Nsize)
{
    printf("allocating memory...\n");
    int i;
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


void free_memory_3d(int Nsize)
{
    if(Nsize>0)
    {
        free(P3d[1]);
        P3d++;
        free(P3d);
    }
}





