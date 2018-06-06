#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "proto.h"
#include "ngbtree3d.h"



void allocate_3d(int Nsize);
void set_particle_pointer(int Nsize, float* x, float* y, float* z, float* vx, float* vy, float* vz, float* mass);
void free_memory_3d(int Nsize);

void kernel_main(double u, double hinv3, double hinv4, double* wk, double* dwk, int mode);

struct particle_3d
{
    float Pos[3];
    float Vel[3];
    float mass;
} **P3d;


// revised call for python calling //
int gasvelocitydispersion(int N_gas,
                      float* x_gas,  float* y_gas,  float* z_gas, 
                      float* vx_gas, float* vy_gas, float* vz_gas,
		      float* m_gas,
                      int DesNgb, float Hmax, float* H_OUT)
{
    float h_guess, h2, xyz[3], vxyz[3], dummy[3], h_guess_0;
    int i, ngbfound;
    allocate_3d(N_gas);
    float *r2list;
    int *ngblist;
    int count, maxcount;
    

    double h, Ngb, Rho, dx, dy, dz, dvx, dvy, dvz, r, r2, u, hinv, hinv3, hinv4, wk, dwk, mj_wk;
    double dvx_mean, dvy_mean, dvz_mean;
    double dvx_stdev, dvy_stdev, dvz_stdev;
    int n, j;
    
    printf("N_gas=%d\n",N_gas);
    printf("Hmax=%g\n",Hmax);
    printf("DesNgb=%d\n",DesNgb);
    
    

    // allocate and build tree for gas particle //
    ngb3d_treeallocate(N_gas, 2*N_gas);
    set_particle_pointer(N_gas, x_gas,y_gas,z_gas, vx_gas,vy_gas,vz_gas, m_gas);
    ngb3d_treebuild((float **)&P3d[1], N_gas, 0, dummy, dummy);
    h_guess = Hmax/150.0e0; h_guess_0=h_guess;
    
    // tree now constructed for gas, loop over gas particles to find their hsml //
    for(i=0;i<N_gas;i++)
    {
        xyz[0]=x_gas[i];             //P3d[i+1]->Pos[0]+1.0e-10;
        xyz[1]=y_gas[i];             //P3d[i+1]->Pos[1]+1.0e-10;
        xyz[2]=z_gas[i];             //P3d[i+1]->Pos[2]+1.0e-10;

	count = 0;
	maxcount = 30;
	do{
	  
            h2=ngb3d_treefind( xyz, DesNgb ,1.04*h_guess, &ngblist, &r2list, Hmax, &ngbfound);    
	    if(ngbfound < (DesNgb-2)) h_guess*=1.26;
            if(ngbfound > (DesNgb-2)) h_guess/=1.26;

	    count += 1;
	}while( ( (ngbfound < (DesNgb-2)) || (ngbfound > (DesNgb+2) ) ) && (count < maxcount)  );

        h = sqrt(h2);
        hinv= 1.0/h;
        hinv3= hinv * hinv * hinv;
        hinv4= hinv3 * hinv;
        
        Ngb = 0.0;
        dvx_mean = 0.0;
        dvy_mean = 0.0;
        dvz_mean = 0.0;

        
        for(n=0; n<ngbfound; n++)
        {
            j = ngblist[n];
            dx = xyz[0] - P3d[j+1]->Pos[0];
            dy = xyz[1] - P3d[j+1]->Pos[1];
            dz = xyz[2] - P3d[j+1]->Pos[2];

            dvx= vxyz[0] = P3d[j+1]->Vel[0];
            dvy= vxyz[1] = P3d[j+1]->Vel[1];
            dvz= vxyz[2] = P3d[j+1]->Vel[2];

            r2 = dx * dx + dy * dy + dz * dz;
        
            if(r2 < h2)
            {
                dvx_mean += dvx;
                dvy_mean += dvy;
                dvz_mean += dvz;
                Ngb  += 1.;
            }
        }
        dvx_mean /= Ngb;
        dvy_mean /= Ngb;
        dvz_mean /= Ngb;

        Ngb = 0.0;
        dvx_stdev = 0.0;
        dvy_stdev = 0.0;
        dvz_stdev = 0.0;
        for(n=0; n<ngbfound; n++)
        {
            j = ngblist[n];
            dx = xyz[0] - P3d[j+1]->Pos[0];
            dy = xyz[1] - P3d[j+1]->Pos[1];
            dz = xyz[2] - P3d[j+1]->Pos[2];
            
            dvx= vxyz[0] = P3d[j+1]->Vel[0];
            dvy= vxyz[1] = P3d[j+1]->Vel[1];
            dvz= vxyz[2] = P3d[j+1]->Vel[2];

            r2 = dx * dx + dy * dy + dz * dz;

            if(r2 < h2)
            {
                dvx_stdev += (dvx - dvx_mean) * (dvx - dvx_mean) ;
                dvy_stdev += (dvy - dvy_mean) * (dvy - dvy_mean) ;
                dvz_stdev += (dvz - dvz_mean) * (dvz - dvz_mean) ;
                Ngb  += 1.;
            }
        }
        dvx_stdev /= Ngb;
        dvy_stdev /= Ngb;
        dvz_stdev /= Ngb;
        dvx_stdev = sqrt( dvx_stdev );
        dvy_stdev = sqrt( dvy_stdev );
        dvz_stdev = sqrt( dvz_stdev );


        
            
        if(!(i%10000))
        {
            printf("i=%d hmax=%g h_guess=%g h=%g xyz=%g|%g|%g ngb=%d \n",
                   i,Hmax,h_guess,sqrt(h2),xyz[0],xyz[1],xyz[2],ngbfound); fflush(stdout);
        }
        H_OUT[i] = sqrt( dvx_stdev * dvx_stdev + dvy_stdev * dvy_stdev + dvz_stdev * dvz_stdev) ;

        h_guess = sqrt(h2); // use this value for next guess, should speed things up //
        //if (h_guess>10.*h_guess_0) h_guess=2.*h_guess_0;
    } 
    
    ngb3d_treefree();
    free_memory_3d(N_gas);
    printf("done\n");
    return 0;
}



void set_particle_pointer(int Nsize, float* x, float* y, float* z, float* vx, float* vy, float* vz, float* mass)
{
    int i;
    float *pos;
    for(i=1;i<=Nsize;i++)
    {
        P3d[i]->Pos[0] = x[i-1];
        P3d[i]->Pos[1] = y[i-1];
        P3d[i]->Pos[2] = z[i-1];
        P3d[i]->Vel[0] = vx[i-1];
        P3d[i]->Vel[1] = vy[i-1];
        P3d[i]->Vel[2] = vz[i-1];
        P3d[i]->mass   = mass[i-1];
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


/* Attention: Here we assume that kernel is only called
 with range 0..1 for u as done in hydra or density !!
 Call with mode 0 to calculate dwk and wk
 Call with mode -1 to calculate only wk
 Call with mode +1 to calculate only dwk */

void kernel_main(double u, double hinv3, double hinv4,
                               double *wk, double *dwk, int mode)
{
    if(u < 0.5)
    {
        if(mode >= 0)
        *dwk = u * (18.0 * u - 12.0);
        if(mode <= 0)
        *wk = (1.0 + 6.0 * (u - 1.0) * u * u);
    }
    else
    {
        double t1 = (1.0 - u);
        double t2 = t1 * t1;
        if(mode >= 0)
        *dwk = -6.0 * t2;
        if(mode <= 0)
        *wk = 2.0 * t2 * t1;
    }
    
    if(mode >= 0)
    *dwk *= NORM * hinv4;
    if(mode <= 0)
    *wk *= NORM * hinv3;
    
    return;
}

