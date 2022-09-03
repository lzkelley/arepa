#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "proto.h"
#include "ngbtree3d.h"



void allocate_3d(int Nsize);
void set_particle_pointer(int Nsize, float* x, float* y, float* z, float* mass, float* rho_gas, float* u_gas);
void free_memory_3d(int Nsize);

void kernel_main(double u, double hinv3, double hinv4, double* wk, double* dwk, int mode);

struct particle_3d
{
    float Pos[3];
    float mass;
    float u_gas;
    float rho_gas;
//    int   my_id;
} **P3d;


// revised call for python calling //
//int stellarhsml(int argc,void *argv[])
int stellargastemperature(int N_gas, int N_star,
                      float* x_gas, float* y_gas, float* z_gas, float* m_gas, float* rho_gas, float* u_gas,
                      float* x_star, float* y_star, float* z_star,
                      int DesNgb, float Hmax, float* H_OUT)
{
    float h_guess, h2, xyz[3], dummy[3], h_guess_0;
    int i, ngbfound;
    allocate_3d(N_gas);
    float *r2list;
    int *ngblist;
    int count, maxcount;
    

    double h, Ngb, Rho, Temp, dx, dy, dz, r, r2, u, hinv, hinv3, hinv4, wk, dwk, mj_wk;
    int n, j;
    
    printf("N_gas=%d\n",N_gas);
    printf("N_star=%d\n",N_star);
    printf("Hmax=%g\n",Hmax);
    printf("DesNgb=%d\n",DesNgb);
    
    

    // allocate and build tree for gas particle //
    ngb3d_treeallocate(N_gas, 2*N_gas);
    set_particle_pointer(N_gas, x_gas,y_gas,z_gas, m_gas, rho_gas, u_gas);
    ngb3d_treebuild((float **)&P3d[1], N_gas, 0, dummy, dummy);
    h_guess = Hmax/150.0e0; h_guess_0=h_guess;
    
    // tree now constructed for gas, loop over star particles to find their hsml //
    for(i=0;i<N_star;i++)
    {
        xyz[0]=x_star[i];             //P3d[i+1]->Pos[0]+1.0e-10;
        xyz[1]=y_star[i];             //P3d[i+1]->Pos[1]+1.0e-10;
        xyz[2]=z_star[i];             //P3d[i+1]->Pos[2]+1.0e-10;
     //   printf("search position for i=%d (out of %d) = (%f|%f|%f)\n",i,N_star, xyz[0], xyz[1], xyz[2]);

	count = 0;
	maxcount = 30;
	do{
	  
            h2=ngb3d_treefind( xyz, DesNgb ,1.04*h_guess, &ngblist, &r2list, Hmax, &ngbfound);    
	    if(ngbfound < (DesNgb-2)) h_guess*=1.26;
            if(ngbfound > (DesNgb-2)) h_guess/=1.26;

	    count += 1;
	    //printf("count %d h_guess=%f ngbfound=%d\n",count, h_guess,ngbfound);
	}while( ( (ngbfound < (DesNgb-2)) || (ngbfound > (DesNgb+2) ) ) && (count < maxcount)  );

	//printf("ngbfound = %d\n", ngbfound);

        h = sqrt(h2);
        hinv= 1.0/h;
        hinv3= hinv * hinv * hinv;
        hinv4= hinv3 * hinv;
        
        Ngb = 0.0;
        Rho = 0.0;
        Temp= 0.0;
        
        for(n=0; n<ngbfound; n++)
        {
            j = ngblist[n];
            dx = xyz[0] - P3d[j+1]->Pos[0];
            dy = xyz[1] - P3d[j+1]->Pos[1];
            dz = xyz[2] - P3d[j+1]->Pos[2];
            //printf("n=%d \t j=%d p3d_id=%d \t  dr=(%f)\n\n\n",n,j,j, sqrt(dx*dx+dy*dy+dz*dz));
	    fflush(stdout);

            r2 = dx * dx + dy * dy + dz * dz;
        
        
            if(r2 < h2)
            {
                r = sqrt(r2);
                u = r * hinv;
                kernel_main(u, hinv3, hinv4, &wk, &dwk, 0);
                mj_wk = P3d[j+1]->mass * wk;
                
                Ngb         += wk;
                Rho         += mj_wk;
                Temp        += mj_wk * P3d[j+1]->u_gas ;
		//printf("mj_wk = %f\n",mj_wk);
            }
        }

	//printf("Rho in routine = %f\n",Rho);
        
            
        if(!(i%10000))
        {
            printf("i=%d hmax=%g h_guess=%g h=%g xyz=%g|%g|%g ngb=%d \n",
                   i,Hmax,h_guess,sqrt(h2),xyz[0],xyz[1],xyz[2],ngbfound); fflush(stdout);
        }
        H_OUT[i] = Temp / Rho;

        h_guess = sqrt(h2); // use this value for next guess, should speed things up //
        //if (h_guess>10.*h_guess_0) h_guess=2.*h_guess_0;
    } 
    
    ngb3d_treefree();
    free_memory_3d(N_gas);
    printf("done\n");
    return 0;
}



void set_particle_pointer(int Nsize, float* x, float* y, float* z, float* mass, float* rho_gas, float* u_gas)
{
    int i;
    float *pos;
    for(i=1;i<=Nsize;i++)
    {
        P3d[i]->Pos[0] = x[i-1];
        P3d[i]->Pos[1] = y[i-1];
        P3d[i]->Pos[2] = z[i-1];
//	printf(" xyz in tree is (%f|%f|%f)\n",x[i-1], y[i-1], z[i-1] );
        P3d[i]->mass   = mass[i-1];
        P3d[i]->rho_gas= rho_gas[i-1];
        P3d[i]->u_gas  = u_gas[i-1];

//float* rho_gas, float* u_gas

//a        P3d[i]->my_id  = i;
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

