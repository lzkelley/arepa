#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// revised call for python calling //
int kde(int N_in, float* x, float* y, float* z,  
        int N_search_min, int N_search_max,
        float rmin, int N_rgrid,
        float* R_GRID, float* N_R_GRID)
{
    double rmax = 0.5;
    double r2min = rmin*rmin, r2max=rmax*rmax;
    double log_rmin = log10( rmin );
    double log_rmax = log10( rmax );
    double dlogr = (log_rmax - log_rmin) / ((double)N_rgrid);
    double dlogr_inv = 1. / dlogr;
    double zeropt = -log_rmin * dlogr_inv;
    long i,j,k, Nmax=N_rgrid-1;
    double N_R_GRID_double[N_rgrid];
    for(i=0;i<N_rgrid;i++)
    {
        R_GRID[i] = log_rmin + dlogr*((double)i);
        N_R_GRID_double[i] = 0;
    }

    printf("N_search_min=%d N_search_max=%d N_in=%d \n",N_search_min,N_search_max,N_in);
    if(N_search_max>N_in) N_search_max=N_in;
    for(i=N_search_min;i<N_search_max-1;i++)
    {
        double x0=x[i], y0=y[i], z0=z[i];
        if(!(i%1000)) {printf("   ...i=%d \n",i); fflush(stdout);}
        for(j=i+1;j<N_in;j++)
        {
            double dx=x[j]-x0, dy=y[j]-y0, dz=z[j]-z0;
            if(dx<-0.5) {dx+=1;} else {if(dx>0.5) {dx-=1;}}
            if(dy<-0.5) {dy+=1;} else {if(dy>0.5) {dy-=1;}}
            if(dz<-0.5) {dz+=1;} else {if(dz>0.5) {dz-=1;}}
            double r2 = dx*dx + dy*dy + dz*dz;
            if(r2>r2max) continue;
            if(r2<r2min) continue;
            double log_r = 0.5*log10(r2);
            k = (int)(log_r*dlogr_inv + zeropt);
            if(k<0) {k=0;} else {if(k>Nmax) {k=Nmax;}}
            N_R_GRID_double[k]++;
        }
    }
    for(i=0;i<N_rgrid;i++)
    {
        N_R_GRID[i] = (float)N_R_GRID_double[i];
    }
    
    return 0;
}
