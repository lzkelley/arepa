#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "proto.h"

static struct NODE
{
    float  center[3],len;  /* center and sidelength of treecubes */
    float  xmin[3],xmax[3];
    int    count;          /* total # of particles in cell      */
    struct NODE *father,*suns[8];
    int    first;          /* index of first particle           */
} *nodes;

static struct COORDINATES
{
    float xyz[3];
} **Pos;

static int    numnodes,MaxNodes;
static int    N;
static int    *next,*previous;   /* Link-List for particle groups */
static int    *ngblist,numngb;
static float  *r2list;
static float searchmin[3],searchmax[3];



float ngb3d_treefind(float xyz[3], int desngb, int desngb_tolerance, float hguess,int **ngblistback, float **r2listback,
                     float hmin, float hmax, int *ngbfound)
{
#define SECFACTOR  1.2
#define  PI        3.1415927
    void ngb3d_treesearch(struct NODE *this, float xyz[3]);
    float selectb(unsigned long k, unsigned long n, float arr[],int ind[]);
    float sr,sr2,h2max;  /* search radius */
    int  i,ind,ni,j,subnode,fak,k,rep=0;
    float dx,dy,dz,r2;
    struct NODE *th,*nn;
    
    
    
    if(hguess>0)
    {
        sr=hguess;
        if(sr>hmax) sr=hmax;
        if(sr<hmin) sr=hmin;
    }
    else
    {
        /* determine estimate of local density */
        th=nodes;
        while(th->count>200)
        {
            for(j=0,subnode=0,fak=1;j<3;j++,fak<<=1)
                if(xyz[j]>th->center[j])
                    subnode+=fak;
            
            if(nn=th->suns[subnode])
                if(nn->count>200)
                    th=nn;
                else
                    break;
                else
                    break;
        }
        
        sr=th->len * pow((3.0/(4*PI)*SECFACTOR)*desngb/((float)(th->count)),1.0/3);
        if(sr>hmax) sr=hmax;
        if(sr<hmin) sr=hmin;
    }
    
    double left = 0.0;
    double right = 1.0e37;
    double fac = 1;
    do
    {
        for(k=0;k<3;k++)
        {
            searchmin[k]=-sr;
            searchmax[k]=+sr;
        }
        sr2=sr*sr;
        
        numngb=0;
        ngb3d_treesearch(nodes, xyz);  rep++;
        
        if(rep > 300)
        {
            printf("...long iteration: numngb = %d sr = %g hmax = %g hmin = %g desnum = %d desnumtol = %d rep = %d xyz=%g/%g/%g \n",
                   numngb,sr,hmax,hmin,desngb,desngb_tolerance,rep,xyz[0],xyz[1],xyz[2]);
        }
        
        if((left<right) && (fabs(right-left) > 0.001*(right+left)) && (rep < 10000))
        {
            if(numngb < desngb - desngb_tolerance)
            {
                if(sr>left) left=sr;
                if(sr<hmax)
                {
                    if(numngb>1)
                    {
                        //sr*=pow((2.1*(float)desngb)/numngb,1.0/3.);
                        fac = pow(((float)desngb)/numngb,1.0/3.0);
                        double flim = 1.+ 0.4 / pow((1.+(float)rep),0.33);
                        if(fac > flim) fac=flim;
                    }
                    else
                    {
                        fac = 1. + 1.0 / pow((1.+(float)rep),0.33);
                    }
                    sr *= fac;
                    if(sr <= left) sr = left * (1. + 0.05/pow((1.+(float)rep),0.25));
                    if(sr >= right) sr = right / (1. + 0.05/pow((1.+(float)rep),0.25));
                    if(sr>hmax) sr=hmax;
                    if(sr<hmin) sr=hmin;
                    continue;
                }
                else
                {
                    for(i=0;i<numngb;i++)
                    {
                        ind=ngblist[i];
                        dx=Pos[ind]->xyz[0]-xyz[0];
                        if(dx > BoxSizeX/2) dx-=BoxSizeX;
                        if(dx < -BoxSizeX/2) dx+=BoxSizeX;
                        dy=Pos[ind]->xyz[1]-xyz[1];
                        if(dy > BoxSizeY/2) dy-=BoxSizeY;
                        if(dy < -BoxSizeY/2) dy+=BoxSizeY;
                        dz=Pos[ind]->xyz[2]-xyz[2];
                        if(dz > BoxSizeZ/2) dz-=BoxSizeZ;
                        if(dz < -BoxSizeZ/2) dz+=BoxSizeZ;
                        r2=dx*dx+dy*dy+dz*dz;
                        r2list[i]=r2;
                    }
                    *ngbfound=numngb;
                    *ngblistback=ngblist;
                    *r2listback=r2list;
                    
                    return hmax*hmax;
                }
                
            }
            
            
            if(numngb > desngb + desngb_tolerance)
            {
                if(sr<right) right=sr;
                if(sr>hmin)
                {
                    if(numngb>1)
                    {
                        //sr*=pow(((1./2.1)*(float)desngb)/numngb,1.0/3);
                        fac = pow(((float)desngb)/numngb,1.0/3.0);
                        double flim = 1. - 0.3 / pow((1.+(float)rep),0.33);
                        if(fac < flim) fac=flim;
                    }
                    else
                    {
                        fac = 1. / (1. + 1.0 / pow((1.+(float)rep),0.33));
                    }
                    sr *= fac;
                    if(sr <= left) sr = left * (1. + 0.05/pow((1.+(float)rep),0.25));
                    if(sr >= right) sr = right / (1. + 0.05/pow((1.+(float)rep),0.25));
                    if(sr>hmax) sr=hmax;
                    if(sr<hmin) sr=hmin;
                    continue;
                }
                else
                {
                    for(i=0;i<numngb;i++)
                    {
                        ind=ngblist[i];
                        dx=Pos[ind]->xyz[0]-xyz[0];
                        if(dx > BoxSizeX/2) dx-=BoxSizeX;
                        if(dx < -BoxSizeX/2) dx+=BoxSizeX;
                        dy=Pos[ind]->xyz[1]-xyz[1];
                        if(dy > BoxSizeY/2) dy-=BoxSizeY;
                        if(dy < -BoxSizeY/2) dy+=BoxSizeY;
                        dz=Pos[ind]->xyz[2]-xyz[2];
                        if(dz > BoxSizeZ/2) dz-=BoxSizeZ;
                        if(dz < -BoxSizeZ/2) dz+=BoxSizeZ;
                        r2=dx*dx+dy*dy+dz*dz;
                        r2list[i]=r2;
                    }
                    *ngbfound=numngb;
                    *ngblistback=ngblist;
                    *r2listback=r2list;
                    
                    return hmin*hmin;
                }
                
            }
        }
        
        
        for(i=0;i<numngb;i++)
        {
            ind=ngblist[i];
            dx=Pos[ind]->xyz[0]-xyz[0];
            if(dx > BoxSizeX/2) dx-=BoxSizeX;
            if(dx < -BoxSizeX/2) dx+=BoxSizeX;
            dy=Pos[ind]->xyz[1]-xyz[1];
            if(dy > BoxSizeY/2) dy-=BoxSizeY;
            if(dy < -BoxSizeY/2) dy+=BoxSizeY;
            dz=Pos[ind]->xyz[2]-xyz[2];
            if(dz > BoxSizeZ/2) dz-=BoxSizeZ;
            if(dz < -BoxSizeZ/2) dz+=BoxSizeZ;
            r2=dx*dx+dy*dy+dz*dz;
            
            r2list[i]=r2;
            //printf("r2 = %g P_xyz=%g/%g/%g 00=%g/%g/%g \n",r2,Pos[ind]->xyz[0],Pos[ind]->xyz[1],Pos[ind]->xyz[2],xyz[0],xyz[1],xyz[2]);
        }

        h2max = sr2;
        //h2max=0; for(i=0;i<numngb;i++) {if(r2list[i]>h2max) h2max=r2list[i];}
        //h2max=selectb(desngb,numngb,r2list-1,ngblist-1);
        //if(h2max<=sr2) break;
        break;
    }
    while(1);
    
    
    /*  printf("numngb: %d   rep: %d\n",numngb,rep);  */
    *ngbfound=numngb;
    *ngblistback=ngblist;
    *r2listback=r2list;
    
    return h2max;
}



void ngb3d_treesearch(struct NODE *this, float xyz[3])
{
    int k,p;
    struct NODE *nn;
    float BoxSize[3];
    BoxSize[0]=BoxSizeX;
    BoxSize[1]=BoxSizeY;
    BoxSize[2]=BoxSizeZ;
    float dd;
    
    if(this->count==1)
    {
        for(k=0,p=this->first;k<3;k++)
        {
            dd = Pos[p]->xyz[k] - xyz[k];
            if(dd> BoxSize[k]/2) dd-=BoxSize[k];
            if(dd<-BoxSize[k]/2) dd+=BoxSize[k];
            if(dd > searchmax[k]) return;
            if(dd < searchmin[k]) return;
        }
        ngblist[numngb++]=p;
    }
    else
    {
        for(k=0;k<3;k++)
        {
            dd = this->xmax[k] - xyz[k];
            if(dd> BoxSize[k]/2) dd-=BoxSize[k];
            if(dd<-BoxSize[k]/2) dd+=BoxSize[k];
            if(dd < searchmin[k]) return;

            dd = this->xmin[k] - xyz[k];
            if(dd> BoxSize[k]/2) dd-=BoxSize[k];
            if(dd<-BoxSize[k]/2) dd+=BoxSize[k];
            if(dd > searchmax[k]) return;
        }
        
        for(k=0;k<3;k++)
        {
            dd = this->xmax[k] - xyz[k];
            if(dd> BoxSize[k]/2) dd-=BoxSize[k];
            if(dd<-BoxSize[k]/2) dd+=BoxSize[k];
            if(dd > searchmax[k]) break;
            
            dd = this->xmin[k] - xyz[k];
            if(dd> BoxSize[k]/2) dd-=BoxSize[k];
            if(dd<-BoxSize[k]/2) dd+=BoxSize[k];
            if(dd < searchmin[k]) break;
        }
        
        if(k>=3)
        {
            /* cell lies completely inside */
            
            p=this->first;
            
            for(k=0;k<this->count;k++)
            {
                ngblist[numngb++]=p;
                p=next[p];
            }
        }
        else
        {
            for(k=0;k<8;k++)
                if(nn=this->suns[k])
                {
                    ngb3d_treesearch(nn, xyz);
                }
        }
    }
}




void ngb3d_treeallocate(int npart,int maxnodes)  /* usually maxnodes=2*npart is suffiecient */
{
    MaxNodes=maxnodes;
    N=npart;
    
    if(!(nodes=malloc(MaxNodes*sizeof(struct NODE))))
    {
        fprintf(stderr,"Failed to allocate %d nodes.\n",MaxNodes);
        exit(0);
    }
    
    if(!(next=malloc(N*sizeof(int))))
    {
        fprintf(stderr,"Failed to allocate %d spaces for next array\n",N);
        exit(0);
    }
    
    if(!(previous=malloc(N*sizeof(int))))
    {
        fprintf(stderr,"Failed to allocate %d spaces for previous array\n",N);
        exit(0);
    }
    
    if(!(ngblist=malloc(N*sizeof(int))))
    {
        fprintf(stderr,"Failed to allocate %d spaces for ngblist array\n",N);
        exit(0);
    }
    
    if(!(r2list=malloc(N*sizeof(float))))
    {
        fprintf(stderr,"Failed to allocate %d spaces for r2list array\n",N);
        exit(0);
    }
}




void ngb3d_treefree(void)
{
    free(nodes);
    free(next);
    free(previous);
    free(ngblist);
    free(r2list);
}










void ngb3d_treebuild(float **pospointer,int Npart,int MinMaxFlag,float XminOpt[3],float XmaxOpt[3])
/* packs the particles 0...Npart-1 in tree */
{
    int i,j,k,subp,subi,p,ni,subnode,fak;
    float xmin[3],xmax[3],len,x;
    
    struct NODE *nfree,*th,*nn;
    
    
    printf("Begin tree construction.\n");
    if(Npart<2)
    {
        fprintf(stderr,"must be at least two particles in tree.\n");
        exit(0);
    }
    
    Pos=(struct COORDINATES **)pospointer;
    N=Npart;
    
    
    if(MinMaxFlag)
    {
        for(j=0;j<3;j++)
        {
            xmin[j]=XminOpt[j];
            xmax[j]=XmaxOpt[j];
        }
    }
    else
    {
        for(j=0;j<3;j++)
            xmin[j]=xmax[j]=Pos[0]->xyz[j];
    }
    
    /* find enclosing rectangle */
    for(i=1;i<Npart;i++)
        for(j=0;j<3;j++)
        {
            if(Pos[i]->xyz[j]>xmax[j])
                xmax[j]=Pos[i]->xyz[j];
            if(Pos[i]->xyz[j]<xmin[j])
                xmin[j]=Pos[i]->xyz[j];
        }
    
    /* determine maxmimum extension */
    for(j=1,len=xmax[0]-xmin[0];j<3;j++)
        if((xmax[j]-xmin[j])>len)
            len=xmax[j]-xmin[j];
    
    len*=1.01;
    
    printf("Enclosing max has length= %g\n",len); fflush(stdout);
    
    
    /* insert particle 1 in root node */
    
    nfree=nodes;
    
    for(j=0;j<3;j++)
        nfree->center[j]=(xmax[j]+xmin[j])/2;
    nfree->len=len;
    
    nfree->father=0;
    for(i=0;i<8;i++)
        nfree->suns[i]=0;
    nfree->first=0;
    
    nfree->count=1;
    
    numnodes=1;
    nfree++;
    
    
    next[0]=-1;
    previous[0]=-1;
    
    
    for(i=1;i<Npart;i++)  /* insert all other particles */
    {
        th=nodes;
        
        /*
         if(!(i%10000))
         {
         printf("%d..",i); fflush(stdout);
         }
         */
        
        while(1)
        {
            th->count++;           /* we're putting a particle in, so increase the count.
                                    note that every branch in get's it's count increase
                                    as we walk into a level */
            
            if(th->count==2)       /* cell was occupied with only one particle */
                break;
            
            for(j=0,subnode=0,fak=1;j<3;j++,fak<<=1)   /* which sector do we go into */
                if(Pos[i]->xyz[j]>th->center[j])
                    subnode+=fak;
            
            if(nn=th->suns[subnode])                  /* if it exists, step into it */
                th=nn;
            else
                break;                                  /* else create it */
        }
        
        
        if(th->count==2)  /* cell was occcupied with one particle */
        {
            while(1)
            {
                p=th->first;
                
                for(j=0,subp=0,fak=1;j<3;j++,fak<<=1)    /* which sector is the original particle */
                    if(Pos[p]->xyz[j]>th->center[j])
                        subp+=fak;
                
                nfree->father=th;                        /* father = this node */
                for(j=0;j<8;j++)
                    nfree->suns[j]=0;                      /* no sons yet */
                
                nfree->len=th->len/2;
                
                for(j=0;j<3;j++)
                    nfree->center[j]=th->center[j];
                
                for(j=0;j<3;j++)
                    if(Pos[p]->xyz[j]>nfree->center[j])
                        nfree->center[j]+=nfree->len/2;
                    else
                        nfree->center[j]-=nfree->len/2;
                
                nfree->first=p;
                nfree->count=1;
                
                th->suns[subp]=nfree;
                
                numnodes++;
                nfree++;
                
                if(numnodes>=MaxNodes)
                {
                    fprintf(stderr,"maximum node number %d in neighbour tree reached.\n",numnodes);
                    fprintf(stderr,"subi= %d  subp= %d\n",subi,subp);
                    fprintf(stderr,"i=%d  xyz= %g  %g  %g\n",i,Pos[i]->xyz[0],Pos[i]->xyz[1],Pos[i]->xyz[2]);
                    fprintf(stderr,"p=%d  xyz= %g  %g  %g\n",p,Pos[p]->xyz[0],Pos[p]->xyz[1],Pos[p]->xyz[2]);
                    fprintf(stderr,"parent len= %g center = %g %g %g\n",th->len,th->center[0],th->center[1],th->center[2]);
                    exit(0);
                }
                
                for(j=0,subi=0,fak=1;j<3;j++,fak<<=1)   /* which sector is the new particle */
                    if(Pos[i]->xyz[j]>th->center[j])
                        subi+=fak;
                
                /*  Add this to randomize sub-node if the particles are approximately on top of each
                 *  other.  If we don't do this then we tend to get caught in a loop where it 
                 *  continually decreases the len, and finally runs out of nodes.*/
                if(th->len < 1.0e-5)
                {
                    subi = (int) (4.0 * rand() / (RAND_MAX + 1.0));
                }
                
                if(subi==subp)
                {
                    th=nfree-1;
                    th->count++;
                }
                else
                    break;
            }
            
        }
        
        
        
        for(j=0,subi=0,fak=1;j<3;j++,fak<<=1)
            if(Pos[i]->xyz[j]>th->center[j])
                subi+=fak;
        
        nfree->father=th;
        
        
        p=th->first;
        for(j=0;j<(th->count-2);j++)
        {
            p=next[p];
            
            /*	  if(p<0)
             {
             printf("error\n");
             exit(0);
             } */
        }
        
        if(next[p]>=0)
        {
            previous[next[p]]=i;
        }
        
        next[i]=next[p];
        previous[i]=p;
        next[p]=i;
        
        
        for(j=0;j<8;j++)
            nfree->suns[j]=0;
        
        nfree->len=th->len/2;
        for(j=0;j<3;j++)
            nfree->center[j]=th->center[j];
        
        for(j=0;j<3;j++)
            if(Pos[i]->xyz[j]>nfree->center[j])
                nfree->center[j]+=nfree->len/2;
            else
                nfree->center[j]-=nfree->len/2;
        
        nfree->count=1;
        
        nfree->first=i;
        th->suns[subi]=nfree;
        
        numnodes++;
        nfree++;
        
        if(numnodes>=MaxNodes)
        {
            fprintf(stderr,"maximum node number %d in neighbour tree reached.\n",numnodes);
            fprintf(stderr,"i=%d  xyz= %g  %g  %g\n",i,Pos[i]->xyz[0],Pos[i]->xyz[1],Pos[i]->xyz[2]);
            exit(0);
        }
    }
    
    
    
    for(ni=0,th=nodes;ni<numnodes;ni++,th++)
    {
        for(k=0;k<3;k++)
        {
            th->xmin[k]=th->center[k]-th->len/2;
            th->xmax[k]=th->center[k]+th->len/2;
        }
    }
    
    
    printf("Tree contruction finished. Number of nodes: %d\n",numnodes);
}


