#include "nrutil.h"

static float BoxSizeX;
static float BoxSizeY;
static float BoxSizeZ;


int loadpositions(char *fname,int allocflag);


//float selectb(unsigned long k, unsigned long n, float arr[],int ind[]);


int set_particles_into_zones(void);

float cooling_rate(float dens,float press);

void treebuild(void);
int  find_neighbors(float xyz[3],float r,int *nl,float *r2l,int maxngb);



void begrun(void);
void endrun(int ierr);

void read_initial_conditions(char *fname);

/* in eos.c */
float cool(float enp1);
void   enttab(int irho, int itemp);

void   eospre(void);
void   eoscor(void);

void   eos2(int iact,int *iter2,int *iterflag);
void   eos3(int iact,int *iter2,int *iterflag);

void   settab(void);

void   getion(float Density,float EnergySpecific,
	      float *temp,float *decol,float *photo,float *MassFracIons);

void   xflux(void);

float phh(float xnu);
float phhe(float xnu);
float phhe2(float xnu);
float pchh(float xnu);
float pchhe(float xnu);
float pchhe2(float xnu);
float xinu(float xnu);
float xinunu(float xnu);
/****/


void increment_file_suffix(void);

void frduni(void);
float z2t(float z,float OmegaMatter,float OmegaLambda,float HubbleToday);
float t2a(float Time,float OmegaMatter,float OmegaLambda,float HubbleToday,float TimeToday);
float ht0(float eta);


void glbout(void);

void init(void);
void read_parameter_file(char *fname);

void forwrd(void);

void bckwrd(void);

void oldnew(void);

void pltout(void);

void prnout(void);

void propag(void);



void restart(int mod);



/*** system.c ***/
float second(void);
float tremain(void);
float dmin(double,double);
float dmax(double,double);
int imin(int,int);
int imax(int,int);
int    nint(double);
/****/

void densty(int iact, int nsph, int *isph);
void trweos(void);


void grvfrc(void);
void trwfrc(void);



void trwngh(void);


void tstep(int startflag);







float zbrentMS(float (*func)(double), float x1, float x2, float tol,int *ierr);
float midinf(float (*funk)(double), float aa, float bb, int n);
void polint(float xa[], float ya[], int n, float x, float *y, float *dy);
float zriddr(float (*func)(double), float x1, float x2, float xacc);




void indexintx(unsigned long n, int arr[], unsigned long indx[]);
void indexshortx(unsigned long n, short arr[], int indx[]);
void indexx(unsigned long n, float arr[], int indx[]);
void rank(unsigned long n, int indx[], int irank[]);



float selip(unsigned long k, unsigned long n, float arr[]);


void ngb3d_treebuild(float **pospointer,int Npart,int MinMaxFlag,float XminOpt[3],float XmaxOpt[3]);
void ngb3d_treefree(void);
void ngb3d_treeallocate(int npart,int maxnodes);  /* usually maxnodes=2*npart is suffiecient */
float ngb3d_treefind(float xyz[3], int desngb, int desngb_tolerance, float hguess,int **ngblistback, float **r2listback,float hmin,float hmax,int *ngbfound);
float ngb3d_treetest(float xyz[3], int desngb, float hguess,int **ngblistback, float **r2listback);
void allocate_ngblists(void);
void free_ngblists(void);

void ngb_treebuild(float **pospointer,int Npart,int MinMaxFlag,float XminOpt[3],float XmaxOpt[3]);

void ngb_treefree(void);
void ngb_treeallocate(int npart,int maxnodes);  /* usually maxnodes=2*npart is suffiecient */

float ngb_treefind(float xyz[3], int desngb, float hguess,int **ngblistback, float **r2listback);


float ngb_treetest(float xyz[3], int desngb, float hguess,int **ngblistback, float **r2listback);

void force_treeallocate(int maxnodes);  /* usually maxnodes=2*npart is suffiecient */
void force_treefree(void);


void force_treebuild(float **pospointer,int Npart,
		     float thetamax,
		     float maxSoftening,
		     float *softeningtable);



void force_treeevaluate(float *targetpart);


void force_testforce(float *targetpart);
