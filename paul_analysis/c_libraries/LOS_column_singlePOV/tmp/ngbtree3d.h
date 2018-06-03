

void ngb3d_treebuild(float **pospointer,long Npart,long MinMaxFlag,float XminOpt[3],float XmaxOpt[3]);

void ngb3d_treefree(void);
void ngb3d_treeallocate(long npart,long maxnodes);  /* usually maxnodes=2*npart is suffiecient */

float ngb3d_treefind(float xyz[3], long desngb, float hguess,long **ngblistback, float **r2listback,float hmax,long *ngbfound);


float ngb3d_treetest(float xyz[3], long desngb, float hguess,long **ngblistback, float **r2listback);
