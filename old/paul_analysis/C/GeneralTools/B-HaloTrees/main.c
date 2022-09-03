#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <signal.h>
#include <hdf5.h>

#include "allvars.h"
#include "proto.h"



int CountUsed, CountSumUsed;

int main(int argc, char **argv)
{
  int num, count;

  if(argc != 2)
    {
      printf("\n  usage: B-HaloTrees <parameterfile>\n");
      printf("  <parameterfile>    see readparmeterfile.\n\n");
      exit(1);
    }

  read_parameter_file(argv[1]);
#ifdef SKIP_SNAP
  read_outputlist(OutputList);
#endif

  Cats = mymalloc(sizeof(struct halo_catalogue) * (LastSnapShotNr + 1));
  FirstHaloInSnap = mymalloc(sizeof(int) * (LastSnapShotNr + 1));

  count_halos();

  printf("TotHalos = %d\n", TotHalos);
  printf("Want to allocate %g GB\n",
	 ((double) TotHalos) * (sizeof(struct halo_data) +
				sizeof(struct halo_aux_data)) / (1024.0 * 1024.0 * 1024.0));

  Halo    = mymalloc(TotHalos * sizeof(struct halo_data));
  HaloAux = mymalloc(TotHalos * sizeof(struct halo_aux_data));


  printf("loading halo catalogues and descendant tree files...\n");
  fflush(stdout);

  for(num = LastSnapShotNr, count = 0; num >= FirstSnapShotNr; num -= SnapSkipFac)
    {
#ifdef SKIP_SNAP
      if(OutputListFlag[num] ==0)
        {
          continue;
        }
#endif

      FirstHaloInSnap[num] = count;
      load_subhalo_catalogue(num);
      count += Cats[num].TotNsubhalos;

      read_snap_header_attributes_in_hdf5(num);
    }

  printf("done.\n");


  printf("setting progenitor pointers...\n");
  fflush(stdout);
  set_progenitor_pointers();
  printf("progenitor pointers done.\n");
  fflush(stdout);

  generate_trees();

  return 0;
}



void generate_trees(void)
{
  int filenr, i, k, maxhalos, treenr;
  int *npertree;
  char buf[500], hdf5_buf[500];


  for(i = 0; i < TotHalos; i++)
    HaloAux[i].UsedFlag = HaloAux[i].HaloFlag = 0;
  
  for(i = 0, filenr = 0; i < Cats[LastSnapShotNr].TotNsubhalos; i++)
    {
      if(HaloAux[i].UsedFlag == 0)
	{
	  walk_it(i, 0, filenr);
	  
	  filenr++;
	  
	  if(filenr >= NumberOfOutputFiles)
	    filenr = 0;
	}
    }

 

  for(filenr=0; filenr < NumberOfOutputFiles; filenr++)
    {

      NtreesPerFile = 0;
      NhalosPerFile = 0;
      
      CountSumUsed = 0;
      
      for(i = 0; i < TotHalos; i++)
	HaloAux[i].UsedFlag = HaloAux[i].HaloFlag = 0;
      
      maxhalos = 0;
      
      for(i = 0; i < Cats[LastSnapShotNr].TotNsubhalos; i++)
	{
	  if(HaloAux[i].UsedFlag == 0 && HaloAux[i].FileNr == filenr)
	    {
	      CountUsed = 0;
	      
	      walk_it(i, 0, filenr);
	      
	      NtreesPerFile += 1;
	      NhalosPerFile += CountUsed;
	      
	      if(CountUsed > maxhalos)
		maxhalos = CountUsed;
	      
	      CountSumUsed += CountUsed;
	    }
	}

      printf("TotHalos=%d   Used=%d   maxhalos=%d\n", TotHalos, CountSumUsed, maxhalos);
      fflush(stdout);
      
      for(i = 0; i < TotHalos; i++)
	HaloAux[i].UsedFlag = HaloAux[i].HaloFlag = 0;
      

      HaloList = mymalloc(maxhalos * sizeof(struct halo_data));
      
      
      sprintf(buf, "%s/treedata/trees_sf%d_%03d.%d", TreeOutputDir, SnapSkipFac, LastSnapShotNr, filenr);
      
      printf("starting: %s\n", buf);
      fflush(stdout);
      

/*==============================================================================================*/
      hid_t hdf5_file = 0, hdf5_headergrp = 0;
      sprintf(hdf5_buf, "%s/treedata/trees_sf%d_%03d.%d.hdf5", TreeOutputDir, SnapSkipFac, LastSnapShotNr, filenr);
      printf("and : %s \n", hdf5_buf);
      fflush(stdout);

      hdf5_file = H5Fcreate(hdf5_buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
/*==============================================================================================*/

      printf("NtreesPerFile=%d  NhalosPerFile=%d\n", NtreesPerFile, NhalosPerFile);
      
      npertree = mymalloc(NtreesPerFile * sizeof(int));
      for(i = 0; i < NtreesPerFile; i++)
	npertree[i] = 0;
      
      treenr = 0;
      
      for(i = 0; i < Cats[LastSnapShotNr].TotNsubhalos; i++)
	{
	  if(HaloAux[i].UsedFlag == 0 && HaloAux[i].FileNr == filenr)
	    {
	      CountUsed = 0;
	      
	      walk_it(i, 1, filenr);
	      
	      for(k = 0; k < CountUsed; k++)
		{
		  if(HaloList[k].Descendant >= 0)
		    HaloList[k].Descendant = HaloAux[HaloList[k].Descendant].TargetIndex;
		  
		  if(HaloList[k].FirstProgenitor >= 0)
		    HaloList[k].FirstProgenitor = HaloAux[HaloList[k].FirstProgenitor].TargetIndex;
		  
		  if(HaloList[k].NextProgenitor >= 0)
		    HaloList[k].NextProgenitor = HaloAux[HaloList[k].NextProgenitor].TargetIndex;
		  
		  if(HaloList[k].FirstHaloInFOFgroup >= 0)
		    HaloList[k].FirstHaloInFOFgroup = HaloAux[HaloList[k].FirstHaloInFOFgroup].TargetIndex;
		  
		  if(HaloList[k].NextHaloInFOFgroup >= 0)
		    HaloList[k].NextHaloInFOFgroup = HaloAux[HaloList[k].NextHaloInFOFgroup].TargetIndex;
		}

	      if(CountUsed > -1) 
		{
		  printf("writing tree treenr=%d and i=%d with %d subhaloes\n",treenr,i,CountUsed);
		  fflush(stdout);
	          my_hdf5_write_tree(HaloList, treenr, CountUsed, hdf5_file);			/* new routine to write TreeData[treenr] */
		}	      

	      
	      npertree[treenr] = CountUsed;
	      treenr++;
	    }
	  
	}
     

/*==============================================================================================*/
      hdf5_headergrp = H5Gcreate(hdf5_file, "/Header", 0);
      write_header_attributes_in_hdf5(hdf5_headergrp, npertree);
      H5Gclose(hdf5_headergrp);
/*==============================================================================================*/

 
      H5Fclose(hdf5_file);
      

      myfree(npertree);
      
      
      myfree(HaloList);
      
      printf("Saved=%d\n", CountSumUsed);

    }
}






void count_halos(void)
{
  int num, ngroups, nids, nFiles, nsubhalos;
  long long totNids;
  char buf[1000];

  TotHalos = 0;

  for(num = LastSnapShotNr; num >= FirstSnapShotNr; num -= SnapSkipFac)
    {
#ifdef SKIP_SNAP
      if(OutputListFlag[num] ==0)
        {
          continue;
        }
#endif

      nFiles = 1;

      get_filename(buf, num, 0, 1);
      read_hdf5_subhalo_header(buf, num, &Cats[num], &nFiles, &totNids, &nids, &nsubhalos, &ngroups);

      TotHalos += Cats[num].TotNsubhalos;
    }
  printf("total number of halos=%d\n", TotHalos);
}





void load_subhalo_catalogue(int num)
{
  int i, ngroups, nids, nFiles, nsubhalos, subcount;
  int groupcount, filenr, ncount;
  int subgr, subnr, gr, nh, sc, gr_nh;
  char buf[1000];
  FILE *fd;
  int *groupNsubs, *descendant_haloindex, *descendant_snapnum, *filenrOfHalo, *subhaloindex;
  float *group_M_Mean200, *group_M_Crit200, *group_M_TopHat200;
  long long totNids; 
  int type;

  int *groupLenType;
  int *subhaloGrNr, *subhaloSnapFileNr, *subhaloLen, *subhaloLenType;
  float *subhaloPos, *subhaloVel, *subhaloSpin, *subhaloVelDisp, *subhaloVmax, *subhaloHalfmassRad, *subhaloMassType, *subhaloMassInRadType, *subhaloHalfmassRadType;
  MyIDType *subhaloIDMostbound, *groupOffsetType, *subhaloOffsetType;

#ifdef SUBFIND_EXTRA
  float *subhaloGasMetalFractions, *subhaloGasMetalFractionsSfr, *subhaloGasMetallicity, *subhaloGasMetallicitySfr, *subhaloSFR, *subhaloStarMetalFractions, *subhaloStarMetallicity, *subhaloBHMass, *subhaloBHMdot, *subhaloSFRinRad, *subhaloStellarPhotometrics;
#endif

  printf("Catalogue num=%d has %d groups and %d subhalos \n", num,Cats[num].TotNgroups, Cats[num].TotNsubhalos);

/*==============================================================*/
/*          ALLOCATE ARRAYS TO READ SUBFIND CATALOGS 		*/
/*==============================================================*/
  descendant_haloindex  = mymalloc(1 * sizeof(int) * Cats[num].TotNsubhalos);
  descendant_snapnum    = mymalloc(1 * sizeof(int) * Cats[num].TotNsubhalos);
  filenrOfHalo          = mymalloc(1 * sizeof(int) * Cats[num].TotNsubhalos);
  subhaloindex          = mymalloc(1 * sizeof(int) * Cats[num].TotNsubhalos);

  subhaloSnapFileNr     = mymalloc(1 * sizeof(int) * Cats[num].TotNsubhalos);
  groupOffsetType       = mymalloc(6 * sizeof(MyIDType) * Cats[num].TotNgroups);
  subhaloOffsetType     = mymalloc(6 * sizeof(MyIDType) * Cats[num].TotNsubhalos);

  groupNsubs 		= mymalloc(1 * sizeof(int) * Cats[num].TotNgroups);
  group_M_Mean200	= mymalloc(1 * sizeof(float) * Cats[num].TotNgroups);
  group_M_Crit200	= mymalloc(1 * sizeof(float) * Cats[num].TotNgroups);
  group_M_TopHat200 	= mymalloc(1 * sizeof(float) * Cats[num].TotNgroups);
  groupLenType       	= mymalloc(6 * sizeof(int) * Cats[num].TotNgroups);

  subhaloLen 	          = mymalloc(1 * sizeof(int) * Cats[num].TotNsubhalos);
  subhaloLenType          = mymalloc(6 * sizeof(int) * Cats[num].TotNsubhalos);
  subhaloPos 	          = mymalloc(3 * sizeof(float) * Cats[num].TotNsubhalos);
  subhaloVel 	          = mymalloc(3 * sizeof(float) * Cats[num].TotNsubhalos);
  subhaloSpin             = mymalloc(3 * sizeof(float) * Cats[num].TotNsubhalos);
  subhaloVelDisp          = mymalloc(1 * sizeof(float) * Cats[num].TotNsubhalos);
  subhaloVmax 	          = mymalloc(1 * sizeof(float) * Cats[num].TotNsubhalos);
  subhaloHalfmassRad      = mymalloc(1 * sizeof(float) * Cats[num].TotNsubhalos);
  subhaloIDMostbound      = mymalloc(1 * sizeof(MyIDType) * Cats[num].TotNsubhalos);
  subhaloMassType         = mymalloc(6 * sizeof(float) * Cats[num].TotNsubhalos);
  subhaloMassInRadType    = mymalloc(6 * sizeof(float) * Cats[num].TotNsubhalos);
  subhaloHalfmassRadType  = mymalloc(6 * sizeof(float) * Cats[num].TotNsubhalos);
  subhaloGrNr             = mymalloc(1 * sizeof(MyIDType) * Cats[num].TotNsubhalos);

#ifdef SUBFIND_EXTRA
  subhaloGasMetallicity       = mymalloc(1 * sizeof(float) * Cats[num].TotNsubhalos);
  subhaloGasMetalFractions    = mymalloc(9 * sizeof(float) * Cats[num].TotNsubhalos);
  subhaloGasMetallicitySfr    = mymalloc(1 * sizeof(float) * Cats[num].TotNsubhalos);
  subhaloGasMetalFractionsSfr = mymalloc(9 * sizeof(float) * Cats[num].TotNsubhalos);
  subhaloStarMetallicity      = mymalloc(1 * sizeof(float) * Cats[num].TotNsubhalos);
  subhaloStarMetalFractions   = mymalloc(9 * sizeof(float) * Cats[num].TotNsubhalos);
  subhaloSFR		      = mymalloc(1 * sizeof(float) * Cats[num].TotNsubhalos);
  subhaloBHMass               = mymalloc(1 * sizeof(float) * Cats[num].TotNsubhalos);
  subhaloBHMdot               = mymalloc(1 * sizeof(float) * Cats[num].TotNsubhalos);
  subhaloSFRinRad             = mymalloc(1 * sizeof(float) * Cats[num].TotNsubhalos);
  subhaloStellarPhotometrics  = mymalloc(8 * sizeof(float) * Cats[num].TotNsubhalos);
#endif

  subcount = 0;
  groupcount = 0;

  nFiles = 1;

  for(filenr = 0; filenr < nFiles; filenr++)
    {
      get_filename(buf, num, filenr, 1);    
      printf("%s\n",buf);
      fflush(stdout);

      read_hdf5_subhalo_header(buf, num, &Cats[num], &nFiles, &totNids, &nids, &nsubhalos, &ngroups);

      if (ngroups > 0)
	{
          read_hdf5_int_dataset(buf,"/Group","GroupNsubs"		,ngroups, 1, &groupNsubs[groupcount]);
          read_hdf5_float_dataset(buf,"/Group","Group_M_Mean200"	,ngroups, 1, &group_M_Mean200[groupcount]);
          read_hdf5_float_dataset(buf,"/Group","Group_M_Crit200"	,ngroups, 1, &group_M_Crit200[groupcount]);
          read_hdf5_float_dataset(buf,"/Group","Group_M_TopHat200"	,ngroups, 1, &group_M_TopHat200[groupcount]);
          read_hdf5_int_dataset(buf,"/Group","GroupLenType"		,ngroups, 6, &groupLenType[6 * groupcount]);
	}

      if (nsubhalos > 0)
	{
          read_hdf5_int_dataset(buf,"/Subhalo","SubhaloLen"	 	  ,nsubhalos, 1, &subhaloLen[subcount]);
          read_hdf5_int_dataset(buf,"/Subhalo","SubhaloLenType"		  ,nsubhalos, 6, &subhaloLenType[6 * subcount]);
          read_hdf5_float_dataset(buf,"/Subhalo","SubhaloPos"		  ,nsubhalos, 3, &subhaloPos[3 * subcount]);
          read_hdf5_float_dataset(buf,"/Subhalo","SubhaloVel"		  ,nsubhalos, 3, &subhaloVel[3 * subcount]);
          read_hdf5_float_dataset(buf,"/Subhalo","SubhaloSpin"		  ,nsubhalos, 3, &subhaloSpin[3 * subcount]);
          read_hdf5_float_dataset(buf,"/Subhalo","SubhaloVelDisp"	  ,nsubhalos, 1, &subhaloVelDisp[subcount]);
          read_hdf5_float_dataset(buf,"/Subhalo","SubhaloVmax"		  ,nsubhalos, 1, &subhaloVmax[subcount]);	
          read_hdf5_float_dataset(buf,"/Subhalo","SubhaloHalfmassRad"	  ,nsubhalos, 1, &subhaloHalfmassRad[subcount]);
          read_hdf5_myidtype_dataset(buf,"/Subhalo","SubhaloIDMostbound"  ,nsubhalos, 1, &subhaloIDMostbound[subcount]);
          read_hdf5_float_dataset(buf,"/Subhalo","SubhaloMassType"	  ,nsubhalos, 6, &subhaloMassType[6 * subcount]);
          read_hdf5_float_dataset(buf,"/Subhalo","SubhaloMassInRadType"   ,nsubhalos, 6, &subhaloMassInRadType[6 * subcount]);
          read_hdf5_float_dataset(buf,"/Subhalo","SubhaloHalfmassRadType" ,nsubhalos, 6, &subhaloHalfmassRadType[6 * subcount]);
	  read_hdf5_int_dataset(buf,"/Subhalo","SubhaloGrNr"              ,nsubhalos, 1, &subhaloGrNr[subcount]);

#ifdef SUBFIND_EXTRA
          read_hdf5_float_dataset(buf,"/Subhalo","SubhaloGasMetallicity"        ,nsubhalos, 1, &subhaloGasMetallicity[subcount]);
          read_hdf5_float_dataset(buf,"/Subhalo","SubhaloGasMetalFractions"     ,nsubhalos, 9, &subhaloGasMetalFractions[subcount]);
          read_hdf5_float_dataset(buf,"/Subhalo","SubhaloGasMetallicitySfr"     ,nsubhalos, 1, &subhaloGasMetallicitySfr[subcount]);
          read_hdf5_float_dataset(buf,"/Subhalo","SubhaloGasMetalFractionsSfr"  ,nsubhalos, 9, &subhaloGasMetalFractionsSfr[subcount]);
          read_hdf5_float_dataset(buf,"/Subhalo","SubhaloStarMetallicity"       ,nsubhalos, 1, &subhaloStarMetallicity[subcount]);
          read_hdf5_float_dataset(buf,"/Subhalo","SubhaloStarMetalFractions"    ,nsubhalos, 9, &subhaloStarMetalFractions[subcount]);
          read_hdf5_float_dataset(buf,"/Subhalo","SubhaloSFR"  			,nsubhalos, 1, &subhaloSFR[subcount]);
          read_hdf5_float_dataset(buf,"/Subhalo","SubhaloBHMass"                ,nsubhalos, 1, &subhaloBHMass[subcount]); 
          read_hdf5_float_dataset(buf,"/Subhalo","SubhaloBHMdot"                ,nsubhalos, 1, &subhaloBHMdot[subcount]);
          read_hdf5_float_dataset(buf,"/Subhalo","SubhaloSFRinRad"              ,nsubhalos, 1, &subhaloSFRinRad[subcount]);
	  read_hdf5_float_dataset(buf,"/Subhalo","SubhaloStellarPhotometrics"   ,nsubhalos, 8, &subhaloStellarPhotometrics[8 * subcount]);
#endif

          for(i=0;i<nsubhalos;i++)
	    subhaloSnapFileNr[subcount+i] = filenr;

          for(subgr = 0; subgr < nsubhalos; subgr++)
  	    filenrOfHalo[subcount + subgr] = filenr;

          for(subgr = 0; subgr < nsubhalos; subgr++)
	    subhaloindex[subcount + subgr] = subcount + subgr;

	}


      subcount += nsubhalos;
      groupcount += ngroups;
    }



    for(type=0; type < 6; type++)
      {
        if(Cats[num].TotNgroups > 0)
	  {
	    groupOffsetType[type] = 0;
	    if(groupNsubs[0] > 0)
	      subhaloOffsetType[type] = 0;
	  }
      }

    long long int count=0;
    subnr = 0;
    for(subgr = 0; subgr < Cats[num].TotNgroups; subgr++)
      {
        for(type=0; type < 6; type++)
	{
	  if(subgr > 0)
	    {
	      groupOffsetType[subgr*6 + type] = groupOffsetType[(subgr-1)*6 + type] + groupLenType[(subgr-1)*6 + type];
	      if(groupNsubs[subgr]>0)
	        subhaloOffsetType[count*6 + type] = groupOffsetType[subgr*6 + type];
	    }

  	  for(subnr=1; subnr < groupNsubs[subgr]; subnr++)
  	    subhaloOffsetType[(count+subnr)*6 + type] = subhaloOffsetType[(count+subnr-1)*6 + type] + subhaloLenType[(count+subnr-1)*6 + type];
	  
	}
        count += groupNsubs[subgr];
      }

  printf("subfind file data is loaded...\n");
  printf("load data from descendant files...\n");


  if(num < LastSnapShotNr)
    {
      sprintf(buf, "%s/treedata/sub_desc_sf%d_%03d", TreeOutputDir, SnapSkipFac, num);
      if(!(fd = fopen(buf, "r")))
	{
	  printf("can't open file `%s'\n", buf);
	  exit(1);
	}

      printf(" Cats[num].TotNsubhalos = %d\n", Cats[num].TotNsubhalos);
      my_fread(&ncount, sizeof(int), 1, fd);
      printf("ncount = %d\n",ncount);
      my_fread(descendant_haloindex, sizeof(int), Cats[num].TotNsubhalos, fd);
      printf("descendant_haloindex[0] = %d\n",descendant_haloindex[0]);
      my_fread(descendant_snapnum, sizeof(int), Cats[num].TotNsubhalos, fd);
      printf("descendant_snapnum[0] = %d\n",descendant_snapnum[0]);

      fclose(fd);
    }
  printf("descendant file data is loaded...\n");

  nh = FirstHaloInSnap[num];
  sc = 0;


  for(gr = 0; gr < Cats[num].TotNgroups; gr++)
    {
      for(subgr = 0, gr_nh = nh; subgr < groupNsubs[gr]; subgr++, sc++, nh++)
	{
	  /* assign properties */
	  Halo[nh].FirstHaloInFOFgroup = gr_nh;

	  if(subgr == groupNsubs[gr] - 1)
	    Halo[nh].NextHaloInFOFgroup = -1;
	  else
	    Halo[nh].NextHaloInFOFgroup = nh + 1;

	  if(num < LastSnapShotNr)
	    {
	      if(descendant_haloindex[sc] >= 0)
		Halo[nh].Descendant = FirstHaloInSnap[descendant_snapnum[sc]] + descendant_haloindex[sc];
	      else
		Halo[nh].Descendant = -1;
	    }
	  else
	    Halo[nh].Descendant = -1;

	  Halo[nh].FirstProgenitor = -1;
	  Halo[nh].NextProgenitor = -1;


	  if(subgr == 0)
	    {
	      Halo[nh].Group_M_Mean200   = group_M_Mean200[gr];
	      Halo[nh].Group_M_Crit200   = group_M_Crit200[gr];
	      Halo[nh].Group_M_TopHat200 = group_M_TopHat200[gr];
	    }
	  else
	    {
	      Halo[nh].Group_M_Mean200   = 0;
	      Halo[nh].Group_M_Crit200   = 0;
	      Halo[nh].Group_M_TopHat200 = 0;
	    }

	  for(i = 0; i < 3; i++)
	    {
	      Halo[nh].SubhaloPos[i]  = subhaloPos[3 * sc + i];
	      Halo[nh].SubhaloVel[i]  = subhaloVel[3 * sc + i];
	      Halo[nh].SubhaloSpin[i] = subhaloSpin[3 * sc + i];
	    }

	  Halo[nh].SubhaloLen          = subhaloLen[sc];
	  Halo[nh].SubhaloVelDisp      = subhaloVelDisp[sc];
	  Halo[nh].SubhaloVmax         = subhaloVmax[sc];
	  Halo[nh].SubhaloIDMostBound  = subhaloIDMostbound[sc];
	  Halo[nh].SnapNum             = num;
	  Halo[nh].FileNr              = filenrOfHalo[sc];
	  Halo[nh].SubhaloIndex        = subhaloindex[sc];
	  Halo[nh].SubhalfMassRadius   = subhaloHalfmassRad[sc];

	  HaloAux[nh].UsedFlag = 0;

	  for(i = 0; i < 6; i++)
	    Halo[nh].SubhaloMassType[i] = subhaloMassType[6 * sc + i];

#ifdef SUBFIND_EXTRA
	  for(i = 0 ; i < 9; i++)
	    { 
	      Halo[nh].SubhaloGasMetalFractions[i]    = subhaloGasMetalFractions[9 * sc + i];
	      Halo[nh].SubhaloGasMetalFractionsSfr[i] = subhaloGasMetalFractionsSfr[9 * sc + i];
	      Halo[nh].SubhaloStarMetalFractions[i]   = subhaloStarMetalFractions[9 * sc + i];
	    }

	  Halo[nh].SubhaloGasMetallicity	= subhaloGasMetallicity[sc];
	  Halo[nh].SubhaloGasMetallicitySfr	= subhaloGasMetallicitySfr[sc];
          Halo[nh].SubhaloStarMetallicity       = subhaloStarMetallicity[sc];
	  Halo[nh].SubhaloSFR			= subhaloSFR[sc];
	  Halo[nh].SubhaloBHMass        	= subhaloBHMass[sc];
	  Halo[nh].SubhaloBHMdot        	= subhaloBHMdot[sc];
	  Halo[nh].SubhaloSFRinRad      	= subhaloSFRinRad[sc];
        
          for(i = 0; i < 6; i++)
            {
              Halo[nh].SubhaloMassInRadType[i] = subhaloMassInRadType[6 * sc + i];
              Halo[nh].SubhaloHalfmassRadType[i] = subhaloHalfmassRadType[6 * sc + i];
            }
  
          for (i = 0; i < 8; i++)
            {
              Halo[nh].SubhaloStellarPhotometrics[i] = subhaloStellarPhotometrics[8 * sc + i];
            }
#endif
	  Halo[nh].SubhaloGrNr			= subhaloGrNr[sc];
	  Halo[nh].SubhaloSnapFileNr		= subhaloSnapFileNr[sc];

	  for(i=0; i<6; i++)
	    {
	      Halo[nh].SubhaloLenType[i]    = subhaloLenType[    sc*6 + i];
	      Halo[nh].SubhaloOffsetType[i] = subhaloOffsetType[ sc*6 + i];

	      Halo[nh].GroupLenType[i]    = groupLenType[subhaloGrNr[sc]*6+i];
	      Halo[nh].GroupOffsetType[i] = groupOffsetType[subhaloGrNr[sc]*6+i];
	    }
	}
    }

  for(gr = 0; gr < nh; gr++)
    {
      if(Halo[gr].NextHaloInFOFgroup == gr)
	{
	  printf("bummer! %d\n", gr);
	}
    }

  myfree(descendant_haloindex);
  myfree(descendant_snapnum);
  myfree(filenrOfHalo);
  myfree(subhaloindex);

  myfree(subhaloSnapFileNr);
  myfree(groupOffsetType);
  myfree(subhaloOffsetType);

  myfree(groupNsubs);
  myfree(group_M_Mean200);
  myfree(group_M_Crit200);
  myfree(group_M_TopHat200);
  myfree(groupLenType);

  myfree(subhaloLen);
  myfree(subhaloLenType);
  myfree(subhaloPos);          
  myfree(subhaloVel);         
  myfree(subhaloSpin);        
  myfree(subhaloVelDisp);
  myfree(subhaloVmax);        
  myfree(subhaloHalfmassRad);
  myfree(subhaloIDMostbound);  
  myfree(subhaloMassType); 
  myfree(subhaloMassInRadType);
  myfree(subhaloHalfmassRadType); 
  myfree(subhaloGrNr);  

#ifdef SUBFIND_EXTRA
  myfree(subhaloGasMetalFractions);
  myfree(subhaloGasMetalFractionsSfr);
  myfree(subhaloGasMetallicity);
  myfree(subhaloGasMetallicitySfr);
  myfree(subhaloSFR);
  myfree(subhaloStarMetalFractions);
  myfree(subhaloStarMetallicity);
  myfree(subhaloBHMass);
  myfree(subhaloBHMdot);
  myfree(subhaloSFRinRad);
  myfree(subhaloStellarPhotometrics);
#endif
}


void set_progenitor_pointers(void)
{
  int i, first, desc;

  printf("setting progenitor pointers...\n");
  fflush(stdout);

  for(i = 0; i < TotHalos; i++)
    {
      if((desc = Halo[i].Descendant) >= 0)
	{
	  if((first = Halo[desc].FirstProgenitor) >= 0)
	    {

	      if(Halo[i].SubhaloLen >= Halo[first].SubhaloLen)
		{
		  Halo[i].NextProgenitor = first;
		  Halo[desc].FirstProgenitor = i;
		}
	      else
		{
		  Halo[i].NextProgenitor = Halo[first].NextProgenitor;
		  Halo[first].NextProgenitor = i;
		}
	    }
	  else
	    {
	      Halo[desc].FirstProgenitor = i;
	    }
	}
    }
}

