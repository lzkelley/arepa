
#ifndef ALLVARS_H
#define ALLVARS_H

#include <stdio.h>


extern int  LastSnapShotNr;
extern int  FirstSnapShotNr;
extern int  SnapSkipFac;

extern char TreeOutputDir[512];
extern char RunOutputDir[512];
extern char SnapshotFileBase[512];

extern int TotHalos;
extern int NtreesPerFile;
extern int NhalosPerFile;

extern double ParticleMass; 

extern int  NumberOfOutputFiles;
extern int  MultipleFilesPerSnap; // -1 if unset, 0 if NumFilesPerSnapshot=1, 1 if true

extern int    *FirstHaloInSnap;
extern struct halo_catalogue
{
      int TotNsubhalos;
      int TotNgroups;
      float redshift;
}
*Cats;

#ifndef LONGIDS
typedef unsigned int MyIDType;
#define H5T_MYIDTYPE H5T_NATIVE_UINT
#else
typedef unsigned long long int MyIDType;
#define H5T_MYIDTYPE H5T_NATIVE_UINT64
#endif


extern struct halo_data
{
  int Descendant;
  int FirstProgenitor;
  int NextProgenitor;
  int FirstHaloInFOFgroup;
  int NextHaloInFOFgroup;
  int SubhaloLen;
  float Group_M_Mean200;
  float Group_M_Crit200;
  float Group_M_TopHat200;
  float SubhaloPos[3];
  float SubhaloVel[3];
  float SubhaloVelDisp;
  float SubhaloVmax;
  float SubhaloSpin[3];
  MyIDType SubhaloIDMostBound;
  int SnapNum;
  int FileNr;
  int SubhaloIndex;
  float SubhalfMassRadius;
  float SubhaloMassType[6];
  float SubhaloMassInRadType[6];  
  float SubhaloHalfmassRadType[6];

  int SubhaloSnapFileNr;
  int SubhaloGrNr;
  int GroupLenType[6];
  int SubhaloLenType[6];

  MyIDType GroupOffsetType[6];
  MyIDType SubhaloOffsetType[6];

#ifdef SUBFIND_EXTRA
  float SubhaloGasMetallicity;
  float SubhaloGasMetalFractions[9];
  float SubhaloGasMetallicitySfr;
  float SubhaloGasMetalFractionsSfr[9];
  float SubhaloStarMetallicity;
  float SubhaloStarMetalFractions[9];
  float SubhaloSFR;
  float SubhaloBHMass;
  float SubhaloBHMdot;
  float SubhaloSFRinRad;
  float SubhaloStellarPhotometrics[8];
#endif

} *Halo, *HaloList;



extern struct halo_aux_data
{
  int UsedFlag;
  int FileNr;
  int TargetIndex;
  int Origin;
  int HaloFlag;

}  *HaloAux;


enum iofields
{
        IO_DESCENDANT,
        IO_FIRSTPROGENITOR,
        IO_NEXTPROGENITOR,
        IO_FIRSTHALOINFOFGROUP,
        IO_NEXTHALOINFOFGROUP,
        IO_LEN,
        IO_M_MEAN200,
        IO_M_CRIT200,
        IO_M_TOPHAT,
        IO_POS,
        IO_VEL,
        IO_VELDISP,
        IO_VMAX,
        IO_SPIN,
        IO_IDMOSTBOUND,
        IO_SnapNum,
        IO_FILENR,
        IO_GROUPNR,
        IO_SUBHALO_NR,
        IO_SUBHALO_SFR,
        IO_SUBHALO_GAS_METALLICITY,
        IO_SUBHALO_GAS_METALLICITY_SFR,
        IO_SUBHALO_STELLAR_METALLICITY,
        IO_SUBHALO_BH_MASS,
        IO_SUBHALO_BH_MDOT,
        IO_SUBHALO_SFR_IN_RAD,
        IO_SUBHALO_STELLAR_PHOTOMETRICS,
        IO_SUBHALO_OFFSET_TYPE,
        IO_SUBHALO_LEN_TYPE,
        IO_SUBHALO_MASS_TYPE,
        IO_SUBHALO_MASSINRAD_TYPE,
        IO_SUBHALO_HALFMASSRAD_TYPE,
        IO_LASTENTRY
};

#endif

#ifdef SKIP_SNAP
#define MAXLEN_OUTPUTLIST 300

extern int OutputListLength;
extern double OutputListTimes[MAXLEN_OUTPUTLIST];
extern int OutputListFlag[MAXLEN_OUTPUTLIST];

extern char OutputList[512];
#endif





