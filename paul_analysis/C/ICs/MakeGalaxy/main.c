#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"

#include "prototypes.h"
#include "globvars.h"



int main(int argc, char *argv[])
{
  char parameter_filename[100];
  char buf[1000];
  FILE *fd;
  int rep;

  if(argc != 2)
    {
      fprintf(stderr, "\n\nwrong argument(s).  Specify a parameterfile.\n\n");
      exit(0);
    }
  strcpy(parameter_filename, argv[1]);

  read_parameterfile(parameter_filename);

  init_units();			/* set system of units */

  structure();			/* determine structure of halo, disk, and bulge */

  init();			/* allocate arrays */

  init_clouds();


  if(M_GAS > 0)
    {
      init_central_densities();
      
      integrate_and_adjust(0);
      
      // XXX How do we know that 5 iterations is enough to converge here?
      for(rep = 0; rep < 5; rep++)
	{
	  printf("Iteration %d:\n",rep);
	  set_dummy_particles();
	  
	  compute_vertical_force_field();
	  
	  integrate_and_adjust(1);

	  /*
	  // for debugging
	  char gsdfile[50]="";
	  sprintf(gsdfile,"gsddump-%d",rep);
	  dump_gas_density(gsdfile);
	  */
	}

      char gsdfile[50]="";
      
      if(strstr(OutputFile,".hdf5"))
	strncpy(gsdfile, OutputFile, strlen(OutputFile)-5);
      if(strstr(OutputFile,".dat"))
	strncpy(gsdfile, OutputFile, strlen(OutputFile)-4);
      strcat(gsdfile, ".gsd");
      dump_gas_density(gsdfile);
    }
  else
    {
      if(M_DISK > 0)
	{
	  set_dummy_particles();
	  compute_vertical_force_field();
	}
    }


  compute_radial_force_field();

  compute_phi_field();

  compute_DF_lookuptable();

  compute_vstream_gas();
  compute_velocity_dispersions_disk();
  compute_velocity_dispersions_halo();
  compute_velocity_dispersions_bulge();

  set_gas_positions();
  set_halo_positions();
  set_disk_positions();
  set_bulge_positions();

  set_halo_velocities();
  set_disk_velocities();
  set_gas_velocities();
  set_bulge_velocities();


  if(M_DISK > 0)
	dump_veldisp_field();

  mkdir(OutputDir, 02755);

#ifdef HDF5
  /* given as *.dat */
  if(strstr(OutputFile,"dat"))
    {
	strcpy(buf, "");
	strncpy(buf, OutputFile, strlen(OutputFile)-3);
	strcat(buf, "hdf5");
	strcpy(OutputFile, buf);
    }

  /* add .hdf5 if it doesn't have one */
  if(!(strstr(OutputFile,"hdf5")))
	strcat(OutputFile,".hdf5");
#endif

  sprintf(buf, "%s/%s", OutputDir, OutputFile);

  save_particles(buf);

/*sprintf(buf, "%s/curve.txt", OutputDir);*/
  strcpy(buf, "");
  if(strstr(OutputFile,".hdf5"))
	snprintf(buf, strlen(OutputFile)-4, "%s", OutputFile);
        /* strncpy(buf, OutputFile, strlen(OutputFile)-5) */;
  if(strstr(OutputFile,".dat"))
        strncpy(buf, OutputFile, strlen(OutputFile)-4);
  strcat(buf, ".vc");
  printf("writing header, vc, toomre Q to file: %s\n",buf);
  if((fd = fopen(buf, "w")))
    {
      printf("writing circular velocity curve + Toomre's Q\n");
      write_header(fd);
      plot_circular_speeds(fd);
      plot_toomre_stability(fd);
      fclose(fd);
    }
  else
    {
      fprintf(stderr, "Can't open file '%s'.\n", buf);
      exit(0);
    }
  printf("done.\n");


  strcpy(buf, "");
  sprintf(buf, "%s/%s", OutputDir, OutputFile);
  strcat(buf, ".parameters");
  if((fd=fopen(buf,"w")))
    {
      write_header(fd);
      fclose(fd);
    }
  else
    {
      fprintf(stderr,"Can't open file '%s'.\n",buf);
      exit(0);
    }
  printf("done.\n");



  printf("Disk scale length: %g\n", H);
  printf("R200: %g\n", R200);

  return 0;
}



void write_header(FILE * fd)
{
  fprintf(fd, "C     \t%g\n", CC);
  fprintf(fd, "V200  \t%g\n", V200);
  fprintf(fd, "R200  \t%g\n", R200);
  fprintf(fd, "M200  \t%g\n", M200);
  fprintf(fd, "RH    \t%g\n", RH);
  fprintf(fd, "LAMBDA\t%g\n", LAMBDA);
  fprintf(fd, "MD    \t%g\n", MD);
  fprintf(fd, "JD    \t%g\n", JD);
  fprintf(fd, "Rd    \t%g\n", H);
  fprintf(fd, "z0    \t%g\n", DiskHeight);
  fprintf(fd, "sig_r \t%g\n", RadialDispersionFactor);
  fprintf(fd, "MB    \t%g\n", MB);
  fprintf(fd, "BSize \t%g\n", BulgeSize);
  fprintf(fd, "a     \t%g\n", A);
  fprintf(fd, "AniR  \t%g\n", AnisotropyRadius);
  fprintf(fd, "f     \t%g\n", GasFraction);
  fprintf(fd, "GasD  \t%d\n", GasDistribution);
  fprintf(fd, "GasAl \t%g\n", GasExpAlpha);
  fprintf(fd, "PLCut \t%g\n", PowerLawCutoff);
  fprintf(fd, "MBH   \t%g\n", MBH);
  fprintf(fd, "N_HALO\t%d\n", N_HALO);
  fprintf(fd, "N_GAS \t%d\n", N_GAS);
  fprintf(fd, "N_DISK\t%d\n", N_DISK);
  fprintf(fd, "N_BULG\t%d\n", N_BULGE);
  fprintf(fd, "h     \t%g\n", H0*10.0);
  fprintf(fd, "Omeg_m\t%g\n", Omega_m0);
  fprintf(fd, "Omeg_L\t%g\n", Omega_L0);
  fprintf(fd, "z     \t%g\n", REDSHIFT);
  fprintf(fd, "t_SF  \t%g\n", MaxSfrTimescale);
  fprintf(fd, "e_SN  \t%g\n", FactorSN);
  fprintf(fd, "A_0   \t%g\n", FactorEVP);
  fprintf(fd, "T_SN  \t%g\n", TempSupernova);
  fprintf(fd, "T_clds\t%g\n", TempClouds);
  fprintf(fd, "q_0   \t%g\n", FactorForSofterEQS);
  fprintf(fd, "epp   \t%g\n", GravSoftening);
  fprintf(fd, "done\n");
}



