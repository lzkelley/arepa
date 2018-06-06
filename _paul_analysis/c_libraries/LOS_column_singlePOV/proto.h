#ifndef INLINE_FUNC
#ifdef INLINE
#define INLINE_FUNC inline
#else
#define INLINE_FUNC
#endif
#endif

/* Paul's functions */
void readsnap(char *filename); 
void printhead(void); 
long allocate_gas(void); 
long allocate_halo(void); 
long allocate_disk(void); 
long allocate_bh(void); 
long allocate_bulge(void);
long allocate_star(void);  
long allocate_(void); 
void free_memory(); 


/* From lineofsight.h */
/* Function declarations */

/* tree functions */
void initialize_the_tree(void);
void build_cell(long cell_ID);
void get_tree_localvals(long cell_ID);
void get_tree_localval_derivs(long cell_ID);
long find_cell_from_scratch(float pos[3]);
long find_cell_from_nearby_cell(float cell_ID_0, float pos[3]);
void carry_ray_across_cell(Ray_struct *ray, long cell_id);
void integrate_ray_to_escape(Ray_struct *ray, long origin_cell_id);

void generate_ray_angles(Ray_struct *R, long N_rays);
void setup_ray(Ray_struct *ray, long nrays, long bh_origin);
void local_vars(float x, float y, float z, LOCALVAL *vals);
void line_of_sight(Ray_struct *ray);
void find_neighbors(float r0[3],long Ntot,long *nbr_list,float *r2_list,float *r2_max);
void integrate_ray(Ray_struct *R);
float simp(float *dx, float *y, long n);
void init_clouds(void);
void set_CGS_units(void);
void set_All_struct_terms(void);
void set_units_sfr(void);
void local_vars(float x, float y, float z, LOCALVAL *vals);
void find_neighbors(float r0[3],long Ntot,long *nbr_list,float *r2_list, float *r2_max);
void get_two_phase_breakdown(double u_tot, double rho_tot, double ne_given, LOCALVAL *local);
void setup_lineofsight_overhead(void);
void free_2d_array(void **array);
void free_3d_array(void ***array);
void** calloc_2d_array(size_t nr, size_t nc, size_t size);
void*** calloc_3d_array(size_t nr, size_t nc, size_t nz, size_t size);


/* From cooling.h */
double AbundanceRatios(double u, double rho, double *ne_guess, double *nH0_pointer);
double convert_u_to_temp(double u, double rho, double *ne_guess);
double CoolingRate(double logT, double rho, double *nelec);
double CoolingRateFromU(double u, double rho, double *ne_guess);
double DoCooling(double u_old, double rho, double dt, double *ne_guess);
double GetCoolingTime(double u_old, double rho,  double *ne_guess);
double DoInstabilityCooling(double m_old, double u, double rho, double dt, double fac, double *ne_guess);
void   find_abundances_and_rates(double logT, double rho, double *ne_guess);
void   InitCool(void);
void   InitCoolMemory(void);
void   IonizeParams(void);
void   IonizeParamsFunction(void);
void   IonizeParamsTable(void);
double INLINE_FUNC LogTemp(double u, double ne);
void   MakeCoolingTable(void);
void   *mymalloc(size_t size);
void   ReadIonizeParams(char *fname);
void   SetZeroIonization(void);
void   TestCool(void);


/* From proto.h */
void   advance(void);
void   allocate_commbuffers(void);
void   allocate_memory(void);
void   allocate_memory_2d(void);
void   begrun(void);
void   close_outputfiles(void);
void   check_omega(void);
void   compute_accelerations(long mode);
void   compute_global_quantities_of_system(void);
void   compute_potential(void);
void   construct_timetree(void);
void   delete_node(long i);
void   density(void);
//void   grid_density(struct particle_data GridP[], struct sph_particle_data SphGridP[], long ngridp);
void   determine_interior(void);
double dmax(double,double);
double dmin(double,double);
void   do_box_wrapping(void);
void   DomainDecomposition(void);
double drand48();
void   endrun(long);
void   energy_statistics(void);
void   ensure_neighbours(long mode);
void   every_timestep_stuff(void);
void   INLINE_FUNC ewald_corr(double dx, double dy, double dz, double *fper);
void   ewald_force(double x[3], double force[3]);
void   ewald_init(void);
double ewald_psi(double x[3]);
double INLINE_FUNC ewald_pot_corr(double dx, double dy, double dz);
long    find_ancestor(long i);
double find_next_outputtime(double time);
void   find_next_time(void);
long    find_next_time_walk(long node);
void   free_memory(void);
void   free_memory_2d(void);
void   find_timesteps(long mode);   
void   gravity_tree(void);
void   hydro_force(void);
long    imax(long,long);
long    imin(long,long);
//void   init(long ngridp);
void   init(void);
void   init_2d(long ngridp);
void   inquire_about_grape_system(void);
void   insert_node(long i);
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE *stream);
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE *stream);
void   open_outputfiles(void);
double INLINE_FUNC periodic(double x);
double pot_integrand(double xx);
void   predict(double time);
void   predict_collisionless_only(double time);
void   predict_sph_particles(double time);
void   read_ic(char *fname);
//void   read_ic(char *fname, long ngridp);
void   read_ic_cluster(char *fname);
void   read_ic_cluster_gas(char *fname);
void   read_ic_cluster_wimp(char *fname);
long    read_outputlist(char *fname);
void   read_parameter_file(char *fname);
void   restart(long mod);
void   run(void);
void   savepositions(long num);
void   savepositions_ioformat1(long num);
double second(void);
void   set_softenings(void);
void   set_sph_kernel(void);
void   set_units(void);
//void   setup_smoothinglengths(long desired_ngb, long ngridp);
void   setup_smoothinglengths(long desired_ngb);
void   setup_smoothinglengths_veldisp(long desired_ngb);
void   statistics(void);
double timediff(double t0,double t1);
void   veldisp(void);
