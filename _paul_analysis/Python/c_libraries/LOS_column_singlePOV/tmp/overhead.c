#include"overhead.h"
struct global_data_all_processes All;
double  Kernel[KERNEL_TABLE+2],
        KernelDer[KERNEL_TABLE+2],
        KernelRad[KERNEL_TABLE+2];
long TOTAL_NUMBER_OF_CELLS;
long ALL_CELL_COUNTER;
struct CELL_STRUCT *CELL;
struct particle_bh *PBH;
struct particle_star *PS;
struct particle_bulge *PB;
struct particle_disk *PD;
struct particle_halo *PH;
struct particle_gas *PG;
struct snap_header header;
