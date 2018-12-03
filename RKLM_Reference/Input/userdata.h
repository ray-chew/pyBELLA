/*******************************************************************************
 File:   userdata.h
 Author: Rupert
 Date:   Thu May 16 16:17:40 WET 2003 
 *******************************************************************************/
#ifndef USERDATA_H
#define USERDATA_H

#include "enumerator.h"
#include "enum_bdry.h"
#include "variable.h"
#include "explicit.h"
#include "limiter.h"
#include "mpv.h"
#include <stdio.h>


/*------------------------------------------------------------------------------
 User_Data
 ------------------------------------------------------------------------------*/
typedef struct {
	
	/* Administration */
	enum Switch write_stdout;
	int write_stdout_period;
	enum Switch write_file;
	int write_file_period;
	enum FileFormat file_format;
	char file_name[FILENAME_MAX];
    char test_dir_name[FILENAME_MAX];
	int ncache;
    int n_time_series;
    
    /* reference quantities */
    double h_ref;
    double t_ref;
    double T_ref;
    double p_ref;
    double rho_ref;
    double u_ref;
    double Nsq_ref;
    double g_ref;
    
	/* Numerical parameters */
	
	/* Time discretization */
	double CFL;
	int stepmax;
	double tout[20];
    double dtfixed0;
	double dtfixed;
    enum TimeIntegrator time_integrator;
    TimeIntegratorParams tips;

    /* Space discretization */
	enum RecoveryOrder recovery_order;
	enum LimiterType limiter_type_scalars;
	enum LimiterType limiter_type_velocity;
	double kp; 
	double kz; 
	double km; 
	double kY; 
	double kZ; 
    
	/* Grid */
	int inx;
	int iny;
	int inz;
	double xmin;
	double xmax;
	double ymin;
	double ymax;
	double zmin;
	double zmax;
	
	/* initial conditions -- highly case-dependent! */
	double wind_speed;
    double wind_shear;
	double hill_height;
	double hill_length_scale;
	
	/* Boundary condition */
	enum BdryType bdrytype_min[3];
	enum BdryType bdrytype_max[3];
	enum Boolean absorber;
	
	/* Thermodynamics and chemistry */
    int nspec;
	double gamm;
	double Rg_over_Rv; 
	double Q; 
	double Ea;
	double B; 
	double Tswitch;
	double Tstar;  
	
	/* Low Mach */
	int acoustic_timestep;
    
    int is_nonhydrostatic;
    double nonhydrostasy;
	int is_compressible;
	double compressibility;
    
    double acoustic_order;
	
    double Msq;
	
	/* Geo-stuff */
	double gravity_strength[3];
	int i_gravity[3];
	int gravity_direction;

    double coriolis_strength[3];
	int i_coriolis[3];
	int coriolis_direction;
    
	/* linear solver stuff */
	double flux_correction_precision;
	double flux_correction_local_precision;
	double second_projection_precision;
	double second_projection_local_precision;
	int flux_correction_max_iterations;
	int second_projection_max_iterations;
    enum Boolean initial_projection;
    enum Boolean initial_impl_Euler;
    enum Boolean column_preconditioner;
    enum Boolean synchronize_nodal_pressure;
    double synchronize_weight;
    /* auxiliary:  effective machine accuracy */
	double eps_Machine;
} User_Data;


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void User_Data_init(User_Data* userdata);


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void Sol_initial(
				 ConsVars* Sol, 
				 const ElemSpaceDiscr* elem,
                 const NodeSpaceDiscr* node);

double stratification(
					  double y);


#endif /* USERDATA_H */


