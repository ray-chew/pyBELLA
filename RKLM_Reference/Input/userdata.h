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
    int write_history;
	enum FileFormat file_format;
	char file_name[FILENAME_MAX];
    char test_dir_name[FILENAME_MAX];
	int ncache;
	int nhist;
	
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
    int no_of_steps_to_CFL;
    int no_of_steps_to_dtfixed;
    enum TimeIntegrator time_integrator;
    TimeIntegratorParams tips;

    /* Space discretization */
    enum Boolean p_flux_correction;
    double latw[3];
    double p_extrapol;

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
	double h;
	
	/* initial conditions -- highly case-dependent! */
	double wind_speed;
    double wind_shear;
	double hill_height;
	double hill_length_scale;
    double epsmu;
	
	/* Boundary condition */
	enum BdryType bdrytype_min[3];
	enum BdryType bdrytype_max[3];
	enum Boolean absorber;
	
	/* Thermodynamics and chemistry */
	double gamm;
    double nspec;
	double Rg_over_Rv; 
	double Q; 
	double Ea;
	double B; 
	double Tswitch;
	double Tstar;  
	
	/* Low Mach */
	int acoustic_timestep;
	int is_compressible;
	double compressibility;
	double Msq;
	double M;
	
	/* Geo-stuff */
	double gravity_strength[3];
	int i_gravity[3];
	int gravity_direction;
    int implicit_gravity_theta;
    int implicit_gravity_press;
    int implicit_gravity_theta2;
    int implicit_gravity_press2;

    double coriolis_strength[3];
	int i_coriolis[3];
	int coriolis_direction;
    
	/* linear solver stuff */
    int which_projection_first;
	enum SolverType Solver;
	enum SolverType Solver_Node;
    enum Boolean precondition;
	double flux_correction_precision;
	double flux_correction_local_precision;
	double second_projection_precision;
	double second_projection_local_precision;
	double implicitness;
	int flux_correction_max_iterations;
	int second_projection_max_iterations;
	int flux_correction_max_MG_cycles;
	int flux_correction_output_period;
	int max_projection_iterations;
	
	/* Multigrid */
	int max_no_of_multigrid_levels;
	int no_of_multigrid_levels;
	
	/* Front related stuff */
	double G_layer_scale;
	double spread_overlap_scale;
	double epsiter;
	double epsVF;  
	double epsVF_source;
	double eps_Dspread; 
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


