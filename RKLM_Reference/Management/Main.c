    /*******************************************************************************
 File:   Main.c
 Author: Rupert 
 Date:   ?
 *******************************************************************************/
#include "Common.h"
#include <stdio.h>
#include <float.h>
#ifdef MACPROFILE
#include <profiler.h> 
#endif
#include "error.h"
#include "math_own.h"
#include "io.h"
#include "memory.h"
#include "explicit.h"
#include "Gasdynamics.h"
#include "thermodynamic.h"
#include "data.h"
#include "main.h"
#include "enumerator.h"
#include "Eos.h"
#include "recovery.h"
#include "boundary.h"
/* #include "space_discretization.h" */
#include "flux_correction.h"
#include "second_projection_bilinear_p.h"
#include "set_ghostcells_p.h"
#include "mpv.h"
#include "boundary.h"
#include "numerical_flux.h"
#include "enumerator.h"

static void (*rotate[])(ConsVars* Sol, double* rhs, double *Sbg, double *buoyS, const enum Direction dir) = {NULL, rotate2D, rotate3D};

int i_OpSplit;

/* ============================================================================= */

int main( void )
{
	/* low Mach */
	extern MPV* mpv;
	
	/* User data */
	extern User_Data ud;
	extern BDRY* bdry;
	
	/* Time discretization */
	extern int step;  
	extern double t;
	extern double dt;
	
	/* Grid and space discretization */
	extern ElemSpaceDiscr* elem;
	extern NodeSpaceDiscr* node;
	
	/* Arrays */
	extern ConsVars* Sol; 
	extern ConsVars* Sol0; 
	extern ConsVars* flux[3];
    extern VectorField* adv_flux;
    extern VectorField* adv_flux0;
    extern double* buoyS;
    extern VectorField* buoy; 
    extern VectorField* buoy0; 
    
    extern double* W0;
    extern double* W1;
    extern double* Sbg;
	
	Speeds a_u_max;
    TimeStepInfo dt_info;
	double lambda, cfl, cfl_ac, cfl_adv;   
	const double* tout = ud.tout;
    const int sequence = 1;  /*  xyyx sequence (good)  1;   yxxy sequence (bad) 0; */
	int Split = 0;
	int output_switch = 0;
	int time_step_switch = 0;
    int stage;
    int which_projection = ud.which_projection_first;
         
    /* 
     const double stepfrac[] = {1.0, 2.0, 0.0, 1.0}; 
     */
    const double stepfrac[] = {1.0, 1.0, 1.0, 1.0};
    
    const enum GravityTimeIntegrator GRAVITY_SWITCH_I  = (stepfrac[1] == 2.0 ? IMPLICIT_TRAPEZOIDAL : EULER_BACKWARD);
    const enum GravityTimeIntegrator GRAVITY_SWITCH_II = (stepfrac[2] == 0.0 ? IMPLICIT_TRAPEZOIDAL : EULER_FORWARD);
    
    int substep;
        
    enum Boolean reset_init_data = CORRECT; /* CORRECT; WRONG; */
    
    FILE *tsfile = NULL; 
    char tsfilename[100];
    sprintf(tsfilename, "timeseries.txt");
    tsfile = fopen(tsfilename, "w+");

    a_u_max.u              = 999999.999;
	a_u_max.u_plus_c       = 888888.888;
	cfl = cfl_ac = cfl_adv = 777777.777;
	
	/* data allocation and initialization */
	Data_init();

    set_wall_massflux(bdry, Sol, elem);
    Set_Explicit_Boundary_Data(Sol, elem, mpv, 1);
    ud.compressibility = compressibility(0);
    
    VectorField* adv_flux_diff = adv_flux0;
    
	if(ud.write_file == ON) 
        putout(Sol, t, 0.0, 0, 0, ud.file_name, "Sol", 1);
    
    if (ud.write_history) {
        WriteTimeHistories(Sol, elem, 0.0, 0, 1);
    }
    
    ConsVars_set(Sol0, Sol, elem->nc);
    	
    double* Zaux_p = W0; 
    double* Zaux_b = W1; 
    
    /* generate divergence-controlled initial data  */
    dt = 1.0;
    mpv->dt = dt;

    ConsVars_set(Sol0, Sol, elem->nc);
    which_projection = 1;

	/* Main loop over the sequence of time values of tout */
	while(t < *tout && step < ud.stepmax) {
		
		/* start major time stepping loop */
		while( t < *tout && step < ud.stepmax ) { 
			            
            /* initialize fluxes in preparation of explicit predictor */
            ConsVars_setzero(flux[0], elem->nfx);
            if(elem->ndim > 1) ConsVars_setzero(flux[1], elem->nfy);
            if(elem->ndim > 2) ConsVars_setzero(flux[2], elem->nfz);            

			/* Timestep computation */
			if(time_step_switch) {
				time_step_switch = 0;
			}
			else {
                dt_info          = dynamic_timestep(Sol, t, *tout, elem, step);
                dt               = dt_info.time_step;
                mpv->dt = dt;
                cfl              = dt_info.cfl;
                cfl_ac           = dt_info.cfl_ac;
                cfl_adv          = dt_info.cfl_adv;
                time_step_switch = dt_info.time_step_switch;
			}
			            
			set_wall_massflux(bdry, Sol0, elem);
                       
            if (ud.time_integrator == OP_SPLIT || ud.time_integrator == OP_SPLIT_MD_UPDATE) {

                if (reset_init_data) {
                    adjust_pi(Sol, mpv, Sol0, elem, 1.0);
                    for (int n=0; n<elem->nc; n++) {
                        Zaux_p[n] = Sol->rhoZ[PRES][n];
                        Zaux_b[n] = Sol->rhoZ[SOLD][n];
                    }
                    ConsVars_set(Sol, Sol0, elem->nc);                
                    for (int n=0; n<elem->nc; n++) {
                        Sol->rhoZ[PRES][n] = Zaux_p[n];
                        Sol->rhoZ[SOLD][n] = Zaux_b[n];
                    }
                    if (step == 0) {
                        reset_init_data = WRONG;
                        step = 0;
                        t = 0.0;
                    }
                } else {
                    adjust_pi(Sol, mpv, Sol0, elem, 1.0);
                    ConsVars_set(Sol0, Sol, elem->nc);            
                }
                
#if OUTPUT_SUBSTEPS
                putout(Sol, t, *tout , step, 0, ud.file_name, "Sol", 1);
#endif

                printf("\nCoriolis 1 ---------------------------------- \n");                
                Explicit_Coriolis(Sol, elem, 0.5*dt);                
#if OUTPUT_SUBSTEPS
                putout(Sol, t, *tout , step, 0, ud.file_name, "Sol", 1);
#endif

                printf("\nnonlinear fluxes ---------------------------- \n");
                
                substep = 0;
                
                /* valid for 2D only */
                assert(elem->ndim < 3);
                if (sequence == 0) {
                    (*rotate[elem->ndim - 1])(Sol, mpv->dp2_cells, Sbg, buoyS, FORWARD);
                }
                
                /* FORWARD gasdynamics */
                stage = 0;
                for(i_OpSplit = 0; i_OpSplit < elem->ndim; i_OpSplit++) {
                    lambda = stepfrac[substep]*ud.tips.dt_frac*dt/elem->dx;
                    Split = sequence * i_OpSplit + (1 - sequence) * ((elem->ndim - 1) - i_OpSplit);
                    Explicit_step_and_flux(Sol, flux[Split], buoyS, buoy, mpv->dp2_cells, mpv->HydroState, lambda, elem->nc, Split, stage, FLUX_INTERNAL, GRAVITY_SWITCH_I, WITH_MUSCL, WITH_GRAVITY);
                    substep++;
                    
#if OUTPUT_SUBSTEPS_PREDICTOR
                    /* TODO: remove necessity of calling b.c. routine for all directions after each splitstep */
                    if (Split == 1)  {
                        (*rotate[elem->ndim - 1])(Sol, mpv->dp2_cells, Sbg, buoyS, BACKWARD);
                    }
                    Set_Explicit_Boundary_Data(Sol, elem, mpv, 0);
                    putout(Sol, t, *tout , step, 0, ud.file_name, "Sol", OUTPUT_SPLITSTEPS);
                    if (Split == 1)  {
                        (*rotate[elem->ndim - 1])(Sol, mpv->dp2_cells, Sbg, buoyS, FORWARD);
                    }
#endif
                    if(i_OpSplit < elem->ndim - 1) (*rotate[elem->ndim - 1])(Sol, mpv->dp2_cells, Sbg, buoyS, FORWARD);
                }
                                
                /* BACKWARD gasdynamics */
                stage = 1;
                for(i_OpSplit = 0; i_OpSplit < elem->ndim; i_OpSplit++) {
                    lambda = stepfrac[substep]*ud.tips.dt_frac*dt/elem->dx;
                    Split = (1 - sequence) * i_OpSplit + sequence * ((elem->ndim - 1) - i_OpSplit);
                    Explicit_step_and_flux(Sol, flux[Split], buoyS, buoy, mpv->dp2_cells, mpv->HydroState, lambda, elem->nc, Split, stage, FLUX_INTERNAL, GRAVITY_SWITCH_II, WITH_MUSCL, WITH_GRAVITY);
                    substep++;
                    
#if OUTPUT_SUBSTEPS_PREDICTOR
                    if (Split == 1)  {
                        (*rotate[elem->ndim - 1])(Sol, mpv->dp2_cells, Sbg, buoyS, BACKWARD);
                    }
                    Set_Explicit_Boundary_Data(Sol, elem, mpv, 0);
                    putout(Sol, t, *tout , step, 0, ud.file_name, "Sol", OUTPUT_SPLITSTEPS);
                    if (Split == 1)  {
                        (*rotate[elem->ndim - 1])(Sol, mpv->dp2_cells, Sbg, buoyS, FORWARD);
                    }
#endif
                    if(i_OpSplit < elem->ndim - 1) (*rotate[elem->ndim - 1])(Sol, mpv->dp2_cells, Sbg, buoyS, BACKWARD);
                }
                
                /* valid for 2D only */
                assert(elem->ndim < 3);
                if (sequence == 0) {
                    (*rotate[elem->ndim - 1])(Sol, mpv->dp2_cells, Sbg, buoyS, BACKWARD);
                }
                
                if (ud.time_integrator == OP_SPLIT_MD_UPDATE) {
#if OUTPUT_SUBSTEPS_PREDICTOR
                    putout(Sol, t, *tout , step, 0, ud.file_name, "Sol", 1);
#endif
                    fullD_explicit_updates(Sol, Sol0, flux, buoyS, buoy, elem, dt, stage);
#if OUTPUT_SUBSTEPS_PREDICTOR
                    putout(Sol, t, *tout , step, 0, ud.file_name, "Sol", 1);
#endif
                }
                
                /* Bracketing Coriolis in Strang Split fashion around full advection-correction,
                 that is after second projection, seems slightly better (see below) */
                printf("\nCoriolis 2 ---------------------------------- \n");                
                Explicit_Coriolis(Sol, elem, 0.5*dt); 
                
                Set_Explicit_Boundary_Data(Sol, elem, mpv, 0);
                
#if OUTPUT_SUBSTEPS
                putout(Sol, t, *tout , step, 0, ud.file_name, "Sol", 1);
#endif
                
                if(PROJECTION1) {
                    
                    flux_correction(flux, adv_flux_diff, buoy, elem, Sol, Sol0, t, dt, ud.implicitness, step);
                    update(Sol, (const ConsVars**)flux, buoyS, buoy, elem, dt);
                    Set_Explicit_Boundary_Data(Sol, elem, mpv, 1);
#if OUTPUT_SUBSTEPS
                    putout(Sol, t, *tout , step, 0, ud.file_name, "Sol", 1);
#endif
                }
                
#if OUTPUT_SUBSTEPS
                putout(Sol, t, *tout , step, 0, ud.file_name, "Sol", 1);
#endif
                
                if(PROJECTION2 == 1) {
                    second_projection(Sol, mpv, (const ConsVars*)Sol0, elem, node, 1.0, t, dt);
                    Set_Explicit_Boundary_Data(Sol, elem, mpv, 0);
#if OUTPUT_SUBSTEPS
                    putout(Sol, t, *tout , step, 0, ud.file_name, "Sol", 1);
#endif
                }
                
                if (ud.absorber) 
                {
                    Absorber(Sol, t, dt); 
                    Set_Explicit_Boundary_Data(Sol, elem, mpv, 0);
                }					                
            }
            else if (ud.time_integrator == SI_MIDPT) {
                
                /* ======================================================================= */
                /* Semi-implicit midpoint taken more seriously                             */
                /* ======================================================================= */
                
                reset_Y_perturbation(Sol, (const MPV*)mpv, elem);
#if OUTPUT_SUBSTEPS
                putout(Sol, t, *tout , step, 0, ud.file_name, "Sol", 1);
#endif

                /* ----------------------------------------------------------------------- */
                /* half-time prediction of advective flux                                  */
                /* ----------------------------------------------------------------------- */
                printf("\nhalf-time prediction of advective flux ----------- \n");

                substep = 0;
                
                /* explicit advection half time step preparing advection flux calculation */
                stage = 0;
                for(Split = 0; Split < elem->ndim; Split++) {
                    lambda = 0.5*dt/elem->dx;
                    Explicit_step_and_flux(Sol, flux[Split], buoyS, buoy, mpv->dp2_cells, mpv->HydroState, lambda, elem->nc, Split, stage, FLUX_INTERNAL, EULER_FORWARD, WITH_MUSCL, WITHOUT_GRAVITY);
                    substep++;
                    
#if OUTPUT_SUBSTEPS_PREDICTOR
                    if (Split == 1)  {
                        (*rotate[elem->ndim - 1])(Sol, mpv->dp2_cells, Sbg, buoyS, BACKWARD);
                    }
                    Set_Explicit_Boundary_Data(Sol, elem, mpv, 0);
                    putout(Sol, t, *tout , step, 0, ud.file_name, "Sol", OUTPUT_SPLITSTEPS);
                    if (Split == 1)  {
                        (*rotate[elem->ndim - 1])(Sol, mpv->dp2_cells, Sbg, buoyS, FORWARD);
                    }
#endif
                    (*rotate[elem->ndim - 1])(Sol, mpv->dp2_cells, Sbg, buoyS, FORWARD);
                    
                }
                fullD_explicit_updates(Sol, Sol0, flux, buoyS, buoy, elem, dt, stage);
#if OUTPUT_SUBSTEPS
                putout(Sol, t, *tout , step, 0, ud.file_name, "Sol", 1);
#endif
                
                /* explicit part of Euler backward gravity over half time step */
                VectorField_setzero(buoy0, elem->nc);
                euler_backward_gravity(Sol, (const MPV*)mpv, elem, 0.5*dt);
                Set_Explicit_Boundary_Data(Sol, elem, mpv, 1);
#if OUTPUT_SUBSTEPS
                putout(Sol, t, *tout , step, 0, ud.file_name, "Sol", 1);
#endif

                /* divergence-controlled advective fluxes at the half time level */
                recompute_advective_fluxes(flux, (const ConsVars*)Sol, elem);
                store_advective_fluxes(adv_flux, (const ConsVars**)flux, elem);
                flux_correction(flux, adv_flux_diff, buoy, elem, Sol, Sol0, t, dt, ud.implicitness, step);                
                update_advective_fluxes(flux, (const VectorField*)adv_flux, elem, node, dt);    
#if OUTPUT_SUBSTEPS
                putout(Sol, t, *tout , step, 0, ud.file_name, "Sol", 1);
#endif
                
                /* reset Sol to state at beginning of time step */
                ConsVars_set(Sol, Sol0, elem->nc);
#if OUTPUT_SUBSTEPS
                putout(Sol, t, *tout , step, 0, ud.file_name, "Sol", 1);
#endif
                
                /* explicit EULER half time step for gravity and pressure gradient */ 
                VectorField_setzero(buoy0, elem->nc);
                euler_forward_non_advective(Sol, (const MPV*)mpv, elem, node, 0.5*dt); 
                
#if OUTPUT_SUBSTEPS
                putout(Sol, t, *tout , step, 0, ud.file_name, "Sol", 1);
#endif
                
                /* explicit full time step advection using div-controlled advective fluxes */
                printf("\nnonlinear advection ------------------------------ \n");
                
                substep = 0;
                VectorField_setzero(buoy, elem->nc);

                stage = 0;
                for(Split = 0; Split < elem->ndim; Split++) {
                    lambda = 0.5*dt/elem->dx;
                    Explicit_step_and_flux(Sol, flux[Split], buoyS, buoy, mpv->dp2_cells, mpv->HydroState, lambda, elem->nc, Split, stage, FLUX_EXTERNAL, EULER_BACKWARD, WITH_MUSCL, WITHOUT_GRAVITY);
                    substep++;
                    
#if OUTPUT_SUBSTEPS_PREDICTOR
                    if (Split == 1)  {
                        (*rotate[elem->ndim - 1])(Sol, mpv->dp2_cells, Sbg, buoyS, BACKWARD);
                    }
                    Set_Explicit_Boundary_Data(Sol, elem, mpv, 0);
                    putout(Sol, t, *tout , step, 0, ud.file_name, "Sol", OUTPUT_SPLITSTEPS);
                    if (Split == 1)  {
                        (*rotate[elem->ndim - 1])(Sol, mpv->dp2_cells, Sbg, buoyS, FORWARD);
                    }
#endif
                    (*rotate[elem->ndim - 1])(Sol, mpv->dp2_cells, Sbg, buoyS, FORWARD);
                }
                
                stage = 1;
                for(i_OpSplit = 0; i_OpSplit < elem->ndim; i_OpSplit++) {

                    (*rotate[elem->ndim - 1])(Sol, mpv->dp2_cells, Sbg, buoyS, BACKWARD);

                    lambda = stepfrac[substep]*ud.tips.dt_frac*dt/elem->dx;
                    Split = (1 - sequence) * i_OpSplit + sequence * ((elem->ndim - 1) - i_OpSplit);
                    Explicit_step_and_flux(Sol, flux[Split], buoyS, buoy, mpv->dp2_cells, mpv->HydroState, lambda, elem->nc, Split, stage, FLUX_EXTERNAL, EULER_FORWARD, WITH_MUSCL, WITHOUT_GRAVITY);
                    substep++;
                    
#if OUTPUT_SUBSTEPS_PREDICTOR
                    if (Split == 1)  {
                        (*rotate[elem->ndim - 1])(Sol, mpv->dp2_cells, Sbg, buoyS, BACKWARD);
                    }
                    Set_Explicit_Boundary_Data(Sol, elem, mpv, 0);
                    putout(Sol, t, *tout , step, 0, ud.file_name, "Sol", OUTPUT_SPLITSTEPS);
                    if (Split == 1)  {
                        (*rotate[elem->ndim - 1])(Sol, mpv->dp2_cells, Sbg, buoyS, FORWARD);
                    }
#endif
                }
                Set_Explicit_Boundary_Data(Sol, elem, mpv, 1);

                euler_backward_gravity(Sol, mpv, elem, 0.5*dt);
#if OUTPUT_SUBSTEPS
                putout(Sol, t, *tout , step, 0, ud.file_name, "Sol", 1);
#endif
                
                second_projection(Sol, mpv, (const ConsVars*)Sol0, elem, node, 1.0, t, dt);
                Set_Explicit_Boundary_Data(Sol, elem, mpv, 1);
#if OUTPUT_SUBSTEPS
                putout(Sol, t, *tout , step, 0, ud.file_name, "Sol", 1);
#endif

                /* This placement is not according to Piotr's rules ... */
                update_SI_MIDPT_buoyancy(Sol, (const ConsVars**)flux, mpv, elem, 0.5*dt);
#if OUTPUT_SUBSTEPS
                putout(Sol, t, *tout , step, 0, ud.file_name, "Sol", 1);
#endif

                /* Bracketing Coriolis in Strang Split fashion around full advection-correction,
                 that is after second projection, seems slightly better (see below) */
                printf("\nCoriolis 2 ---------------------------------- \n");                
                Explicit_Coriolis(Sol, elem, 0.5*dt); 
                Set_Explicit_Boundary_Data(Sol, elem, mpv, 0);
                
#if OUTPUT_SUBSTEPS
                putout(Sol, t, *tout , step, 0, ud.file_name, "Sol", 1);
#endif
                                
                if (ud.absorber) 
                {
                    Absorber(Sol, t, dt); 
                    Set_Explicit_Boundary_Data(Sol, elem, mpv, 0);
                }                                    
            }

            
			t += dt;
			step++;
			            
			if((ud.write_stdout == ON && step % ud.write_stdout_period == 0) || output_switch)
			{
				printf("\nstep %d done,  t=%f,  dt=%f,  cfl=%f, cfl_ac=%f, cfl_adv=%f, umax = %f, upcmax = %f\n", step, t, dt, cfl, cfl_ac, cfl_adv, a_u_max.u, a_u_max.u_plus_c );
			}
			
            ud.compressibility = compressibility(t);
            
			if((ud.write_file == ON && (step % ud.write_file_period  == 0)) || output_switch) 
				putout(Sol, t, *tout , step, 0, ud.file_name, "Sol", 1);
			
			if(ud.write_stdout == ON && step % ud.write_stdout_period == 0)
			{
				printf("\n-----------------------------------------------------\n-----------------------------------------------------\n\n");
			}
            
            if (ud.write_history) {
                WriteTimeHistories(Sol, elem, t, step, 0);
            }
		}
        
		if(ud.write_file == ON) putout(Sol, t, *tout , step, 0, ud.file_name, "Sol", 1);
		tout++;
	}
	
	/* data release */
	Data_free();
	Explicit_free(); 
	recovery_free();
	
    fclose(tsfile);
    
    if (ud.write_history) {
        WriteTimeHistories(Sol, elem, 1000000.0, 1000000, -1);
    }
    
    
	return(1);
}

/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: Main.c,v $
 Revision 1.2  1998/03/07 09:56:44  nicola
 Added flux computation and multiple pressure variables.
 
 Revision 1.1  1998/03/01 18:43:32  nicola
 This is the initial revision of 3d. It comes out after two weeks of work on
 Matthias' version of Rupert's F3D code. It is the starting point for imple
 menting a low Mach number extension.
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
