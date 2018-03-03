/*******************************************************************************
 File:   Main.c
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
    
    extern double* Sbg;
    
    enum FluxesFrom advec_flux;
	
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
         
    const double stepfrac[] = {1.0, 1.0, 1.0, 1.0};
        
    int substep;
            
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
                       
            /* ======================================================================= */
            /* Semi-implicit discretization of non-advective terms aka EULAG           */
            /* ======================================================================= */
            
            reset_Y_perturbation(Sol, (const MPV*)mpv, elem);
            ConsVars_set(Sol0, Sol, elem->nc);            
            
#if OUTPUT_SUBSTEPS  /* 1 */
            putout(Sol, t, *tout , step, 0, ud.file_name, "Sol", 1);
#endif
            
            /* ----------------------------------------------------------------------- */
            /* half-time prediction of advective flux                                  */
            /* ----------------------------------------------------------------------- */
            printf("\nhalf-time prediction of advective flux ----------- \n");
            
            substep = 0;
            
#ifdef HALF_STEP_FLUX_EXTERNAL
            recompute_advective_fluxes(flux, (const ConsVars*)Sol, elem);
            advec_flux = FLUX_EXTERNAL;
#else
            advec_flux = FLUX_INTERNAL;
#endif
            
            printf("\nCoriolis 0 ---------------------------------- \n");                
            Explicit_Coriolis(Sol, elem, 0.5*dt); 
            Set_Explicit_Boundary_Data(Sol, elem, mpv, 0);
            
            /* explicit advection half time step preparing advection flux calculation */
            stage = 0;
            for(Split = 0; Split < elem->ndim; Split++) {
                lambda = 0.5*dt/elem->dx;
                Explicit_step_and_flux(Sol, flux[Split], lambda, elem->nc, Split, stage, advec_flux, WITH_MUSCL);
                substep++;
                
#if OUTPUT_SUBSTEPS_PREDICTOR  /* 2, 3 */
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
            fullD_explicit_updates(Sol, Sol0, flux, elem, dt, stage);
#if OUTPUT_SUBSTEPS /* 4 */
            putout(Sol, t, *tout , step, 0, ud.file_name, "Sol", 1);
#endif
            
            /* explicit part of Euler backward gravity over half time step */
            euler_backward_gravity(Sol, (const MPV*)mpv, elem, 0.5*dt);
            Set_Explicit_Boundary_Data(Sol, elem, mpv, 1);
#if OUTPUT_SUBSTEPS /* 5 */
            putout(Sol, t, *tout , step, 0, ud.file_name, "Sol", 1);
#endif
            
            /* divergence-controlled advective fluxes at the half time level */
            recompute_advective_fluxes(flux, (const ConsVars*)Sol, elem);
            store_advective_fluxes(adv_flux, (const ConsVars**)flux, elem);
            flux_correction(flux, adv_flux_diff, buoy, elem, Sol, Sol0, t, dt, ud.implicitness, step);                
            update_advective_fluxes(flux, (const VectorField*)adv_flux, elem, node, dt);    
#if OUTPUT_SUBSTEPS  /* 6 */
            putout(Sol, t, *tout , step, 0, ud.file_name, "Sol", 1);
#endif
            
            /* reset Sol to state at beginning of time step */
            ConsVars_set(Sol, Sol0, elem->nc);
            
            /* ----------------------------------------------------------------------- */
            /* time step starts here                                                   */
            /* ----------------------------------------------------------------------- */
            
            printf("\nCoriolis 1 ---------------------------------- \n");                
            Explicit_Coriolis(Sol, elem, 0.5*dt); 
            Set_Explicit_Boundary_Data(Sol, elem, mpv, 0);
            
#if OUTPUT_SUBSTEPS  /* 7 */
            putout(Sol, t, *tout , step, 0, ud.file_name, "Sol", 1);
#endif
            
            /* explicit EULER half time step for gravity and pressure gradient */ 
            euler_forward_non_advective(Sol, (const MPV*)mpv, elem, node, 0.5*dt); 
            
#if OUTPUT_SUBSTEPS  /* 8 */
            putout(Sol, t, *tout , step, 0, ud.file_name, "Sol", 1);
#endif
            
            /* explicit full time step advection using div-controlled advective fluxes */
            printf("\nnonlinear advection ------------------------------ \n");
            
            substep = 0;
            VectorField_setzero(buoy, elem->nc);
            
            stage = 0;
            for(Split = 0; Split < elem->ndim; Split++) {
                lambda = 0.5*dt/elem->dx;
                Explicit_step_and_flux(Sol, flux[Split], lambda, elem->nc, Split, stage, FLUX_EXTERNAL, WITH_MUSCL);
                substep++;
                
#if OUTPUT_SUBSTEPS_PREDICTOR  /* 9, 10 */
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
                Explicit_step_and_flux(Sol, flux[Split], lambda, elem->nc, Split, stage, FLUX_EXTERNAL, WITH_MUSCL);
                substep++;
                
#if OUTPUT_SUBSTEPS_PREDICTOR  /* 11, 12 */
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
#if OUTPUT_SUBSTEPS  /* 13 */
            putout(Sol, t, *tout , step, 0, ud.file_name, "Sol", 1);
#endif
            
            second_projection(Sol, mpv, (const ConsVars*)Sol0, elem, node, 1.0, t, dt);
            Set_Explicit_Boundary_Data(Sol, elem, mpv, 1);
#if OUTPUT_SUBSTEPS  /* 14 */
            putout(Sol, t, *tout , step, 0, ud.file_name, "Sol", 1);
#endif
            
            /* This placement is not according to Piotr's rules ... */
            update_SI_MIDPT_buoyancy(Sol, (const ConsVars**)flux, mpv, elem, 0.5*dt);
#if OUTPUT_SUBSTEPS  /* 15 */
            putout(Sol, t, *tout , step, 0, ud.file_name, "Sol", 1);
#endif
            
            /* Bracketing Coriolis in Strang Split fashion around full advection-correction,
             that is after second projection, seems slightly better (see below) */
            printf("\nCoriolis 2 ---------------------------------- \n");                
            Explicit_Coriolis(Sol, elem, 0.5*dt); 
            
            if (ud.is_compressible) {
                adjust_pi_cells(mpv, Sol, elem); 
            }
            Set_Explicit_Boundary_Data(Sol, elem, mpv, 0);
            
#if OUTPUT_SUBSTEPS  /* 16 */
            putout(Sol, t, *tout , step, 0, ud.file_name, "Sol", 1);
#endif
            
            if (ud.absorber) 
            {
                Absorber(Sol, t, dt); 
                Set_Explicit_Boundary_Data(Sol, elem, mpv, 0);
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
