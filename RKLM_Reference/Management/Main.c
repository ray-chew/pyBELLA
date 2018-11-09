/*******************************************************************************
 File:   Main.c
 *******************************************************************************/
#include "Common.h"
#include <stdio.h>
#include <float.h>

#ifdef MACPROFILE
#include <profiler.h> 
#include <math.h>
#include <stdlib.h>
#include <string.h>
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
#include "second_projection.h"
#include "mpv.h"
#include "boundary.h"
#include "numerical_flux.h"
#include "enumerator.h"

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
    extern double* force[3];
    extern ConsVars* flux[3];
        
    TimeStepInfo dt_info;
	const double* tout = ud.tout;
	int output_switch = 0;
        
	/* data allocation and initialization */
	Data_init();
        
    if (ud.n_time_series > 0) {
        initialize_time_series();
    }

    ud.nonhydrostasy   = nonhydrostasy(0);
    ud.compressibility = compressibility(0);
    
    set_wall_massflux(bdry, Sol, elem);
    Set_Explicit_Boundary_Data(Sol, elem);
    /* This pre-projection is beneficial for getting a nodal pressure that
       is free of the multipole perturbations due to basic divergence errors
       that come simply from the divergence approximation on the nodal grid.
     */
    if (ud.initial_projection == CORRECT) {
        assert(0);  /* MAKE SURE THE TIME STEP SIZE PARAMETERS ARE OK IN THE FOLLOWING CALLS!! */
        euler_backward_non_advective_expl_part(Sol, mpv, elem, 5.0);
        euler_backward_non_advective_impl_part(Sol, mpv, (const ConsVars*)Sol0, elem, node, 0.0, 10.0);
    }

	if(ud.write_file == ON) 
        putout(Sol, ud.file_name, "Sol", elem, node, 1);
        
    ConsVars_set(Sol0, Sol, elem->nc);
    	    
    /* generate divergence-controlled initial data  */
    dt_info.time_step_switch = 0;

    ConsVars_set(Sol0, Sol, elem->nc);

	/* Main loop over the sequence of time values of tout */
	while(t < *tout && step < ud.stepmax) {
		
		/* start major time stepping loop */
		while( t < *tout && step < ud.stepmax ) { 
            
            /* Timestep computation */
            dynamic_timestep(&dt_info, mpv, Sol, t, *tout, elem, step);
            dt = dt_info.time_step;

            /* model and numerics controls */
            ud.nonhydrostasy   = nonhydrostasy(t);
            ud.compressibility = compressibility(t);
            ud.acoustic_order  = acoustic_order(t, dt);

            /* initialize fluxes in preparation of explicit predictor */
            ConsVars_setzero(flux[0], elem->nfx);
            if(elem->ndim > 1) ConsVars_setzero(flux[1], elem->nfy);
            if(elem->ndim > 2) ConsVars_setzero(flux[2], elem->nfz);            
			            
			set_wall_massflux(bdry, Sol0, elem);
                       
            /* ======================================================================= */
            /* Semi-implicit discretization of non-advective terms a la EULAG          */
            /* ======================================================================= */
            
            ConsVars_set(Sol0, Sol, elem->nc);            
                        
            printf("\n\n-----------------------------------------------------------------------------------------");
            printf("\nhalf-time prediction of advective flux");
            printf("\n-----------------------------------------------------------------------------------------\n");
                                          
#ifdef CORIOLIS_EXPLICIT
            /* First order splitting for Corilis - just for the advection flux prediction */ 
             Explicit_Coriolis(Sol, elem, 0.5*dt);
#endif
            
            recompute_advective_fluxes(flux, (const ConsVars*)Sol, elem);
#ifdef ADVECTION
            advect(Sol, flux, force, 0.5*dt, elem, FLUX_EXTERNAL, WITH_MUSCL, SINGLE_STRANG_SWEEP, step%2);
            // reset_rhoY(Sol, Sol0, elem);
#endif
            /* divergence-controlled advective fluxes at the half time level */
            euler_backward_non_advective_expl_part(Sol, (const MPV*)mpv, elem, 0.5*dt);            
            recompute_advective_fluxes(flux, (const ConsVars*)Sol, elem);
            flux_correction(flux, Sol, Sol0, elem, node, t, 0.5*dt, step);        

            ConsVars_set(Sol, Sol0, elem->nc);
            cell_pressure_to_nodal_pressure(mpv, elem, node, 2.0-ud.acoustic_order);
          
            printf("\n\n-----------------------------------------------------------------------------------------");
            printf("\nfull time step with predicted advective flux");
            printf("\n-----------------------------------------------------------------------------------------\n");
            
#ifdef CORIOLIS_EXPLICIT
            /* Strang splitting for Coriolis, first step */
             Explicit_Coriolis(Sol, elem, 0.5*dt);  
#endif
            
            /* explicit EULER half time step for gravity and pressure gradient */ 
            euler_forward_non_advective(Sol, mpv, (const ConsVars*)Sol0, elem, node, 0.5*dt, WITH_PRESSURE);

                        
#ifdef ADVECTION
            /* explicit full time step advection using div-controlled advective fluxes */
            advect(Sol, flux, force, 1.0*dt, elem, FLUX_EXTERNAL, WITH_MUSCL, DOUBLE_STRANG_SWEEP, step%2);
#endif
            
            /* implicit EULER half time step for gravity and pressure gradient */ 
            euler_backward_non_advective_expl_part(Sol, mpv, elem, 0.5*dt);
            euler_backward_non_advective_impl_part(Sol, mpv, (const ConsVars*)Sol0, elem, node, t, 0.5*dt);
                        
#ifdef CORIOLIS_EXPLICIT
            /* Strang splitting for Coriolis, second step */
             Explicit_Coriolis(Sol, elem, 0.5*dt);  
#endif
            
#if 0
            if((ud.write_file == ON && ((step+1) % ud.write_file_period  == 0)) || output_switch) 
                putout(Sol, ud.file_name, "Sol", elem, node, 1);
#endif
            
            synchronize_variables(mpv, Sol, elem, node, ud.synchronize_nodal_pressure);

                        
            if (ud.absorber) {
                Absorber(Sol, (const ElemSpaceDiscr*)elem, (const double)t, (const double)dt); 
            }                                    
                        
            if (ud.n_time_series > 0) {
                store_time_series_entry(Sol, elem, step);
            }            
            
            t += dt;
            step++;
#if 1
            if((ud.write_file == ON && (step % ud.write_file_period  == 0)) || output_switch) 
                putout(Sol, ud.file_name, "Sol", elem, node, 1);
#endif        
            if((ud.write_stdout == ON && step % ud.write_stdout_period == 0) || output_switch) {
                printf("\n############################################################################################");
                printf("\nstep %d done,  t=%f,  dt=%f,  cfl=%f, cfl_ac=%f, cfl_adv=%f", step, t, dt, dt_info.cfl, dt_info.cfl_ac, dt_info.cfl_adv);
                printf("\n############################################################################################\n");
            }
		}  
        
        if(ud.write_file == ON) {
            putout(Sol, ud.file_name, "Sol", elem, node, 1);
            dt_info.time_step_switch = 0;
        }
		tout++;
	}
	
	/* data release */
	Data_free();
	    
    if (ud.n_time_series > 0) {
        close_time_series();
    }

	return(1);
}
