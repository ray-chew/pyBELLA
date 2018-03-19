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
    extern VectorField* adv_flux;
    extern ConsVars* flux[3];
    
    TimeStepInfo dt_info;
	const double* tout = ud.tout;
	int output_switch = 0;
                             	
	/* data allocation and initialization */
	Data_init();

    set_wall_massflux(bdry, Sol, elem);
    Set_Explicit_Boundary_Data(Sol, elem);
    /* This pre-projection is beneficial for getting a nodal pressure that
       is free of the multipole perturbations due to basic divergence errors
       that come simply from the divergence approximation on the nodal grid.
    second_projection(Sol, mpv, (const ConsVars*)Sol0, elem, node, 1.0, 0.0, 1.0);
    cell_pressure_to_nodal_pressure(mpv, elem, node);
     */
    ud.compressibility = compressibility(0);
        
	if(ud.write_file == ON) 
        putout(Sol, ud.file_name, "Sol", 1);
        
    ConsVars_set(Sol0, Sol, elem->nc);
    	    
    /* generate divergence-controlled initial data  */
    dt_info.time_step_switch = 0;

    ConsVars_set(Sol0, Sol, elem->nc);

	/* Main loop over the sequence of time values of tout */
	while(t < *tout && step < ud.stepmax) {
		
		/* start major time stepping loop */
		while( t < *tout && step < ud.stepmax ) { 
			            
            /* initialize fluxes in preparation of explicit predictor */
            ConsVars_setzero(flux[0], elem->nfx);
            if(elem->ndim > 1) ConsVars_setzero(flux[1], elem->nfy);
            if(elem->ndim > 2) ConsVars_setzero(flux[2], elem->nfz);            

            /* Timestep computation */
            dynamic_timestep(&dt_info, mpv, Sol, t, *tout, elem, step);
            dt = dt_info.time_step;
			            
			set_wall_massflux(bdry, Sol0, elem);
                       
            /* ======================================================================= */
            /* Semi-implicit discretization of non-advective terms aka EULAG           */
            /* ======================================================================= */
            
            reset_Y_perturbation(Sol, (const MPV*)mpv, elem);
            ConsVars_set(Sol0, Sol, elem->nc);            
                        
            printf("\n\n-----------------------------------------------------------------------------------------");
            printf("\nhalf-time prediction of advective flux");
            printf("\n-----------------------------------------------------------------------------------------\n");
                                    
            /* First order splitting for Corilis - just for the advection flux prediction */
            Explicit_Coriolis(Sol, elem, 0.5*dt); 
            
            /* explicit advection half time step preparing advection flux calculation 
             advect(Sol, flux, adv_flux, 0.5*dt, elem, FLUX_INTERNAL, WITH_MUSCL, SINGLE_STRANG_SWEEP, step%2);
             advect(Sol, flux, adv_flux, 0.5*dt, elem, FLUX_INTERNAL, WITH_MUSCL, DOUBLE_STRANG_SWEEP, step%2);
             advect(Sol, flux, adv_flux, 0.5*dt, elem, FLUX_EXTERNAL, WITH_MUSCL, DOUBLE_STRANG_SWEEP, step%2);
             - symmetrized splitting instead of simple splitting does not improve vortex symmetry  
             */
            recompute_advective_fluxes(flux, (const ConsVars*)Sol, elem);
            advect(Sol, flux, adv_flux, 0.5*dt, elem, FLUX_EXTERNAL, WITH_MUSCL, SINGLE_STRANG_SWEEP, step%2);

            /* explicit part of Euler backward gravity over half time step */
            euler_backward_gravity(Sol, (const MPV*)mpv, elem, 0.5*dt);
            
            /* divergence-controlled advective fluxes at the half time level */
            recompute_advective_fluxes(flux, (const ConsVars*)Sol, elem);
            store_advective_fluxes(adv_flux, (const ConsVars**)flux, elem);
            flux_correction(flux, elem, Sol, Sol0, t, dt, step);            
#ifndef CORRECT_FLUX_RIGHT_AWAY
            update_advective_fluxes(flux, (const VectorField*)adv_flux, elem, node, dt);    
#endif
            ConsVars_set(Sol, Sol0, elem->nc);
            // if (step == 0) cell_pressure_to_nodal_pressure(mpv, elem, node);
            // if (1) cell_pressure_to_nodal_pressure(mpv, elem, node);

            printf("\n\n-----------------------------------------------------------------------------------------");
            printf("\nfull time step with predicted advective flux");
            printf("\n-----------------------------------------------------------------------------------------\n");
            
            /* Strang splitting for Coriolis, first step */
            Explicit_Coriolis(Sol, elem, 0.5*dt); 
                        
            /* explicit EULER half time step for gravity and pressure gradient */ 
            euler_forward_non_advective(Sol, (const MPV*)mpv, elem, node, 0.5*dt); 
                        
            /* explicit full time step advection using div-controlled advective fluxes */
            advect(Sol, flux, adv_flux, dt, elem, FLUX_EXTERNAL, WITH_MUSCL, DOUBLE_STRANG_SWEEP, step%2);
                        
            /* implicit EULER half time step for gravity and pressure gradient */ 
            euler_backward_gravity(Sol, mpv, elem, 0.5*dt);
            second_projection(Sol, mpv, (const ConsVars*)Sol0, elem, node, 1.0, t, dt);
            
            /* auxiliary buoyancy update by vertical advection of background stratification */
            update_SI_MIDPT_buoyancy(Sol, (const ConsVars**)flux, mpv, elem, 0.5*dt);
            
            /* Strang splitting for Coriolis, second step */
            Explicit_Coriolis(Sol, elem, 0.5*dt); 
            
            if (ud.is_compressible) {
                adjust_pi_cells(mpv, Sol, elem); 
            }
                        
            if (ud.absorber) {
                Absorber(Sol, (const ElemSpaceDiscr*)elem, (const double)t, (const double)dt); 
            }                                    
            
            
			t += dt;
			step++;
			            			
            ud.compressibility = compressibility(t);
            
			if((ud.write_file == ON && (step % ud.write_file_period  == 0)) || output_switch) 
				putout(Sol, ud.file_name, "Sol", 1);
			
            if((ud.write_stdout == ON && step % ud.write_stdout_period == 0) || output_switch) {
                printf("\n############################################################################################");
                printf("\nstep %d done,  t=%f,  dt=%f,  cfl=%f, cfl_ac=%f, cfl_adv=%f", step, t, dt, dt_info.cfl, dt_info.cfl_ac, dt_info.cfl_adv);
                printf("\n############################################################################################\n");
            }
		}  
        
		if(ud.write_file == ON) putout(Sol, ud.file_name, "Sol", 1);
		tout++;
	}
	
	/* data release */
	Data_free();
	    
	return(1);
}
