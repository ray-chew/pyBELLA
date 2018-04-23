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
    extern double* force[3];
    extern ConsVars* flux[3];
        
    TimeStepInfo dt_info;
	const double* tout = ud.tout;
	int output_switch = 0;
                             	
	/* data allocation and initialization */
	Data_init();
    
    ud.nonhydrostasy   = nonhydrostasy(0);
    ud.compressibility = compressibility(0);
    

#ifdef FORCES_UNDER_OPSPLIT
    /* FORCES_UNDER_OPSPLIT not yet implemented for implicit gravity */
    assert(ud.g_ref < sqrt(DBL_EPSILON));  
#endif

    set_wall_massflux(bdry, Sol, elem);
    Set_Explicit_Boundary_Data(Sol, elem);
    /* This pre-projection is beneficial for getting a nodal pressure that
       is free of the multipole perturbations due to basic divergence errors
       that come simply from the divergence approximation on the nodal grid.
     */
    euler_backward_gravity(Sol, mpv, elem, 5.0);
    second_projection(Sol, mpv, (const ConsVars*)Sol0, elem, node, 0.0, 10.0);
    cell_pressure_to_nodal_pressure(mpv, elem, node);

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

            double flux_weight_old, flux_weight_new;
            
            /* initialize fluxes in preparation of explicit predictor */
            ConsVars_setzero(flux[0], elem->nfx);
            if(elem->ndim > 1) ConsVars_setzero(flux[1], elem->nfy);
            if(elem->ndim > 2) ConsVars_setzero(flux[2], elem->nfz);            

            /* Timestep computation */
            dynamic_timestep(&dt_info, mpv, Sol, t, *tout, elem, step);
            dt = dt_info.time_step;
			            
			set_wall_massflux(bdry, Sol0, elem);
                       
            /* ======================================================================= */
            /* Semi-implicit discretization of non-advective terms a la EULAG          */
            /* ======================================================================= */
            
            reset_Y_perturbation(Sol, (const MPV*)mpv, elem); // 4
            ConsVars_set(Sol0, Sol, elem->nc);            
                        
            printf("\n\n-----------------------------------------------------------------------------------------");
            printf("\nhalf-time prediction of advective flux");
            printf("\n-----------------------------------------------------------------------------------------\n");
                                    
            /* First order splitting for Corilis - just for the advection flux prediction */
            Explicit_Coriolis(Sol, elem, 0.5*dt);  // 5
            
            /* explicit advection half time step preparing advection flux calculation 
             advect(Sol, flux, force, 0.5*dt, elem, FLUX_INTERNAL, WITHOUT_FORCES, WITH_MUSCL, SINGLE_STRANG_SWEEP, step%2);
             advect(Sol, flux, force, 0.5*dt, elem, FLUX_INTERNAL, WITHOUT_FORCES, WITH_MUSCL, DOUBLE_STRANG_SWEEP, step%2);
             advect(Sol, flux, force, 0.5*dt, elem, FLUX_EXTERNAL, WITHOUT_FORCES, WITH_MUSCL, DOUBLE_STRANG_SWEEP, step%2);
             - symmetrized splitting instead of simple splitting does not improve vortex symmetry  
             */
            flux_weight_old = 0.0;
            flux_weight_new = 1.0;

#ifdef FULL_TIME_STEP_ADV_PREDICTOR
            recompute_advective_fluxes(flux, (const ConsVars*)Sol, elem, flux_weight_old, flux_weight_new);
            advect(Sol, flux, force, dt, elem, FLUX_EXTERNAL, WITHOUT_FORCES, WITH_MUSCL, SINGLE_STRANG_SWEEP, step%2);

            /* explicit part of Euler backward gravity over half time step */
            euler_backward_gravity(Sol, (const MPV*)mpv, elem, dt);
            
            /* divergence-controlled advective fluxes at the half time level */
            flux_weight_old = 0.5;
            flux_weight_new = 0.5;
            recompute_advective_fluxes(flux, (const ConsVars*)Sol, elem, flux_weight_old, flux_weight_new);
            flux_weight_old = 0.0;
            flux_weight_new = 1.0;
#else
            recompute_advective_fluxes(flux, (const ConsVars*)Sol, elem, flux_weight_old, flux_weight_new);
            advect(Sol, flux, force, 0.5*dt, elem, FLUX_EXTERNAL, WITHOUT_FORCES, WITH_MUSCL, SINGLE_STRANG_SWEEP, step%2); // 6

            /* explicit part of Euler backward gravity over half time step */
            euler_backward_gravity(Sol, (const MPV*)mpv, elem, 0.5*dt); // 7
            
            /* divergence-controlled advective fluxes at the half time level */
            recompute_advective_fluxes(flux, (const ConsVars*)Sol, elem, flux_weight_old, flux_weight_new);
#endif
            flux_correction(flux, Sol, Sol0, elem, node, t, dt, step);   // 8         

            ConsVars_set(Sol, Sol0, elem->nc);
            // if (step == 0) cell_pressure_to_nodal_pressure(mpv, elem, node);
            // if (1) cell_pressure_to_nodal_pressure(mpv, elem, node);

            printf("\n\n-----------------------------------------------------------------------------------------");
            printf("\nfull time step with predicted advective flux");
            printf("\n-----------------------------------------------------------------------------------------\n");
            
            /* Strang splitting for Coriolis, first step */
            Explicit_Coriolis(Sol, elem, 0.5*dt);  // 9
                        
            /* explicit EULER half time step for gravity and pressure gradient */ 
#ifdef FORCES_UNDER_OPSPLIT
            pressure_gradient_forces(force, Sol, (const MPV*)mpv, elem, node);
#else 
            euler_forward_non_advective(Sol, (const MPV*)mpv, elem, node, 0.5*dt, WITH_PRESSURE); // 10
#endif
                        
            /* explicit full time step advection using div-controlled advective fluxes */
            advect(Sol, flux, force, dt, elem, FLUX_EXTERNAL, WITH_FORCES, WITH_MUSCL, DOUBLE_STRANG_SWEEP, step%2); // 11
                        
            /* implicit EULER half time step for gravity and pressure gradient */ 
            euler_backward_gravity(Sol, mpv, elem, 0.5*dt); // 12
            second_projection(Sol, mpv, (const ConsVars*)Sol0, elem, node, t, dt); // 13
            
            /* auxiliary buoyancy update by vertical advection of background stratification */
            update_SI_MIDPT_buoyancy(Sol, (const ConsVars**)flux, mpv, elem, 0.5*dt); // 14
            
            /* Strang splitting for Coriolis, second step */
            Explicit_Coriolis(Sol, elem, 0.5*dt);  // 15
            
            if (ud.is_compressible) {
                adjust_pi_cells(mpv, Sol, elem); 
            }
                        
            if (ud.absorber) {
                Absorber(Sol, (const ElemSpaceDiscr*)elem, (const double)t, (const double)dt); 
            }                                    
            
            
			t += dt;
			step++;
			            			
            ud.nonhydrostasy   = nonhydrostasy(t);
            ud.compressibility = compressibility(t);
            
			if((ud.write_file == ON && (step % ud.write_file_period  == 0)) || output_switch) 
				putout(Sol, ud.file_name, "Sol", elem, node, 1);
			
            if((ud.write_stdout == ON && step % ud.write_stdout_period == 0) || output_switch) {
                printf("\n############################################################################################");
                printf("\nstep %d done,  t=%f,  dt=%f,  cfl=%f, cfl_ac=%f, cfl_adv=%f", step, t, dt, dt_info.cfl, dt_info.cfl_ac, dt_info.cfl_adv);
                printf("\n############################################################################################\n");
            }
		}  
        
		if(ud.write_file == ON) 
            putout(Sol, ud.file_name, "Sol", elem, node, 1);
		tout++;
	}
	
	/* data release */
	Data_free();
	    
	return(1);
}
