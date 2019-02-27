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
    extern double* W0;
    extern enum Boolean W0_in_use;

    TimeStepInfo dt_info;
	const double* tout = ud.tout;
	int output_switch = 0;
        
    double dt_factor;
    
	/* data allocation and initialization */
	Data_init();
        
    /* An initial implicit Euler step in the second projection removes noise 
     in some cases. It effectively acts like a filter on the initial data. */
    dt_factor = (ud.initial_impl_Euler == CORRECT ? 0.5 : 1.0);

    if (ud.n_time_series > 0) {
        initialize_time_series();
    }

	if(ud.write_file == ON) 
        putout(Sol, ud.file_name, "Sol", elem, node, 1);
            	    
    /* generate divergence-controlled initial data  */
    dt_info.time_step_switch = 0;

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
			            
            ConsVars_set(Sol0, Sol, elem->nc);            
            set_wall_rhoYflux(bdry, Sol0, mpv, elem);
           
            /* ======================================================================= */
            /* Semi-implicit discretization of non-advective terms a la EULAG          */
            /* ======================================================================= */
            
            printf("\n\n-----------------------------------------------------------------------------------------");
            printf("\nhalf-time prediction of advective flux");
            printf("\n-----------------------------------------------------------------------------------------\n");
                                                      
            recompute_advective_fluxes(flux, (const ConsVars*)Sol, bdry, elem, 0.5*dt);
            advect(Sol, flux, Sol0, 0.5*dt, elem, FLUX_EXTERNAL, WITH_MUSCL, SINGLE_STRANG_SWEEP, ud.advec_time_integrator, step%2);

            /* divergence-controlled advective fluxes at the half time level */
            for (int nn=0; nn<node->nc; nn++) mpv->p2_nodes0[nn] = mpv->p2_nodes[nn];
            euler_backward_non_advective_expl_part(Sol, (const MPV*)mpv, elem, 0.5*dt); 
            euler_backward_non_advective_impl_part(Sol, mpv, (const ConsVars*)Sol0, elem, node, t, 0.5*dt, 1.0);
            recompute_advective_fluxes(flux, (const ConsVars*)Sol, bdry, elem, 0.5*dt);
            for (int nn=0; nn<node->nc; nn++) mpv->p2_nodes[nn] = mpv->p2_nodes0[nn];
            
            ConsVars_set(Sol, Sol0, elem->nc);
            
            printf("\n\n-----------------------------------------------------------------------------------------");
            printf("\nfull time step with predicted advective flux");
            printf("\n-----------------------------------------------------------------------------------------\n");

            /* explicit EULER half time step for gravity and pressure gradient */ 
            euler_forward_non_advective(Sol, mpv, (const ConsVars*)Sol0, elem, node, (dt_factor-0.5)*dt, WITH_PRESSURE);
                        
            /* explicit full time step advection using div-controlled advective fluxes */
            advect(Sol, flux, Sol0, dt_factor*dt, elem, FLUX_EXTERNAL, WITH_MUSCL, DOUBLE_STRANG_SWEEP, ud.advec_time_integrator, step%2); 
            
            /* implicit EULER half time step for gravity and pressure gradient */ 
            euler_backward_non_advective_expl_part(Sol, mpv, elem, 0.5*dt);
            euler_backward_non_advective_impl_part(Sol, mpv, (const ConsVars*)Sol0, elem, node, t, 0.5*dt, 2.0);
                                                
            synchronize_variables(mpv, Sol, elem, node, ud.synchronize_nodal_pressure);

            if (ud.absorber) {
                Absorber(Sol, (const ElemSpaceDiscr*)elem, (const double)t, (const double)dt); 
            }                                    
                        
            if (ud.n_time_series > 0) {
                store_time_series_entry(Sol, elem, step);
            }            
            
            t += dt;
            step++;
            dt_factor = 1.0;

            if((ud.write_file == ON && (step % ud.write_file_period  == 0)) || output_switch) 
                putout(Sol, ud.file_name, "Sol", elem, node, 1);

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
