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
#include "molecular_transport.h"
#include "time.h"

/* ============================================================================= */

int main( void )
{
	/* low Mach */
	extern MPV* mpv;
	
	// /* User data */
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
    extern double* diss;

    TimeStepInfo dt_info;
	const double* tout = ud.tout;
	int output_switch = 0;

    // FILE *pfluxfile = NULL;
    char fn[120], fieldname[90];
    int swtch = 0;
        
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
        putout(Sol, ud.file_name, "Sol", elem, node, 1, 0, "ic");
            	    
    /* generate divergence-controlled initial data  */
    dt_info.time_step_switch = 0;

    time_t tic = time(NULL);
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

            FILE *tmpfile = NULL;
            int tmp_step = 0;
            printf("current tmpstep = %d", tmp_step);
            char fn[120], fieldname[90];
            if (step < 10) {
                sprintf(fn, "%s/flux_x/rhoY_00%d.hdf", ud.file_name, tmp_step);
            } else if(step < 100) {
                sprintf(fn, "%s/flux_x/rhoY_0%d.hdf", ud.file_name, tmp_step);
            } else {
                sprintf(fn, "%s/flux_x/rhoY_%d.hdf", ud.file_name, tmp_step);
            }
            sprintf(fieldname, "rhoY");
            WriteHDF(tmpfile, elem->icx, elem->icy, elem->icz, elem->ndim, Sol->rhoY, fn, fieldname);
            
            if (step < 10) {
                sprintf(fn, "%s/flux_x/rhou_00%d.hdf", ud.file_name, tmp_step);
            } else if(step < 100) {
                sprintf(fn, "%s/flux_x/rhou_0%d.hdf", ud.file_name, tmp_step);
            } else {
                sprintf(fn, "%s/flux_x/rhou_%d.hdf", ud.file_name, tmp_step);
            }
            sprintf(fieldname, "rhou");
            WriteHDF(tmpfile, elem->icy, elem->icx, elem->icz, elem->ndim, Sol->rhou, fn, fieldname);

            if (step < 10) {
                sprintf(fn, "%s/flux_x/rho_00%d.hdf", ud.file_name, tmp_step);
            } else if(step < 100) {
                sprintf(fn, "%s/flux_x/rho_0%d.hdf", ud.file_name, tmp_step);
            } else {
                sprintf(fn, "%s/flux_x/rho_%d.hdf", ud.file_name, tmp_step);
            }
            sprintf(fieldname, "rho");
            WriteHDF(tmpfile, elem->icy, elem->icx, elem->icz, elem->ndim, Sol->rho, fn, fieldname);

            if (swtch == 1){
                int tmp_step = 1;
                FILE *pfluxfile = NULL;
                if (step < 10) {
                    sprintf(fn, "%s/flux_x/rhoYu_00%d_before_advect.hdf", ud.file_name, tmp_step);
                } else if(step < 100) {
                    sprintf(fn, "%s/flux_x/rhoYu_0%d_before_advect.hdf", ud.file_name, tmp_step);
                } else {
                    sprintf(fn, "%s/flux_x/rhoYu_%d_before_advect.hdf", ud.file_name, tmp_step);
                }
                sprintf(fieldname, "rhoYu");
                WriteHDF(pfluxfile, elem->ifx, elem->icy, elem->icz, elem->ndim, flux[0]->rhoY, fn, fieldname);

                if (step < 10) {
                    sprintf(fn, "%s/flux_y/rhoYv_00%d_before_advect.hdf", ud.file_name, tmp_step);
                } else if(step < 100) {
                    sprintf(fn, "%s/flux_y/rhoYv_0%d_before_advect.hdf", ud.file_name, tmp_step);
                } else {
                    sprintf(fn, "%s/flux_y/rhoYv_%d_before_advect.hdf", ud.file_name, tmp_step);
                }
                sprintf(fieldname, "rhoYv");
                WriteHDF(pfluxfile, elem->ify, elem->icx, elem->icz, elem->ndim, flux[1]->rhoY, fn, fieldname);
                tmp_step++;
                // swtch = 0;
            }
            if (swtch == 1){
                putout(Sol, ud.file_name, "Sol", elem, node, 1, step, "before_advect");
                // swtch = 0;
            }            
            advect(Sol, flux, Sol0, 0.5*dt, elem, FLUX_EXTERNAL, WITH_MUSCL, DOUBLE_STRANG_SWEEP, ud.advec_time_integrator, step%2);

            if (swtch == 1){
                putout(Sol, ud.file_name, "Sol", elem, node, 1, step, "after_advect");
                // swtch = 0;
            }
          
            /* divergence-controlled advective fluxes at the half time level */
            for (int nn=0; nn<node->nc; nn++) mpv->p2_nodes0[nn] = mpv->p2_nodes[nn];
            euler_backward_non_advective_expl_part(Sol, (const MPV*)mpv, elem, 0.5*dt); 

            if (swtch == 1){
                putout(Sol, ud.file_name, "Sol", elem, node, 1, step, "after_ebnaexp");
                // swtch = 0;
            }        
            euler_backward_non_advective_impl_part(Sol, mpv, (const ConsVars*)Sol0, elem, node, t, 0.5*dt, 1.0, step);

            if (swtch == 1){
                putout(Sol, ud.file_name, "Sol", elem, node, 1, step, "after_ebnaimp");
                // swtch = 0;
            }
            recompute_advective_fluxes(flux, (const ConsVars*)Sol, bdry, elem, 0.5*dt);
            for (int nn=0; nn<node->nc; nn++) mpv->p2_nodes[nn] = mpv->p2_nodes0[nn];

            ConsVars_set(Sol, Sol0, elem->nc);
            if (swtch == 1){
                putout(Sol, ud.file_name, "Sol", elem, node, 1, step, "after_half_step");
                // swtch = 0;
            }

            if (swtch == 1){
                int tmp_step = 5;
                FILE *pfluxfile = NULL;
                if (step < 10) {
                    sprintf(fn, "%s/flux_x/rhoYu_00%d_after_half_step.hdf", ud.file_name, tmp_step);
                } else if(step < 100) {
                    sprintf(fn, "%s/flux_x/rhoYu_0%d_after_half_step.hdf", ud.file_name, tmp_step);
                } else {
                    sprintf(fn, "%s/flux_x/rhoYu_%d_after_half_step.hdf", ud.file_name, tmp_step);
                }
                sprintf(fieldname, "rhoYu");
                WriteHDF(pfluxfile, elem->ifx, elem->icy, elem->icz, elem->ndim, flux[0]->rhoY, fn, fieldname);

                if (step < 10) {
                    sprintf(fn, "%s/flux_y/rhoYv_00%d_after_half_step.hdf", ud.file_name, tmp_step);
                } else if(step < 100) {
                    sprintf(fn, "%s/flux_y/rhoYv_0%d_after_half_step.hdf", ud.file_name, tmp_step);
                } else {
                    sprintf(fn, "%s/flux_y/rhoYv_%d_after_half_step.hdf", ud.file_name, tmp_step);
                }
                sprintf(fieldname, "rhoYv");
                WriteHDF(pfluxfile, elem->ify, elem->icx, elem->icz, elem->ndim, flux[1]->rhoY, fn, fieldname);
                tmp_step++;
                // swtch = 0;
            }

            if (swtch == 1){
                // printf("### output temporary Sol ###, swtch = %d", swtch);
                // putout(Sol, ud.file_name, "Sol", elem, node, 1);

                // printf('\n\n#### output flux arrays ####');
                FILE *pfluxfile = NULL;
                if (step < 10) {
                    sprintf(fn, "%s/flux_x/rhoYu_00%d.hdf", ud.file_name, tmp_step);
                } else if(step < 100) {
                    sprintf(fn, "%s/flux_x/rhoYu_0%d.hdf", ud.file_name, tmp_step);
                } else {
                    sprintf(fn, "%s/flux_x/rhoYu_%d.hdf", ud.file_name, tmp_step);
                }
                sprintf(fieldname, "rhoYu");
                WriteHDF(pfluxfile, elem->ifx, elem->icy, elem->icz, elem->ndim, flux[0]->rhoY, fn, fieldname);

                if (step < 10) {
                    sprintf(fn, "%s/flux_y/rhoYv_00%d.hdf", ud.file_name, tmp_step);
                } else if(step < 100) {
                    sprintf(fn, "%s/flux_y/rhoYv_0%d.hdf", ud.file_name, tmp_step);
                } else {
                    sprintf(fn, "%s/flux_y/rhoYv_%d.hdf", ud.file_name, tmp_step);
                }
                sprintf(fieldname, "rhoYv");
                WriteHDF(pfluxfile, elem->ify, elem->icx, elem->icz, elem->ndim, flux[1]->rhoY, fn, fieldname);
                tmp_step++;
            }

            printf("\n\n-----------------------------------------------------------------------------------------");
            printf("\nfull time step with predicted advective flux");
            printf("\n-----------------------------------------------------------------------------------------\n");

            /* add effect of heat conduction and diffusion as computed in euler_backward_...() to rhoY */
            diss_to_rhoY(Sol, diss, elem, node);

            /* explicit EULER half time step for gravity and pressure gradient */ 
            euler_forward_non_advective(Sol, mpv, (const ConsVars*)Sol0, elem, node, (dt_factor-0.5)*dt, WITH_PRESSURE);

            if (swtch == 1){
                putout(Sol, ud.file_name, "Sol", elem, node, 1, step, "after_efna");
                // swtch = 0;
            }
            /* explicit full time step advection using div-controlled advective fluxes */
            advect(Sol, flux, Sol0, dt_factor*dt, elem, FLUX_EXTERNAL, WITH_MUSCL, DOUBLE_STRANG_SWEEP, ud.advec_time_integrator, step%2); 
            if (swtch == 1){
                putout(Sol, ud.file_name, "Sol", elem, node, 1, step, "after_full_advect");
                // swtch = 0;
            }

            /* add effect of heat conduction and diffusion as computed in euler_backward_...() to rhoY */
            diss_to_rhoY(Sol, diss, elem, node);        

            /* implicit EULER half time step for gravity and pressure gradient */ 
            euler_backward_non_advective_expl_part(Sol, mpv, elem, 0.5*dt);
            if (swtch == 1){
                putout(Sol, ud.file_name, "Sol", elem, node, 1, step, "after_full_ebnaexp");
                // swtch = 0;
            }
            euler_backward_non_advective_impl_part(Sol, mpv, (const ConsVars*)Sol0, elem, node, t, 0.5*dt, 2.0, step);

            synchronize_variables(mpv, Sol, elem, node, ud.synchronize_nodal_pressure);
           
            if (ud.absorber) {
                Absorber(Sol, (const ElemSpaceDiscr*)elem, (const double)t, (const double)dt); 
            }                                    
                        
            if (ud.n_time_series > 0) {
                store_time_series_entry(Sol, elem, step);
            }

            putout(Sol, ud.file_name, "Sol", elem, node, 1, step, "after_full_step");


            t += dt;
            step++;
            dt_factor = 1.0;


            // if((ud.write_file == ON && (step % ud.write_file_period  == 0)) || output_switch) 
            //     putout(Sol, ud.file_name, "Sol", elem, node, 1);

            // if((ud.write_stdout == ON && step % ud.write_stdout_period == 0) || output_switch) {
            printf("\n############################################################################################");
            printf("\nstep %d done,  t=%f,  dt=%.16f,  cfl=%f, cfl_ac=%f, cfl_adv=%f, dt_factor=%f", step, t, dt, dt_info.cfl, dt_info.cfl_ac, dt_info.cfl_adv, dt_factor);
            printf("\n############################################################################################\n");
            // }            
		}  
        
        // if(ud.write_file == ON) {
        //     putout(Sol, ud.file_name, "Sol", elem, node, 1);
        //     dt_info.time_step_switch = 0;
        // }
		tout++;
	}

    time_t toc = time(NULL);
    printf("\n################");
    printf("\nTime spent = %ds", (toc - tic));
    printf("\n################");
    printf("\n");
	
	/* data release */
	Data_free();
    
    if (ud.n_time_series > 0) {
        close_time_series();
    }

	return(1);
}
