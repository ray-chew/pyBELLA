/*******************************************************************************
 File:   userdata.c
 Author: Nicola
 Date:   Fri Feb 27 09:34:03 CET 1998
 *******************************************************************************/
#include <math.h>
#include <float.h>
#include <stdlib.h>

#include "Common.h"
#include "userdata.h"
#include "enumerator.h"
#include "time_discretization.h"
#include "error.h"
#include "variable.h"
#include "enum_bdry.h"
#include "math_own.h"
#include "space_discretization.h"
#include "thermodynamic.h"
#include "Eos.h"
#include "boundary.h"
#include "memory.h"
#include "Hydrostatics.h"

/*
 Monitor configurations:
 
 OpSplit:
 compressible case 
 unless otherwise listed, I tested options at fixed: 
       CFL           0.45, 
       grid          160x80, 
       t_end         10.0; 
       tol           1.0e-07, 
       limiters      NONE, 
       ud.p_extrapol 1.0
       ccenterweight 2.0
       HYDROSTATES_USING_dp2 off
 
 Option tests:
 - Diagonal Five Point in second projection (P2_DIAGONAL_FIVE_POINT = 1.0 in ProjectionType.h)
 - ud->p_flux_correction = CORRECT
 - ud->p_average = WRONG
 - Diagonal fivepoint (exact proj) or full ninepoint (approx proj) ; egal
 - DP_AVERAGE off in ProjectionType.h
 - DPDT_AVERAGED_IN_HELMHOLTZ off in ProjectionType.h
 
 Results:
 -- CFL 0.99                         -> very good
 -- dt_fixed = 0.004                 -> very good    (0.4 sec for cfl_ac = 1.11 / decreased tol to 1.0e-09) 
 -- CFL 0.99, dt_fixed 0.08, 320x160 -> very good    (kept tol at 1.0e-09)
 -- HYDROSTATES_USING_dp2 on 160x080 -> very good    (dp2_factor 0.5, CFL 0.99 & 0.45 & dt 0.004)  
 - diag-fivepoint with these options -> very good

 sound-proof case: 
 -- last option above, CFL 0.99      -> NO           (first projection doesnt converge in step 2 -- ??????????????)
 -- HYDROSTATES_USING_dp2 off        -> very good  
 
 PREV:   NONE
 NEXT:   TRAVELLING VORTEX
 
 */ 

double pressure_function(double r, double p0, double S0, double u_theta, double Msq, double Gamma);

double rho_function(double psi);


void User_Data_init(User_Data* ud) {
	
	int i;
	
	/* ================================================================================== */
	/* =====  PROBLEM SET UP  =========================================================== */
	/* ================================================================================== */
	
    /* Earth */	
	double grav  = 10.0;                               /* [m/s^2]                         */
	double omega = 0*2*PI*sin(0.25*PI)/(24.0*3600.0);    /*  [s^-1]                         */

	/* references for non-dimensionalization */
	double h_ref = 10000;            /* [m]                             */
	double t_ref = 100;              /* [s]                             */
	double T_ref = 300;              /* [K]                             */
	double p_ref = 8.61 * 10e4;           /* [Pa]                            */
	double u_ref = h_ref/t_ref;      /* Strouhal No == 1 always assumed */
		
	/* thermodynamics and chemistry */
	double R_gas = 287.4;            /* [J/kg/K]                        */
	double R_vap = 461.00;           /* [J/kg/K]                        */
    double Q_vap = 2.53e+06;         /* [J]                             */
	double gamma = 1.4;              /* dimensionless                   */

    double viscm  = 0.0;            /* [m^2/s]                         */
    double viscbm = 0.0;             /* [m^2/s]                         */
    double visct  = 0.0;             /* [m^2/s]                         */
    double viscbt = 0.0;             /* [m^2/s]                         */
    double cond   = 0.0;            /* [m^2/s]                         */

    /* references for non-dimensionalization */
    // double h_ref   = 10000;            /* [m]                             */
    // double t_ref   = 100;              /* [s]                             */
    // double T_ref   = 300;              /* [K]                             */
    // double p_ref   = 101325;           /* [Pa]                            */
    // double u_ref   = h_ref/t_ref;      /* Strouhal No == 1 always assumed */
    // double rho_ref = p_ref / (R_gas*T_ref); /* [kg/m^3]  */
    double rho_ref = p_ref / (R_gas * T_ref);
             
    double Nsq_ref = grav*1.3e-05;     /* [] */

    ud->h_ref       = h_ref;
    ud->t_ref       = t_ref;
    ud->T_ref       = T_ref;
    ud->p_ref       = p_ref;
    ud->rho_ref     = rho_ref;
    ud->u_ref       = u_ref;
    ud->Nsq_ref     = Nsq_ref;
    ud->g_ref       = grav;
    ud->gamm        = gamma;  
    ud->Rg_over_Rv  = R_gas/R_vap;  
    ud->Q           = Q_vap/(R_gas*T_ref);  
    
    /* number of advected species */
    ud->nspec       = NSPEC;  

    /*FULL_MOLECULAR_TRANSPORT, STRAKA_DIFFUSION_MODEL, NO_MOLECULAR_TRANSPORT */
    ud->mol_trans   = NO_MOLECULAR_TRANSPORT; 
    ud->viscm       = viscm  * t_ref/(h_ref*h_ref);
    ud->viscbm      = viscbm * t_ref/(h_ref*h_ref);
    ud->visct       = visct  * t_ref/(h_ref*h_ref);
    ud->viscbt      = viscbt * t_ref/(h_ref*h_ref);
    ud->cond        = cond * t_ref/(h_ref*h_ref*R_gas);

    /* Low Mach */
    ud->is_nonhydrostatic =  1;    /* 0: hydrostatic;  1: nonhydrostatic;  -1: transition (see nonhydrostasy()) */
    ud->is_compressible   =  1;    /* 0: psinc;  1: comp;  -1: psinc-comp-transition (see compressibility()) */
    ud->acoustic_timestep =  0;    /* advective time step -> 0;  acoustic time step -> 1; */
    ud->Msq =  u_ref*u_ref / (R_gas*T_ref);
	
	/* geo-stuff */
	for(i=0; i<3; i++) {
        ud->gravity_strength[i]  = 0.0;   /* corresponds to  M^2 / Fr^2  =  g href / R Tref */
        ud->coriolis_strength[i] = 0.0;   
    }  
	ud->gravity_strength[1]  = grav * h_ref / (R_gas * T_ref); 
	ud->coriolis_strength[0] = omega * t_ref; 
	ud->coriolis_strength[2] = omega * t_ref; 
	
	/* integer gravity indicator */
	for (i=0; i<3; i++){
		ud->i_gravity[i]  = 0;
		ud->i_coriolis[i] = 0;
		if (ud->gravity_strength[i] > ud->eps_Machine || i==1) {
			ud->i_gravity[i] = 1;
			ud->gravity_direction = i;
		}
		if (ud->coriolis_strength[i] > ud->eps_Machine) {
			ud->i_coriolis[i] = 1;
		}
	}

    /* low Froude */
    // ud->implicit_gravity_theta  = 0;
    // ud->implicit_gravity_press  = 0; /* should be on for compressible option */
    // ud->implicit_gravity_theta2 = 0;
    // ud->implicit_gravity_press2 = 0; /* should this, too, be on for compressible option ?  */
	    
	/* flow domain */
	ud->xmin = - 1.0;  
	ud->xmax =   1.0;  
	ud->ymin =   0.0;
	ud->ymax =   1.0; 
	ud->zmin = - 1.0;
	ud->zmax =   1.0;

	/* boundary/initial conditions */
	ud->wind_speed  =  0.0;            /* 0.1535; wind_speed * u_ref  = 10 m/s */ /*  BREAKING WAVE CHANGE */
	ud->wind_shear  = -0.0;              /* velocity in [u_ref/h_ref] */             
    ud->hill_shape  = AGNESI;            /* AGNESI, SCHLUTOW */
    ud->hill_height = 0.0 * 0.096447; /* 0.096447; */ /* hill_height * l_ref = 0.628319 km *//*  BREAKING WAVE CHANGE */
	ud->hill_length_scale = 0.1535;   /* hill_length * l_ref = 1.0 km */
	
	ud->bdrytype_min[0] = PERIODIC; /* DIRICHLET; PERIODIC; WALL; */
	ud->bdrytype_min[1] = WALL; /* SLANTED_WALL; */
	ud->bdrytype_min[2] = WALL;
	ud->bdrytype_max[0] = PERIODIC; /* DIRICHLET; PERIODIC; WALL; */  
	ud->bdrytype_max[1] = WALL;  
	ud->bdrytype_max[2] = WALL;
	
	ud->absorber = WRONG; /* CORRECT;  WRONG; */ /*  BREAKING WAVE CHANGE */
    ud->bottom_theta_bc = BOTTOM_BC_DEFAULT;
	
	/* ================================================================================== */
	/* =====  NUMERICS  ================================================================= */
	/* ================================================================================== */
	
    /* time discretization */
    ud->time_integrator      = SI_MIDPT;
    ud->advec_time_integrator = STRANG;
	ud->CFL                  = 1.0; /* 0.45; 0.9; 0.8; */
    ud->dtfixed0             = 17.0 / ud->t_ref;
    ud->dtfixed              = 17.0 / ud->t_ref; /* 0.0052; */ /*  0.004; */ 
    // ud->no_of_steps_to_CFL   = 1;
    // ud->no_of_steps_to_dtfixed = 1;
    
    set_time_integrator_parameters(ud);
    
	/* Grid and space discretization */
	ud->inx = 160+1; /* 641; 321; 161; 129; 81; */    
	ud->iny =  80+1; /* 321; 161;  81;  65; 41;  */
	ud->inz =  1;
	// ud->h    = MIN_own((ud->xmax-ud->xmin)/(ud->inx),MIN_own((ud->ymax-ud->ymin)/(ud->iny),(ud->zmax-ud->zmin)/(ud->inz)));
	
	/* explicit predictor step */
	/* Recovery */
	ud->recovery_order = SECOND;
	ud->limiter_type_scalars  = NONE; /*  RUPE; NONE; MONOTONIZED_CENTRAL; MINMOD; VANLEER; SWEBY_MUNZ; SUPERBEE; */
	ud->limiter_type_velocity = NONE; /*  RUPE; NONE; MONOTONIZED_CENTRAL; MINMOD; VANLEER; SWEBY_MUNZ; SUPERBEE; */
	
    /* first correction */
	// ud->p_flux_correction = WRONG; /* CORRECT, WRONG; */
    // if (ud->time_integrator == OP_SPLIT || ud->time_integrator == OP_SPLIT_MD_UPDATE) {
    //     ud->latw[0] = ud->latw[2] = 0.125; ud->latw[1] = 0.75; ud->p_extrapol = 1.0; 
    //     /* ud->latw[0] = ud->latw[2] = 0.125; ud->latw[1] = 0.75; ud->p_extrapol = 1.25; */
    //     /* ud->latw[0] = ud->latw[2] = 0.25; ud->latw[1] = 0.5; ud->p_extrapol = 1.0; */
    // } else {
    //     /* ud->latw[0] = ud->latw[2] = 0.125; ud->latw[1] = 0.75; ud->p_extrapol = 1.25;*/
    //     ud->latw[0] = ud->latw[2] = 0.25; ud->latw[1] = 0.5; ud->p_extrapol = 1.5; 
    // }
    
    /* parameters for SWEBY_MUNZ limiter family */
    ud->kp = 1.4;
	ud->kz = 1.4; /* Entro abused for velocity in split-step-aligned velocity ! */
	ud->km = 1.4;
	ud->kY = 1.4;
	ud->kZ = 1.4; /* 2.0 */
	
	ud->ncache =  300; /* (ud->inx+3); */
	
    /* linear solver-stuff */
    double tol                            = 1.e-8 * (ud->is_compressible == 1 ? 0.01 : 1.0);
    ud->flux_correction_precision         = tol;
    ud->flux_correction_local_precision   = tol;    /* 1.e-05 should be enough */
    ud->second_projection_precision       = tol;
    ud->second_projection_local_precision = tol;  /* 1.e-05 should be enough */
    ud->flux_correction_max_iterations    = 6000;
    ud->second_projection_max_iterations  = 6000;
    ud->initial_projection                = WRONG;   /* WRONG;  CORRECT; */
    ud->initial_impl_Euler                = CORRECT;   /* WRONG;  CORRECT; */
    
    ud->column_preconditioner             = CORRECT; /* WRONG; CORRECT; */
    ud->synchronize_nodal_pressure        = WRONG;   /* WRONG; CORRECT; */
    ud->synchronize_weight                = 0.0;    /* relevant only when prev. option is "CORRECT"
                                                     Should ultimately be a function of dt . */  

	/* numerics parameters */
	ud->eps_Machine = sqrt(DBL_EPSILON);
		
	/* ================================================================================== */
	/* =====  FLOW CONTROL  ============================================================= */
	/* ================================================================================== */

    /* output times */  
	ud->tout[0] =  1000.0 / ud->t_ref;      
	ud->tout[1] =  1000.0 / ud->t_ref;      
	// ud->tout[2] =  6.0;      
	// ud->tout[3] =  8.0;      
    // ud->tout[4] = 10.0;      
    // ud->tout[5] = 10.5;      
	// ud->tout[6] = -1.0;      

    ud->stepmax = 10000;
    
	ud->write_stdout = ON;
	ud->write_stdout_period = 1;
	ud->write_file = ON;
	ud->write_file_period = 10000;
	ud->file_format = HDF;
    
    // ud->write_history = 0;


    {
        char *OutputBaseFolder      = "/home/ray/git-projects/RKLM_Reference/RKLM_Reference/output_rising_bubble/";
        char *OutputFolderNamePsinc = "low_Mach_gravity_psinc";
        char *OutputFolderNameComp  = "low_Mach_gravity_comp";
        if (ud->is_compressible == 0) {
            sprintf(ud->file_name, "%s%s", OutputBaseFolder, OutputFolderNamePsinc);
        } else {
            sprintf(ud->file_name, "%s%s", OutputBaseFolder, OutputFolderNameComp);
        }
    }
}

/* ================================================================================== */

void Sol_initial(ConsVars* Sol, 
                 ConsVars* Sol0, 
                 MPV* mpv,
                 BDRY* bdry,
                 const ElemSpaceDiscr* elem,
                 const NodeSpaceDiscr* node) {
	
	extern Thermodynamic th;
	extern User_Data ud;
    extern double *Yinvbg;
    
	const double u0 = ud.wind_speed;
	const double v0 = 0.0;
	const double w0 = 0.0;
	const double delth = 2.0; /* pot temp perturbation in K */
	
    const double y0 = 0.2;
    const double r0 = 0.2;
    	
	const int icx = elem->icx;
	const int icy = elem->icy;
	const int icz = elem->icz;
	const int igy = elem->igy;
	const int igz = elem->igz;
	
    const int icxn = node->icx;
    const int icyn = node->icy;
    const int iczn = node->icz;
    
	int i, j, k, l, m, n;
	double x, y, z;
	double rho, u, v, w, p, rhoY;
	
    double g;
	    
	g = ud.gravity_strength[1];
    
    Hydrostatics_State(mpv, elem, node);
		    
    /* data in the bulk of the domain */
	for(k = igz; k < icz - igz; k++) {l = k * icx * icy; 
		z = elem->z[k];
                
        for(j = igy; j < icy - igy; j++) {m = l + j * icx;
            y = elem->y[j];
			
			for(i = 0; i < icx; i++) {n = m + i;
                double r;
                
				x = elem->x[i];
				
                u   = u0;
                v   = v0;
                w   = w0;

				r = sqrt(x*x + (y-y0)*(y-y0) + (elem->ndim==3)*z*z) / r0;
                p     = mpv->HydroState->p0[j];
                rhoY  = mpv->HydroState->rhoY0[j];
				rho  = rhoY / (stratification(y) + (r > 1.0 ? 0.0 : (delth/300)*cos(0.5*PI*r)*cos(0.5*PI*r)));
				
				Sol->rho[n]  = rho;
				Sol->rhou[n] = rho * u;
				Sol->rhov[n] = rho * v;
				Sol->rhow[n] = rho * w;
				Sol->rhoe[n] = rhoe(rho, u, v, w, p);
				Sol->rhoY[n] = rhoY;
				// Sol->geopot[n] = g * y;

#ifdef THERMCON
                mpv->p2_cells[n] = (p/rhoY) / ud.Msq;
#else
                mpv->p2_cells[n] = (p - mpv->HydroState->p0[j]) / ud.Msq;
#endif
			}
		}        
	}  
	
    /*set nodal pressures to hydrostatic values */
    for(k = 0; k < iczn; k++) {l = k * icxn * icyn;   
        
        for(j = 0; j < icyn; j++) {m = l + j * icxn;                
            double p    = mpv->HydroState_n->p0[j];
            double rhoY = mpv->HydroState_n->rhoY0[j];
            
            for(i = 0; i < icxn; i++) {n = m + i;
                mpv->p2_nodes[n] = ((p-mpv->HydroState_n->p0[j])/rhoY) / ud.Msq;
            }
        }
    }   
    
    ud.nonhydrostasy   = nonhydrostasy(0);
    ud.compressibility = compressibility(0);
    
    set_wall_rhoYflux(bdry, Sol, mpv, elem);
    Set_Explicit_Boundary_Data(Sol, elem);
    
    ConsVars_set(Sol0, Sol, elem->nc);
    
    /* the initial projection should ensure the velocity field is discretely
     divergence-controlled when sound-free initial data are required.
     This can mean vanishing divergence in a zero-Mach flow or the 
     pseudo-incompressible divergence constraint in an atmospheric flow setting
     */ 
    if (ud.initial_projection == CORRECT) {
        /* initial velocity div-control not implemented for the test case yet. */
        assert(0);
    }

}

/* ================================================================================== */

double pressure_function(double r, double p0, double S0, double u_theta, double Msq, double Gamma){
	return pow((pow(p0,Gamma) + Gamma*S0*Msq*u_theta*u_theta*(1.0 - pow((1.0-r),5.0)*(5.0*r+1.0))/30.0), (1.0/Gamma));
}

/* ================================================================================== */

double rho_function(double psi){
	
	double rhomax = 1.0;
	double percent = 0.0;
	double delta   = 0.1;
	
	double smooth, cosine, sgncos;
	
	cosine = cos(PI*(fabs(psi) - (0.5 - 0.5*delta)) / delta);
	sgncos = SIGNnull(cosine);
	
	smooth = 0.5*(1.0 - sgncos*sqrt(sqrt(sqrt(sgncos*cosine))) );
	
	/* return(1.0 - 0.2 * (psi*psi)); */
	/* return(rhomax * (1.0 - percent * (fabs(psi) < 0.5 ? 0.0 : 1.0 )) ); */
	return(rhomax * (1.0 - percent * (fabs(psi) < 0.5 - 0.5*delta ?  
									  0.0 : 
									  (fabs(psi) > 0.5 + 0.5*delta ? 
                                       1.0 : 
                                       smooth
									   ) 
									  )
					 ) 
		   ); 
}

/* ================================================================================== */

double stratification(
					  double y) {
	
	extern User_Data ud;	
	return( 1.0 );
}



/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: userdata.c,v $
 Revision 1.2  1998/03/07 09:56:47  nicola
 Added flux computation and multiple pressure variables.
 
 Revision 1.1  1998/03/01 18:50:44  nicola
 Added csr.c and sod1d.c (user data for 3d code)
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
