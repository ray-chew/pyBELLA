/*******************************************************************************
 File:   userdata.c
 Author: Nicola
 Date:   Fri Feb 27 09:34:03 CET 1998
 *******************************************************************************/
#include "enumerator.h"
#include "Common.h"
#include "ProjectionType.h"
#include "userdata.h"
#include "time_discretization.h"
#include "error.h"
#include "variable.h"
#include "enum_bdry.h"
#include "math_own.h"
#include "space_discretization.h"
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include "thermodynamic.h"
#include "Eos.h"
#include "set_ghostcells_p.h"
#include "set_ghostnodes_p.h"
#include "boundary.h"
#include "memory.h"

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


static void (*rotate[])(ConsVars* Sol, double* rhs, double *Yinvbg, const enum Direction dir) = {NULL, rotate2D, rotate3D};

double pressure_function(double r, double p0, double S0, double u_theta, double Msq, double Gamma);

double rho_function(double psi);


void User_Data_init(User_Data* ud) {
	
	int i, max_no_of_levels;
	
	/* ================================================================================== */
	/* =====  PROBLEM SET UP  =========================================================== */
	/* ================================================================================== */
	
    /* Earth */	
	double grav  = 9.81;                               /* [m/s^2]                         */
	double omega = 2*PI*sin(0.25*PI)/(24.0*3600.0);    /*  [s^-1]                         */

	/* references for non-dimensionalization */
	double h_ref = 10000;            /* [m]                             */
	double t_ref = 100;              /* [s]                             */
	double T_ref = 300;              /* [K]                             */
	double p_ref = 101325;           /* [Pa]                            */
	double u_ref = h_ref/t_ref;      /* Strouhal No == 1 always assumed */
		
	/* thermodynamics and chemistry */
	double R_gas = 287.4;            /* [J/kg/K]                        */
	double R_vap = 461.00;           /* [J/kg/K]                        */
    double Q_vap = 2.53e+06;         /* [J]                             */
	double gamma = 1.4;              /* dimensionless                   */
	ud->gamm = gamma;  
    double Nsq   = grav*1.3e-05;     /* [] */

    ud->h_ref       = h_ref;
    ud->t_ref       = t_ref;
    ud->T_ref       = T_ref;
    ud->p_ref       = p_ref;
    ud->u_ref       = u_ref;
    ud->Nsq_ref     = Nsq;
    ud->g_ref       = grav;
    ud->gamm        = gamma;  
    ud->Rg_over_Rv  = R_gas/R_vap;  
    ud->Q           = Q_vap/(R_gas*T_ref);  
    
    /* number of advected species */
    ud->nspec       = NSPEC;  

	/* Low Mach */
    ud->is_compressible = 0;
    ud->acoustic_timestep =  0; /* 0;  1; */
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
    ud->implicit_gravity_theta  = 0;
    ud->implicit_gravity_press  = 0; /* should be on for compressible option */
    ud->implicit_gravity_theta2 = 0;
    ud->implicit_gravity_press2 = 0; /* should this, too, be on for compressible option ?  */
    ud->thermcon = 0;
	    
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
	ud->hill_height = 0.0 * 0.096447; /* 0.096447; */ /* hill_height * l_ref = 0.628319 km *//*  BREAKING WAVE CHANGE */
	ud->hill_length_scale = 0.1535;   /* hill_length * l_ref = 1.0 km */
	
	ud->bdrytype_min[0] = WALL; /* DIRICHLET; PERIODIC; WALL; */
	ud->bdrytype_min[1] = WALL; /* SLANTED_WALL; */
	ud->bdrytype_min[2] = WALL;
	ud->bdrytype_max[0] = WALL; /* DIRICHLET; PERIODIC; WALL; */  
	ud->bdrytype_max[1] = WALL;  
	ud->bdrytype_max[2] = WALL;
	
	ud->absorber = WRONG; /* CORRECT;  WRONG; */ /*  BREAKING WAVE CHANGE */
	
	
	/* ================================================================================== */
	/* =====  NUMERICS  ================================================================= */
	/* ================================================================================== */
	
    /* time discretization */
    ud->time_integrator      = OP_SPLIT;  /* OP_SPLIT, HEUN, EXPL_MIDPT */
	ud->CFL                  = 0.96; /* 0.45; 0.9; 0.8; */
    ud->dtfixed0             = 0.05;
    ud->dtfixed              = 0.05; /* 0.0052; */ /*  0.004; */ 
    ud->no_of_steps_to_CFL   = 1;
    ud->no_of_steps_to_dtfixed = 1;
    
    set_time_integrator_parameters(ud);
    
	/* Grid and space discretization */
	ud->inx = 160+1; /* 641; 321; 161; 129; 81; */    
	ud->iny =  80+1; /* 321; 161;  81;  65; 41;  */
	ud->inz =  1;
	ud->h    = MIN_own((ud->xmax-ud->xmin)/(ud->inx),MIN_own((ud->ymax-ud->ymin)/(ud->iny),(ud->zmax-ud->zmin)/(ud->inz)));
	
	/* explicit predictor step */
	/* Recovery */
	ud->recovery_order = SECOND;
	ud->limiter_type_scalars  = NONE; /*  RUPE; NONE; MONOTONIZED_CENTRAL; MINMOD; VANLEER; SWEBY_MUNZ; SUPERBEE; */
	ud->limiter_type_velocity = NONE; /*  RUPE; NONE; MONOTONIZED_CENTRAL; MINMOD; VANLEER; SWEBY_MUNZ; SUPERBEE; */
	
    /* first correction */
	ud->p_flux_correction = CORRECT; /* CORRECT, WRONG; */
    if (ud->time_integrator == OP_SPLIT || ud->time_integrator == OP_SPLIT_MD_UPDATE) {
        ud->latw[0] = ud->latw[2] = 0.125; ud->latw[1] = 0.75; ud->p_extrapol = 1.0; 
        /* ud->latw[0] = ud->latw[2] = 0.125; ud->latw[1] = 0.75; ud->p_extrapol = 1.25; */
        /* ud->latw[0] = ud->latw[2] = 0.25; ud->latw[1] = 0.5; ud->p_extrapol = 1.0; */
    } else {
        /* ud->latw[0] = ud->latw[2] = 0.125; ud->latw[1] = 0.75; ud->p_extrapol = 1.25;*/
        ud->latw[0] = ud->latw[2] = 0.25; ud->latw[1] = 0.5; ud->p_extrapol = 1.5; 
    }
    
    /* parameters for SWEBY_MUNZ limiter family */
    ud->kp = 1.4;
	ud->kz = 1.4; /* Entro abused for velocity in split-step-aligned velocity ! */
	ud->km = 1.4;
	ud->kY = 1.4;
	ud->kZ = 1.4; /* 2.0 */
	
	ud->ncache =  300; /* (ud->inx+3); */
	
	/* linear solver-stuff */
    ud->which_projection_first = 1;
	ud->Solver = BICGSTAB_PRECON;        /* options:   JACOBI, BICGSTAB, BICGSTAB_PRECON */
	ud->Solver_Node = BICGSTAB_PRECON;   /* options:   JACOBI, BICGSTAB, BICGSTAB_PRECON */
    double tol = 1.e-8;
	ud->flux_correction_precision = tol;
	ud->flux_correction_local_precision = tol;   /* 1.e-05 should be enough */
	ud->second_projection_precision = tol;  
	ud->second_projection_local_precision = tol;   /* 1.e-05 should be enough */
	ud->implicitness = 1.0;   
	ud->flux_correction_max_MG_cycles = 100;
	ud->flux_correction_output_period = 50;
	ud->max_projection_iterations = 1;
	ud->flux_correction_max_iterations = 6000;
	ud->second_projection_max_iterations = 6000;
	
    max_no_of_levels = 1;

	ud->max_no_of_multigrid_levels = max_no_of_levels; 
	ud->no_of_multigrid_levels     = max_no_of_levels-1;    /* optimal for BICGSTAB with MG:  5 (128x128-Grid) */

	/* numerics parameters */
	ud->eps_Machine = sqrt(DBL_EPSILON);
		
	/* ================================================================================== */
	/* =====  FLOW CONTROL  ============================================================= */
	/* ================================================================================== */

    /* output times */  
	ud->tout[0] =  2.0;      
	ud->tout[1] =  4.0;      
	ud->tout[2] =  6.0;      
	ud->tout[3] =  8.0;      
    ud->tout[4] = 10.0;      
    ud->tout[5] = 10.5;      
	ud->tout[6] = -1.0;      

    ud->stepmax = 10000;
    
	ud->write_stdout = ON;
	ud->write_stdout_period = 1;
	ud->write_file = ON;
	ud->write_file_period = 10000;
	ud->file_format = HDF;
    
    ud->write_history = 0;


    {
        char *OutputBaseFolder      = "/Users/rupert/Documents/Computation/RKLM_Reference/";
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

void Sol_initial(ConsVars* Sol, const ElemSpaceDiscr* elem, const NodeSpaceDiscr* node) {
	
	extern Thermodynamic th;
	extern User_Data ud;
    extern MPV* mpv;
    
	const double u0 = ud.wind_speed;
	const double v0 = 0.0;
	const double w0 = 0.0;
	const double rho0 = 1.0;
	const double delth = 2.0; /* pot temp perturbation in K */
	
    const double y0 = 0.2;
    const double r0 = 0.2;
    
	const double Gamma = th.gm1 / th.gamm;
	const double gamm  = th.gamm;
    const double Gamma_inv = 1.0/Gamma;
	
	const double M_LH_sq = ud.Msq;
	
	const int icx = elem->icx;
	const int icy = elem->icy;
	const int icz = elem->icz;
	const int igy = elem->igy;
	const int igz = elem->igz;
	
	int i, j, k, l, m, n;
	double x, y, z;
	double rho, u, v, w, p, rhoY, p0;
	
    double S_integral, p2_o, g;
	double elevation_m, elevation_p, S_m, S_p, p_hydro, rho_hydro;
	    
	g = ud.gravity_strength[1];
	
	/* Hydrostates in bottom dummy cells */
	y = elem->y[igy];
	elevation_m = y - (elem->y[elem->igy] - 0.5*elem->dy);
	S_m         = 1.0/stratification(elevation_m);
	S_integral  = 0.5 * elevation_m * (S_m + 1.0/stratification(0.0));
	p0 = pow(rho0 * stratification(0.0), gamm);
	
	for(j = igy; j >= 0; j--) {
		y = elem->y[j];
		elevation_m = y - (elem->y[elem->igy] - 0.5*elem->dy);
		
		S_m         = 1.0/stratification(elevation_m);
		p_hydro     = pow( pow(p0,Gamma) - Gamma*g*S_integral ,  Gamma_inv);
		rho_hydro   = S_m * pow(p_hydro,1.0/th.gamm);
		
		mpv->HydroState->geopot[j] = g * elevation_m;
		mpv->HydroState->rho0[j] = rho_hydro;
		mpv->HydroState->p0[j]   = p_hydro;
		mpv->HydroState->S0[j]   = S_m;
		mpv->HydroState->S10[j]  = 0.0;
		mpv->HydroState->Y0[j]  = stratification(elevation_m);
		
		elevation_p = elevation_m - elem->dy;
		S_p         = 1.0/stratification(elevation_p);
		S_integral -= 0.5*elem->dy*(S_m + S_p);
	}
	    
    /* data in the bulk of the domain */
	for(k = igz; k < icz - igz; k++) {l = k * icx * icy; 
		z = elem->z[k];
        
		y = elem->y[igy];
		elevation_p = y - (elem->y[elem->igy] - 0.5*elem->dy);
		S_p         = 1.0/stratification(elevation_p);
		S_integral  = 0.5 * elevation_p * (S_p + 1.0/stratification(0.0));
        
        for(j = igy; j < icy - igy; j++) {m = l + j * icx;
			
            elevation_m = elevation_p;
            S_m         = S_p;
            
            p_hydro     = pow( pow(p0,Gamma) - Gamma*g*S_integral ,  Gamma_inv);
            rho_hydro   = S_m * pow(p_hydro,1.0/th.gamm);
			
            if (k==igz) {
                mpv->HydroState->geopot[j] = g * elevation_m;
                mpv->HydroState->rho0[j] = rho_hydro;
                mpv->HydroState->p0[j]   = p_hydro;
                mpv->HydroState->S0[j]   = S_m;
                mpv->HydroState->S10[j]  = 0.0;
                mpv->HydroState->Y0[j]  = stratification(elevation_m);
            }

			for(i = 0; i < icx; i++) {n = m + i;
                double ystar, r;
                
				x = elem->x[i];
				
                u   = u0;
                v   = v0;
                w   = w0;
                ystar = elevation_m;

				r = sqrt(x*x + (ystar-y0)*(ystar-y0) + (elem->ndim==3)*z*z) / r0;
                p    = p_hydro;
				rhoY = pow(p_hydro,1.0/th.gamm);
				rho  = rhoY / (stratification(ystar) + (r >= 1.0 ? 0.0 : (delth/300.0)*cos(0.5*PI*r)*cos(0.5*PI*r)));
				
				Sol->rho[n]  = rho;
				Sol->rhou[n] = rho * u;
				Sol->rhov[n] = rho * v;
				Sol->rhow[n] = rho * w;
				Sol->rhoe[n] = rhoe(rho, u, v, w, p, g*elevation_m);
				Sol->rhoY[n] = rhoY;
				Sol->geopot[n] = g * elevation_m;

#define P_RHO
#ifdef P_RHO
				mpv->p2_cells[n] = pow(rhoY,th.gamm) / ud.Msq;
#else /* P_RHO */
				mpv->p2_cells[n] = th.Gammainv * pow(rhoY,th.gm1) / ud.Msq;
#endif /* P_RHO */
				Sol->rhoZ[n]     = rho * mpv->p2_cells[n];				
			}
            
            elevation_p = elevation_m + elem->dy;
            S_p         = 1.0/stratification(elevation_p);
            S_integral += 0.5*elem->dy*(S_m + S_p);
			
		}
        
        /* Hydrostates in top dummy cells */
        for(j = icy-igy; j < icy; j++) {
            y = elem->y[j];
            elevation_m = y - (elem->y[elem->igy] - 0.5*elem->dy);
            
            S_m         = 1.0 / stratification(elevation_m);
            p_hydro     = pow( pow(p0,Gamma) - Gamma*g*S_integral ,  Gamma_inv);
            rho_hydro   = S_m * pow(p_hydro,1.0/th.gamm);
            
            mpv->HydroState->geopot[j] = g * elevation_m;
            mpv->HydroState->rho0[j] = rho_hydro;
            mpv->HydroState->p0[j]   = p_hydro;
            mpv->HydroState->S0[j]   = S_m;
            mpv->HydroState->S10[j]  = 0.0;
            mpv->HydroState->Y0[j]  = stratification(elevation_m);
            
            elevation_p = elevation_m + elem->dy;
            S_p         = 1.0 / stratification(elevation_p);
            S_integral += 0.5*elem->dy*(S_m + S_p);
        }
        
		
		/* balanced "elliptic pressure" w.r.t. background state */
		p2_o = 0.0;
		
        for(j = igy; j < icy - igy; j++) {m = l + j * icx; 
			double p2_n;            
            mpv->HydroState->p20[j] = p2_o;
#define P_RHO
#ifdef P_RHO
            p2_n = (mpv->HydroState->p0[j+1] - mpv->HydroState->p0[igy]) / M_LH_sq;
#else /* P_RHO */
		    p2_n   = p2_o - (g/M_LH_sq) * 0.5 * elem->dy * 
                            (  (1.0/stratification(elem->y[j+1]) - ALPHA_BG) 
                             + (1.0/stratification(elem->y[j]) - ALPHA_BG) );
#endif
			
            p2_o = p2_n;
		}
		
		
		/* set all dummy cells */
		/* geopotential in bottom and top dummy cells */
		for(j = 0; j < igy; j++) {m = l + j * icx;  
			y = elem->y[j];
			for(i = 0; i < icx; i++) {n = m + i;
				Sol->geopot[n] = g * y;
			}
		}
		
		for(j = icy-igy; j < icy; j++) {m = l + j * icx;  
			y = elem->y[j];
			for(i = 0; i < icx; i++) {n = m + i;
				Sol->geopot[n] = g * y;
			}
		}
	}  
	
	
    /*            
    for(SplitStep = 0; SplitStep < ndim; SplitStep++) { 
        const double lambda = 1.0;
        Bound(Sol, mpv->HydroState,lambda, elem->nc, SplitStep); 
        if(SplitStep < ndim - 1) (*rotate[ndim - 1])(Sol, mpv->Level[0]->rhs, FORWARD);
    }         
    for(SplitStep = ndim-1; SplitStep > 0; SplitStep--) {
        (*rotate[ndim - 1])(Sol, mpv->Level[0]->rhs, BACKWARD);
    }	
     */
    
    /* put p2_cells into Z for advection */
    for(i=0; i<elem->nc; i++) {
        Sol->rhoZ[i]  = mpv->p2_cells[i] * Sol->rho[i];
    }
}


double pressure_function(double r, double p0, double S0, double u_theta, double Msq, double Gamma){
	return pow((pow(p0,Gamma) + Gamma*S0*Msq*u_theta*u_theta*(1.0 - pow((1.0-r),5.0)*(5.0*r+1.0))/30.0), (1.0/Gamma));
}

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

double stratification(
					  double y) {
	
	extern User_Data ud;
	
	/* const double delS_slope = 0.0651466; */  /* 0.1 ClarkeFarley; 0.0651466 Breaking Wave */ /*  BREAKING WAVE CHANGE */
	
	/* Clarke-Farley-Bubble test
	 return(4.0 * ( (exp(-delS_slope * y) - 1.0) / (1.0 - exp(-delS_slope * ud.ymax)) )); 
	 */
	
	/* breaking waves 
	 return( exp(delS_slope * y) );
	 */
	/* breaking waves */
	return( 1.0 );
}



/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: userdata.c,v $
 Revision 1.2  1998/03/07 09:56:47  nicola
 Added flux computation and multiple pressure variables.
 
 Revision 1.1  1998/03/01 18:50:44  nicola
 Added csr.c and sod1d.c (user data for 3d code)
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
