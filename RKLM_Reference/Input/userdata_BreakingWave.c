/*******************************************************************************
 File:   userdata.c
 Author: Nicola
 Date:   Fri Feb 27 09:34:03 CET 1998
 *******************************************************************************/
#include <math.h>
#include <float.h>
#include <stdlib.h>

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
#include "thermodynamic.h"
#include "Eos.h"
#include "set_ghostcells_p.h"
#include "set_ghostnodes_p.h"
#include "boundary.h"
#include "memory.h"
#include "Hydrostatics.h"

double pressure_function(double r, double p0, double S0, double u_theta, double Msq, double Gamma);

double rho_function(double psi);


void User_Data_init(User_Data* ud) {
	
	int i, max_no_of_levels;
	
	/* ================================================================================== */
	/* =====  PROBLEM SET UP  =========================================================== */
	/* ================================================================================== */
	
    /* Earth */	
	double grav  = 10.0;             /* [m/s^2]                         */
    double omega = 0.0; /* 2*PI*sin(0.25*PI)/(24.0*3600.0); [s^-1] */

	/* references for non-dimensionalization */
	double h_ref = 6515;             /* [m]                             */
	double t_ref = 100;              /* [s]                             */
	double T_ref = 227;              /* [K]                             */
	double p_ref = 101325;           /* [Pa]                            */
	double u_ref = h_ref/t_ref;      /* Strouhal No == 1 always assumed */
		
	/* thermodynamics and chemistry */
	double R_gas = 287;              /* [J/kg/K]                        */
    double R_vap = 461.00;           /* [J/kg/K]                        */
    double Q_vap = 2.53e+06;         /* [J]                             */
	double Nsq_ref   = 1e-5 * grav;      /* [s^-2]                          */
	double c     = grav*grav / (R_gas*T_ref*Nsq_ref);
	double gamma = c / (c-1);        /* breaking wave test choice yields exp. pressure fct. */
	ud->gamm     = gamma;  
	
    ud->h_ref       = h_ref;
    ud->t_ref       = t_ref;
    ud->T_ref       = T_ref;
    ud->p_ref       = p_ref;
    ud->u_ref       = u_ref;
    ud->Nsq_ref     = Nsq_ref;
    ud->g_ref       = grav;
    ud->gamm        = gamma;
    ud->Rg_over_Rv  = R_gas/R_vap;
    ud->Q           = Q_vap/(R_gas*T_ref);

    /* number of advected species */
    ud->nspec       = NSPEC;

	/* Low Mach */
    ud->is_compressible   = 0;   /* 0:psinc; 1:comp;  -1:psinc-comp-trans -> compressibility() */
    ud->acoustic_timestep = 0; /* 0;  1; */
    ud->Msq =  u_ref*u_ref / (R_gas*T_ref);
    
    /* geo-stuff */
    for(i=0; i<3; i++) {
        ud->gravity_strength[i]  = 0.0;   /* corresponds to  M^2 / Fr^2  =  g href / R Tref */
        ud->coriolis_strength[i] = 0.0;
    }
    ud->gravity_strength[1] = grav * h_ref / (R_gas * T_ref);
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
	ud->xmin = - 60000/h_ref;  
	ud->xmax =   60000/h_ref;  
	ud->ymin =   0.0;
	ud->ymax =   60000/h_ref; 
	ud->zmin = - 1.0;
	ud->zmax =   1.0;

	/* boundary/initial conditions */
	ud->wind_speed        = 10/u_ref;        /* velocity in [m/s] */             
    ud->wind_shear        = -0.0;            /* velocity in [u_ref/h_ref] */
	ud->hill_height       = 628.319/h_ref;   /* height   in [m]   */ 
	ud->hill_length_scale = 1000.0/h_ref;    /* length   in [m]   */   
	
	ud->bdrytype_min[0] = PERIODIC; /* DIRICHLET; */
	ud->bdrytype_min[1] = WALL; /* SLANTED_WALL; */
	ud->bdrytype_min[2] = WALL;
	ud->bdrytype_max[0] = PERIODIC; /* DIRICHLET; */  
	ud->bdrytype_max[1] = WALL;  
	ud->bdrytype_max[2] = WALL;
	
    ud->absorber = WRONG; /* CORRECT; */ 
	
	/* ================================================================================== */
	/* =====  NUMERICS  ================================================================= */
	/* ================================================================================== */
	
    /* Time discretization */
    ud->time_integrator        = OP_SPLIT; /* OP_SPLIT, OP_SPLIT_MD_UPDATE, HEUN, EXPL_MIDPT */
    ud->CFL                    = 0.96; /* 0.9; 0.8; */
    ud->dtfixed0               = 0.31;
    ud->dtfixed                = 0.31;  
    ud->stepmax                = 10000;
    ud->no_of_steps_to_CFL     = 1;
    ud->no_of_steps_to_dtfixed = 1;
    
    set_time_integrator_parameters(ud);
    
	/* Grid and space discretization */
	ud->inx = 240+1; /* 641; 321; 161; 129; 81; */    
	ud->iny = 120+1; /* 321; 161;  81;  65; 41;  */
	ud->inz = 1;
	ud->h    = MIN_own((ud->xmax-ud->xmin)/(ud->inx),MIN_own((ud->ymax-ud->ymin)/(ud->iny),(ud->zmax-ud->zmin)/(ud->inz)));

	/* explicit predictor step */
	/* Recovery */
	ud->recovery_order = SECOND;
    ud->limiter_type_scalars  = RUPE; 
    ud->limiter_type_velocity = RUPE; 
    /* Limiter options:  RUPE; NONE; MONOTONIZED_CENTRAL; MINMOD; VANLEER; SWEBY_MUNZ; SUPERBEE; */
	
    ud->kp = 1.4;
	ud->kz = 1.4; /* Entro abused for velocity in split-step-aligned velocity ! */
	ud->km = 1.4;
	ud->kY = 1.4;
	ud->kZ = 1.4; /* 2.0 */
	
    /* first correction */
    ud->p_flux_correction = CORRECT; /* CORRECT, WRONG; */
    if (ud->time_integrator == OP_SPLIT || ud->time_integrator == OP_SPLIT_MD_UPDATE) {
        ud->latw[0] = ud->latw[2] = 0.125; ud->latw[1] = 0.75; ud->p_extrapol = 1.0; /* REFERENCE */
        /* ud->latw[0] = ud->latw[2] = 0.0; ud->latw[1] = 1.0; ud->p_extrapol = 1.0; */   
        /* ud->latw[0] = ud->latw[2] = 0.25; ud->latw[1] = 0.5; ud->p_extrapol = 1.0; */  
        /* ud->latw[0] = ud->latw[2] = 0.125; ud->latw[1] = 0.75; ud->p_extrapol = 1.25; */
        /* ud->latw[0] = ud->latw[2] = 0.2; ud->latw[1] = 0.6; ud->p_extrapol = 1.5;  */
    } else {
        /* ud->latw[0] = ud->latw[2] = 0.125; ud->latw[1] = 0.75; ud->p_extrapol = 1.25;*/
        ud->latw[0] = ud->latw[2] = 0.25; ud->latw[1] = 0.5; ud->p_extrapol = 1.5;
    }

    ud->ncache =  300; /* (ud->inx+3); */
	/* ud->ncache = (ud->inx+3)*(ud->iny == 1 ? 1 : ud->iny+3)*(ud->inz == 1 ? 1 : ud->inz+3); */
	/* 5000+4 for 64x64; 10000 + 4; for 60x120 -- 31000+4 for 120x240*//*  BREAKING WAVE CHANGE */
	
    /* linear solver-stuff */
    ud->which_projection_first = 1;
    ud->Solver = BICGSTAB_PRECON;        /* options:   JACOBI, BICGSTAB, BICGSTAB_PRECON */
    ud->Solver_Node = BICGSTAB_PRECON;   /* options:   JACOBI, BICGSTAB, BICGSTAB_PRECON */
    ud->precondition = WRONG;            /* options:   CORRECT, WRONG */
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

    /* output times  */
    ud->tout[0] =  9000.0/t_ref;             /* times in [s]    */
    ud->tout[1] =  9900.0/t_ref;
    ud->tout[2] = 10800.0/t_ref;
    ud->tout[3] = 11700.0/t_ref;
    ud->tout[4] = 12600.0/t_ref;
    ud->tout[5] = -1.0;

    ud->stepmax = 10000;

    ud->write_stdout = ON;
    ud->write_stdout_period = 1;
    ud->write_file = ON;
    ud->write_file_period = 1;
    ud->file_format = HDF;
    
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
    const double delth = 0.0;
	
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
	double x, y, z, y_p, y_m;
	double rho, u, v, w, p, rhoY, p0;
	
    double S_integral_m, S_integral_n, S_integral_p, p2_o, g;
	double S_m, S_p, p_hydro, rhoY_hydro, rhoY_hydro_n;
	        
	g = ud.gravity_strength[1];
	
    Hydrostatics_State(mpv, elem);

	for(k = igz; k < icz - igz; k++) {l = k * icx * icy; 
		z = elem->z[k];
        
        for(j = igy; j < icy - igy; j++) {m = l + j * icx;
			            
			/*
            mpv->HydroState->geopot[j]
            mpv->HydroState->rho0[j]  
            mpv->HydroState->p0[j]    
            mpv->HydroState->S0[j]    
            mpv->HydroState->S10[j]   
            mpv->HydroState->Y0[j]    
            mpv->HydroState->rhoY0[j] 
            */
            
            y = elem->y[j];

			for(i = 0; i < icx; i++) {n = m + i;
                double r;
                
				x = elem->x[i];
				
                u     = u0;
                v     = v0;
                w     = w0;
				r     = sqrt(x*x + (y-0.2)*(y-0.2)) / 0.2;
                p     = mpv->HydroState->p0[j];
                rhoY  = mpv->HydroState->rhoY0[j];
				rho   = mpv->HydroState->rhoY0[j] / (stratification(y) + (r >= 1.0 ? 0.0 : (delth/300.0)*cos(0.5*PI*r)*cos(0.5*PI*r)));
				
				Sol->rho[n]  = rho;
				Sol->rhou[n] = rho * u;
				Sol->rhov[n] = rho * v;
				Sol->rhow[n] = rho * w;
				Sol->rhoe[n] = rhoe(rho, u, v, w, p, g*y);
				Sol->rhoY[n] = rhoY;
				Sol->geopot[n] = g*y;

                mpv->p2_cells[n] = p / ud.Msq;
				Sol->rhoZ[n]     = rho * mpv->p2_cells[n];				
			}			
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
	
	const double delS_slope = 0.0651466;  /* 0.1 ClarkeFarley; 0.0651466 Breaking Wave */ /*  BREAKING WAVE CHANGE */
	
	/* Clarke-Farley-Bubble test
	 return(4.0 * ( (exp(-delS_slope * y) - 1.0) / (1.0 - exp(-delS_slope * ud.ymax)) )); 
	 */
	
	/* breaking waves */
	 return( exp(delS_slope * y) );
	 
	/* rising bubble 
	return( 1.0 );
	 */
}



/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: userdata.c,v $
 Revision 1.2  1998/03/07 09:56:47  nicola
 Added flux computation and multiple pressure variables.
 
 Revision 1.1  1998/03/01 18:50:44  nicola
 Added csr.c and sod1d.c (user data for 3d code)
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
