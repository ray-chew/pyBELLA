/*******************************************************************************
 File:   userdata.c
 Author: Nicola
 Date:   Fri Feb 27 09:34:03 CET 1998
 *******************************************************************************/
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include "math_own.h"
#include "enumerator.h"
#include "Common.h"
#include "userdata.h"
#include "time_discretization.h"
#include "error.h"
#include "variable.h"
#include "enum_bdry.h"
#include "thermodynamic.h"
#include "Eos.h"
#include "set_ghostcells_p.h"
#include "set_ghostnodes_p.h"
#include "boundary.h"
#include "memory.h"
#include "Hydrostatics.h"

/*
 Monitor configurations:
 
 OpSplit:
 compressible case 
 unless otherwise listed, I tested options at fixed: 
 CFL           0.96, 
 grid          100x100, 
 t_end         4.0; 
 tol           1.0e-09, 
 limiters      NONE, 
 ud.p_extrapol 1.0
 ccenterweight 2.0
 HYDROSTATES_USING_dp2 on (switched off later)
 Diagonal Five Point in second projection
 ud->p_flux_correction = CORRECT
 ud->p_average = WRONG
 DP_AVERAGE off in ProjectionType.h
 DPDT_AVERAGED_IN_HELMHOLTZ off in ProjectionType.h
 
 Results:
 -- CFL 0.96                         -> very distorted
 -- HYDROSTATES_USING_dp2 off        -> that is not it ... pressure fluctuations stay put
 -- HYDROSTATES_USING_dp2 on         -> that is not it ... pressure fluctuations stay put same thing
 -- ud->p_flux_correction = WRONG    -> that is not it ... pressure fluctuations stay put same thing
 
 Found a bug in the initial data (power of p' when computing rhoY)
 -- this was a huge improvement, yet not the entire story.  
    ->  try running compressibility from the start !
 -- that did the trick. 
 
 -- 400x400 works, too.
 
 So this version of the code runs: 
   * the rising bubble (albeit with HYDROSTATES_USING_dp2 off )
   * the travelling vortex
 
 PREV:   RISING BUBBLE
 NEXT:   STRAKA

 */ 

double pressure_function(double r, double p0, double S0, double u_theta, double Msq, double Gamma);

double rho_function(double psi);


void User_Data_init(User_Data* ud) {
	
	int i, max_no_of_levels;
	
	/* ================================================================================== */
	/* =====  PROBLEM SET UP  =========================================================== */
	/* ================================================================================== */
	
    /* Earth */	
	double grav  = 0.0; /* [m/s^2]                                 */
	double omega = 0.0; /* 2*PI*sin(0.25*PI)/(24.0*3600.0); [s^-1] */

	/* references for non-dimensionalization */
	double h_ref = 10;               /* [m]                             */
	double t_ref = 10;               /* [s]                             */
	double T_ref = 300;              /* [K]                             */
    double p_ref = 1.0e+05;          /* [Pa]                            */
	double u_ref = h_ref/t_ref;      /* Strouhal No == 1 always assumed */
		
	/* thermodynamics and chemistry */
	double R_gas = 287.04;           /* [J/kg/K]                        */
	double R_vap = 461.00;           /* [J/kg/K]                        */
    double Q_vap = 2.53e+06;         /* [J]                             */
	double gamma = 1.4;              /* special choice for breaking wave test to get exp pressure fct. */
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
    ud->naux        = NAUX;  

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
    ud->implicit_gravity_theta2 = 0;

	/* flow domain */
	ud->xmin = - 0.5;  
	ud->xmax =   0.5;  
	ud->ymin = - 0.5;
	ud->ymax =   0.5; 
	ud->zmin =   0.0;
	ud->zmax =   0.1;

	/* boundary/initial conditions */
	ud->wind_speed        =  1.0;              /* velocity in [u_ref] */
	ud->wind_shear        = -0.0;              /* velocity in [u_ref/h_ref] */             
	ud->hill_height       =  0.0;              /* height   in [h_ref]   */ 
	ud->hill_length_scale =  99999.9;          /* width    in [h_ref]   */   
	
	ud->bdrytype_min[0] = PERIODIC; /* DIRICHLET; */
	ud->bdrytype_min[1] = PERIODIC; /* SLANTED_WALL; */
	ud->bdrytype_min[2] = WALL;
	ud->bdrytype_max[0] = PERIODIC; /* DIRICHLET; */  
	ud->bdrytype_max[1] = PERIODIC;  
	ud->bdrytype_max[2] = WALL;
	
	ud->absorber = WRONG; /* CORRECT; */ 
	
	/* ======================================================================== */
	/* =====  NUMERICS  ======================================================= */
	/* ======================================================================== */
	
    /* time discretization */
    ud->time_integrator      = OP_SPLIT;  /* OP_SPLIT, OP_SPLIT_MD_UPDATE, HEUN, EXPL_MIDPT; RK3_SKAMA; RK3_TEST */
	ud->CFL                  = 0.96;       /* 0.9; 0.8; 1.1*sqrt(0.5) */
    ud->dtfixed0             = 10000.999;
    ud->dtfixed              = 10000.999;     /* ud->dtfixed = 0.5/ud->t_ref;  */
    ud->no_of_steps_to_CFL   = 1;
    ud->no_of_steps_to_dtfixed = 1;

    set_time_integrator_parameters(ud);
    
	/* Grid and space discretization */
	ud->inx =  128+1; /*  */
	ud->iny =  128+1; /*  */
	ud->inz =  1;
	ud->h   = MIN_own((ud->xmax-ud->xmin)/(ud->inx),MIN_own((ud->ymax-ud->ymin)/(ud->iny),(ud->zmax-ud->zmin)/(ud->inz)));

	/* explicit predictor step */
	/* Recovery */
	ud->recovery_order        = SECOND;
	ud->limiter_type_scalars  = VANLEER; /*  RUPE; NONE; MONOTONIZED_CENTRAL; MINMOD; VANLEER; VANLEERSmooth; SWEBY_MUNZ; SUPERBEE; NO_SLOPE;*/
	ud->limiter_type_velocity = VANLEER; /*  RUPE; NONE; MONOTONIZED_CENTRAL; MINMOD; VANLEER; VANLEERSmooth; SWEBY_MUNZ; SUPERBEE; NO_SLOPE;*/

    /* first correction */
	ud->p_flux_correction = WRONG; /* CORRECT, WRONG; */
    if (ud->time_integrator == OP_SPLIT || ud->time_integrator == OP_SPLIT_MD_UPDATE) {
        ud->latw[0] = ud->latw[2] = 0.125; ud->latw[1] = 0.75; ud->p_extrapol = 1.0;
        /* ud->latw[0] = ud->latw[2] = 0.125; ud->latw[1] = 0.75; ud->p_extrapol = 1.25; */
        /* ud->latw[0] = ud->latw[2] = 0.25; ud->latw[1] = 0.5; ud->p_extrapol = 1.0; */
        /* ud->latw[0] = ud->latw[2] = 0.2; ud->latw[1] = 0.6; ud->p_extrapol = 1.5;  */ 
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
	
	ud->ncache =  201; /* (ud->inx+3); */
	
	/* linear solver-stuff */
    ud->which_projection_first = 1;
	ud->Solver = BICGSTAB_PRECON;        /* options:   JACOBI, BICGSTAB, BICGSTAB_PRECON */
	ud->Solver_Node = BICGSTAB_PRECON;   /* options:   JACOBI, BICGSTAB, BICGSTAB_PRECON */
    ud->precondition = WRONG;
    double tol = 1.e-6;
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
		
	/* ========================================================================== */
	/* =====  FLOW CONTROL  ===================================================== */
	/* ========================================================================== */
    
    /* output times and max no of time steps 
    ud->tout[0]  = 0.2;      
	ud->tout[1]  = 0.4;      
	ud->tout[2]  = 0.6;
	ud->tout[3]  = 0.8;
	ud->tout[4]  = 1.0;
	ud->tout[5]  = 1.2;
	ud->tout[6]  = 1.4;
	ud->tout[7]  = 1.6;
	ud->tout[8]  = 1.8;
	ud->tout[9]  = 2.0;
	ud->tout[10] = -1.0;
     */
    /* output times and max no of time steps */
    ud->tout[0]  = 0.4;      
    ud->tout[1]  = 0.8;      
    ud->tout[2]  = 1.2;
    ud->tout[3]  = 1.6;
    ud->tout[4]  = 2.0;
    ud->tout[5]  = 2.4;
    ud->tout[6]  = 2.8;
    ud->tout[7]  = 3.2;
    ud->tout[8]  = 3.6;
    ud->tout[9]  = 4.0;
    ud->tout[10] = -1.0;

    
    ud->stepmax = 10000;

	ud->write_stdout = ON;
	ud->write_stdout_period = 1;
	ud->write_file = ON;
	ud->write_file_period = 10000;
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

void Sol_initial(ConsVars* Sol,
                 const ElemSpaceDiscr* elem,
                 const NodeSpaceDiscr* node) {
	
	extern Thermodynamic th;
	extern User_Data ud;
    extern MPV* mpv;
    extern double *Yinvbg;
    
	const double u0    = 1.0*ud.wind_speed;
	const double v0    = 1.0*ud.wind_speed;
	const double w0    = 0.0;
    
    const double rotdir = -1.0;
    
	const double rho0    = 0.5;
    const double del_rho = 0.5;  /* 0.0; for homentropic */
    const double R0      = 0.4;
    const double fac     = 1*1024.0; /* 4*1024.0 */
			
	const int icx  = elem->icx;
	const int icy  = elem->icy;
	const int icz  = elem->icz;
	const int igy  = elem->igy;
	const int igz  = elem->igz;

    const int icxn = node->icx;
    const int icyn = node->icy;
    const int iczn = node->icz;
	    
	int i, j, k, l, m, n;
	double x, y, z;
	double rho, u, v, w, rhoY, qv, qc, qr, theta, T, p_hydro;
    double qvs;
	
    double g;
	                
    g = ud.gravity_strength[1];

    Hydrostatics_State(mpv, Yinvbg, elem);
    
    /* Initial data and hydro-states in the flow domain */
	for(k = igz; k < icz - igz; k++) {l = k * icx * icy; 
		z = elem->z[k];
        
        for(j = igy; j < icy - igy; j++) {m = l + j * icx;
			
            y = elem->y[j];
			
            for(i = 0; i < icx; i++) {n = m + i;
                double r, uth;
                
                x       = elem->x[i];
                r       = sqrt(x*x + y*y);
                uth     = rotdir * (r < R0 ? fac * pow( 1.0-r/R0, 6) * pow( r/R0, 6) : 0.0);
                
                u       = u0 + uth * (-y/r);
                v       = v0 + uth * (+x/r);
                w       = w0;
                p_hydro = mpv->HydroState->p0[j];
                rhoY    = mpv->HydroState->rhoY0[j];
                qc      = 0.0;
                qr      = 0.0;
                theta   = stratification(y);
                rho     =  (r < R0 ? (rho0 + del_rho*pow( 1-(r/R0)*(r/R0) , 6)) : rho0);
                T       = T_from_p_rho(p_hydro,rho);
                qvs     = 1.0;
                qv      = 0.0;
                
                Sol->rho[n]  = rho;
                Sol->rhou[n] = rho * u;
                Sol->rhov[n] = rho * v;
                Sol->rhow[n] = rho * w;
                Sol->rhoe[n] = rhoe(rho, u, v, w, p_hydro, g*y);
                Sol->rhoY[n] = rhoY;
                
                if (ud.nspec == 3) {
                    Sol->rhoX[QV][n] = rho*qv;
                    Sol->rhoX[QC][n] = rho*qc;
                    Sol->rhoX[QR][n] = rho*qr;                    
                }
                
                Sol->geopot[n] = g * y;
                
                if ( r/R0 < 1.0 ) {
                    
                    int ii;
                    double coe[25];
                                        
                    coe[0]  =     1.0 / 12.0;
                    coe[1]  = -  12.0 / 13.0;
                    coe[2]  =     9.0 /  2.0;
                    coe[3]  = - 184.0 / 15.0;
                    coe[4]  =   609.0 / 32.0;
                    coe[5]  = - 222.0 / 17.0;
                    coe[6]  = -  38.0 /  9.0; 
                    coe[7]  =    54.0 / 19.0;
                    coe[8]  =   783.0 / 20.0;
                    coe[9]  = - 558.0 /  7.0;
                    coe[10] =  1053.0 / 22.0;
                    coe[11] =  1014.0 / 23.0;
                    coe[12] = -1473.0 / 16.0;
                    coe[13] =   204.0 /  5.0;
                    coe[14] =   510.0 / 13.0;
                    coe[15] = -1564.0 / 27.0;
                    coe[16] =   153.0 /  8.0;
                    coe[17] =   450.0 / 29.0;
                    coe[18] = - 269.0 / 15.0; /* Kadioglu et al.: 259; Papke: 269 */
                    coe[19] =   174.0 / 31.0;
                    coe[20] =    57.0 / 32.0;
                    coe[21] = -  74.0 / 33.0;
                    coe[22] =    15.0 / 17.0;
                    coe[23] = -   6.0 / 35.0;
                    coe[24] =     1.0 / 72.0;
                    
                    mpv->p2_cells[n] = 0.0;
                    
                    for (ii = 0; ii < 25; ii++)
                    {
                        mpv->p2_cells[n] += coe[ii] * (pow(r/R0 ,12+ii) - 1.0) * rotdir * rotdir;
                    }
                }
                else {
                    mpv->p2_cells[n] = 0.0;
                }
                
                /* regular thermodynamic pressure */
                /* mpv->p2_cells[n] = p_hydro/ud.Msq + fac*fac*mpv->p2_cells[n]; */
#ifdef THERMCON
                /* pseudo-incompressible option */
                mpv->p2_cells[n]  = th.Gamma*fac*fac*mpv->p2_cells[n];
                mpv->p2_cells[n]   = ud.is_compressible*pow(p_hydro+mpv->p2_cells[n]*ud.Msq, th.Gamma)/ud.Msq \
                                    + (1-ud.is_compressible)*mpv->p2_cells[n]/mpv->HydroState->rhoY0[j];
                Sol->rhoZ[PRES][n] = mpv->p2_cells[n];
                Sol->rhoY[n]       = ud.is_compressible*pow(p_hydro+mpv->p2_cells[n]*ud.Msq,1.0/th.gamm) \
                                    + (1-ud.is_compressible)*Sol->rhoY[n];
#else
                mpv->p2_cells[n] = fac*fac*mpv->p2_cells[n];
                assert(0);
#endif
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

    /*set nodal pressures */
    for(int k = 0; k < iczn; k++) {int l = k * icxn * icyn;   
        
        for(int j = 0; j < icyn; j++) {int m = l + j * icxn;                
#ifdef THERMCON
            double p    = mpv->HydroState->p0[j];
            double rhoY = mpv->HydroState->rhoY0[j];
            
            for(int i = 0; i < icxn; i++) {int n = m + i;
                mpv->p2_nodes[n] = (p/rhoY) / ud.Msq;
            }
#else
            double p    = mpv->HydroState->p0[j];
            for(int i = 0; i < icxn; i++) {int n = m + i;
                mpv->p2_nodes[n] = p / ud.Msq;
            }
#endif
        }
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

double stratification(
					  double y) {
	
	extern User_Data ud;
			
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
