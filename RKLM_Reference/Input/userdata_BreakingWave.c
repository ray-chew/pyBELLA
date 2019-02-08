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
#include "second_projection.h"

double pressure_function(double r, double p0, double S0, double u_theta, double Msq, double Gamma);

double rho_function(double psi);


void User_Data_init(User_Data* ud) {
	
	int i;
	
	/* ================================================================================== */
	/* =====  PROBLEM SET UP  =========================================================== */
	/* ================================================================================== */
	
    /* Earth */	
	double grav  = 10.0;             /* [m/s^2]                         */
    double omega = 0.0; /* 2*PI*sin(0.25*PI)/(24.0*3600.0); [s^-1] */

    /* thermodynamics and chemistry */
    double R_gas = 287;              /* [J/kg/K]                        */
    double R_vap = 461.00;           /* [J/kg/K]                        */
    double Q_vap = 2.53e+06;         /* [J]                             */

    /* references for non-dimensionalization */
	double h_ref = 6515;             /* [m]                             */
	double t_ref = 100;              /* [s]                             */
	double T_ref = 227;              /* [K]                             */
	double p_ref = 101325;           /* [Pa]                            */
	double u_ref = h_ref/t_ref;      /* Strouhal No == 1 always assumed */
    double rho_ref = p_ref / (R_gas*T_ref); /* [kg/m^3]          */

    /* reference stratification as (buoyancy frequency)^2 */
    double Nsq_ref = 1e-5 * grav;      /* [s^-2]                          */
    double c     = grav*grav / (R_gas*T_ref*Nsq_ref);
    double gamma = c / (c-1);        /* breaking wave test choice yields exp. pressure fct. */
	
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

	/* Low Mach */
    ud->is_nonhydrostatic = 1; /* 0: hydrostatic;  1: nonhydrostatic;  -1: transition (see nonhydrostasy()) */
    ud->is_compressible   = 0; /* 0:psinc; 1:comp;  -1:psinc-comp-trans -> compressibility() */
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
	
    ud->absorber = CORRECT; /* CORRECT; */ 
	
	/* ================================================================================== */
	/* =====  NUMERICS  ================================================================= */
	/* ================================================================================== */
	
    /* Time discretization */
    ud->time_integrator        = SI_MIDPT;  /* this code version has only one option */
    ud->advec_time_integrator  = STRANG; /* HEUN; EXPL_MIDPT;   best tested: STRANG; */
    ud->CFL                    = 0.96; /* 0.9; 0.8; */
    ud->dtfixed0               = 0.31;
    ud->dtfixed                = 0.31;  
    
    set_time_integrator_parameters(ud);
    
	/* Grid and space discretization */
	ud->inx = 240+1; /* 641; 321; 161; 129; 81; */    
	ud->iny = 120+1; /* 321; 161;  81;  65; 41;  */
	ud->inz = 1;

	/* explicit predictor step */
	/* Recovery */
	ud->recovery_order = SECOND;
    ud->limiter_type_scalars  = NONE; 
    ud->limiter_type_velocity = NONE; 
    /* Limiter options:  RUPE; NONE; MONOTONIZED_CENTRAL; MINMOD; VANLEER; SWEBY_MUNZ; SUPERBEE; */
	
    ud->kp = 1.4;
	ud->kz = 1.4; /* Entro abused for velocity in split-step-aligned velocity ! */
	ud->km = 1.4;
	ud->kY = 1.4;
	ud->kZ = 1.4; /* 2.0 */	

    /* all explicit predictor operations are done on ncache-size data worms to save memory */ 
    ud->ncache =  300; /* (ud->inx+3); */
	
    /* linear solver-stuff */
    double tol = 1.e-10 * (ud->is_compressible == 1 ? 0.01 : 1.0);
    ud->flux_correction_precision         = tol;
    ud->flux_correction_local_precision   = tol;    /* 1.e-05 should be enough */
    ud->second_projection_precision       = tol;
    ud->second_projection_local_precision = tol;  /* 1.e-05 should be enough */
    ud->flux_correction_max_iterations    = 6000;
    ud->second_projection_max_iterations  = 6000;
    ud->initial_projection                = CORRECT;   /* WRONG;  CORRECT; */
    ud->initial_impl_Euler                = CORRECT;   /* WRONG;  CORRECT; */
    
    ud->column_preconditioner             = CORRECT; /* WRONG; CORRECT; */
    ud->synchronize_nodal_pressure        = WRONG; /* WRONG; CORRECT; */
    ud->synchronize_weight                = 1.0;    /* relevant only when prev. option is "CORRECT"
                                                     Should ultimately be a function of dt . */  

	/* numerics parameters */
	ud->eps_Machine = sqrt(DBL_EPSILON);
		
	/* ================================================================================== */
	/* =====  CODE FLOW CONTROL  ======================================================== */
	/* ================================================================================== */

    /* output times  */
    ud->tout[0] =  9000.0/t_ref;             /* times in [s]    */
    ud->tout[1] = -9900.0/t_ref;
    ud->tout[2] = 10800.0/t_ref;
    ud->tout[3] = 11700.0/t_ref;
    ud->tout[4] = 12600.0/t_ref;
    ud->tout[5] = -1.0;

    ud->stepmax = 10000;

    ud->write_stdout = ON;
    ud->write_stdout_period = 1;
    ud->write_file = ON;
    ud->write_file_period = 10;
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
                 ConsVars* Sol0, 
                 MPV* mpv,
                 BDRY* bdry,
                 const ElemSpaceDiscr* elem,
                 const NodeSpaceDiscr* node) 
{
    
    extern Thermodynamic th;
    extern User_Data ud;

    const double u0 = ud.wind_speed;
    const double v0 = 0.0;
    const double w0 = 0.0;
    const double delth = 0.0;
    
    const int icxe = elem->icx;
    const int icye = elem->icy;
    const int icze = elem->icz;
    const int igye = elem->igy;
    const int igze = elem->igz;
    
    const int icxn = node->icx;
    const int icyn = node->icy;
    const int iczn = node->icz;
    
    int i, j, k, l, m, n;
    double x, y, z;
    double rho, u, v, w, p, rhoY;	
        
    Hydrostatics_State(mpv, elem, node);
    
    for(k = igze; k < icze - igze; k++) {l = k * icxe * icye; 
        z = elem->z[k];
        
        for(j = igye; j < icye - igye; j++) {m = l + j * icxe;
                        
            y = elem->y[j];
            
            for(i = 0; i < icxe; i++) {n = m + i;
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
                Sol->rhoe[n] = rhoe(rho, u, v, w, p);
                Sol->rhoY[n] = rhoY;
                
                mpv->p2_cells[n]   = (p/rhoY) / ud.Msq;
                Sol->rhoX[BUOY][n] = Sol->rho[n] * ( Sol->rho[n]/Sol->rhoY[n] - mpv->HydroState->S0[j]);
		
            }			
        }                
    }  

    /*set nodal pressures */
    for(k = 0; k < iczn; k++) {l = k * icxn * icyn;   
        
        for(j = 0; j < icyn; j++) {m = l + j * icxn;                
            double p    = mpv->HydroState_n->p0[j];
            double rhoY = mpv->HydroState_n->rhoY0[j];
            
            for(i = 0; i < icxn; i++) {n = m + i;
                mpv->p2_nodes[n] = (p/rhoY) / ud.Msq;
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
        int is_compressible    = ud.is_compressible;
        double compressibility = ud.compressibility;
        ud.is_compressible = 0;
        ud.compressibility = 0.0;
        double *p2aux = (double*)malloc(node->nc*sizeof(double));
        for (int nn=0; nn<node->nc; nn++) {
            p2aux[nn] = mpv->p2_nodes[nn];
        }
        for (int nc=0; nc<elem->nc; nc++) {
            Sol->rhou[nc] -= u0*Sol->rho[nc];
            Sol->rhov[nc] -= v0*Sol->rho[nc];
        }
        
        //euler_backward_non_advective_expl_part(Sol, mpv, elem, ud.dtfixed);
        euler_backward_non_advective_impl_part(Sol, mpv, (const ConsVars*)Sol0, elem, node, 0.0, ud.dtfixed);
        for (int nn=0; nn<node->nc; nn++) {
            mpv->p2_nodes[nn] = p2aux[nn];
            mpv->dp2_nodes[nn] = 0.0;
        }
        free(p2aux);
        
        for (int nc=0; nc<elem->nc; nc++) {
            Sol->rhou[nc] += u0*Sol->rho[nc];
            Sol->rhov[nc] += v0*Sol->rho[nc];
        }
        
        ud.is_compressible = is_compressible;
        ud.compressibility = compressibility;
    }
}

/* ====================================================================== */ 

double pressure_function(double r, double p0, double S0, double u_theta, double Msq, double Gamma){
	return pow((pow(p0,Gamma) + Gamma*S0*Msq*u_theta*u_theta*(1.0 - pow((1.0-r),5.0)*(5.0*r+1.0))/30.0), (1.0/Gamma));
}

/* ====================================================================== */ 

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
