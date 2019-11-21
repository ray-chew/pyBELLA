/*******************************************************************************
 File:   userdata.c
 *******************************************************************************/
#include <math.h>
#include <float.h>
#include <stdlib.h>

#include "enumerator.h"
#include "Common.h"
#include "userdata.h"
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

void User_Data_init(User_Data* ud) {
	
	int i;
	
	/* ================================================================================== */
	/* =====  PROBLEM SET UP  =========================================================== */
	/* ================================================================================== */
	
    /* Earth */	
	double grav  = 0.0; /* [m/s^2]                                 */
	double omega = 0.0; /* 2*PI*sin(0.25*PI)/(24.0*3600.0); [s^-1] */

    /* thermodynamics and chemistry */
    double R_gas    = 287.4;            /* [J/kg/K]               */
    double R_vap    = 461.00;           /* [J/kg/K]               */
    double Q_vap    = 2.53e+06;         /* [J]                    */
    double gamma    = 1.4;              /* dimensionless          */
    
    double viscm  = 0.0;             /* [m^2/s]                         */
    double viscbm = 0.0;             /* [m^2/s]                         */
    double visct  = 0.0;             /* [m^2/s]                         */
    double viscbt = 0.0;             /* [m^2/s]                         */
    double cond   = 0.0;             /* [m^2/s]                         */

    /* references for non-dimensionalization */
    double h_ref    = 10000;                 /* [m]               */
    double t_ref    = 100;                   /* [s]               */
    double T_ref    = 300.00;                /* [K]               */
    double p_ref    = 1e+5;                 /* [Pa]              */
    double u_ref    = h_ref/t_ref;           /* [m/s]; Sr = 1     */
    double rho_ref  = p_ref / (R_gas*T_ref); /* [kg/m^3]          */

    /* reference stratification as (buoyancy frequency)^2 */
    double Nsq_ref  = 0.0;                   /* [1/s^2]           */

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
    ud->is_nonhydrostatic = 1;
    ud->is_compressible   = 1;
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
    
	/* flow domain */
    ud->xmin = - 0.5;  
    ud->xmax =   0.5;  
    ud->ymin = - 1.0;
    ud->ymax =   1.0; 
    ud->zmin = - 0.5;
    ud->zmax =   0.5;

	/* boundary/initial conditions */
	ud->wind_speed        =  1.0;              /* velocity in [u_ref] */
	ud->wind_shear        = -0.0;              /* velocity in [u_ref/h_ref] */             
    ud->hill_shape        = AGNESI;            /* AGNESI, SCHLUTOW */
	ud->hill_height       =  0.0;              /* height   in [h_ref]   */ 
	ud->hill_length_scale =  99999.9;          /* width    in [h_ref]   */   
	
	ud->bdrytype_min[0] = WALL; /* DIRICHLET; */
	ud->bdrytype_min[1] = PERIODIC; /* SLANTED_WALL; */
	ud->bdrytype_min[2] = PERIODIC;
	ud->bdrytype_max[0] = WALL; /* DIRICHLET; */  
	ud->bdrytype_max[1] = PERIODIC;  
	ud->bdrytype_max[2] = PERIODIC;
	
	ud->absorber        = WRONG; /* CORRECT; */ 
    ud->bottom_theta_bc = BOTTOM_BC_DEFAULT;

	/* ======================================================================== */
	/* =====  NUMERICS  ======================================================= */
	/* ======================================================================== */
	
    /* time discretization */
    ud->time_integrator       = SI_MIDPT;  /* this code version has only one option */
    ud->advec_time_integrator = STRANG; /* HEUN; EXPL_MIDPT;   best tested: STRANG; */
	ud->CFL                   = 0.9/2.0;       
    ud->dtfixed0              = 2.1*1.200930e-02;
    ud->dtfixed               = 2.1*1.200930e-02;   
    
    set_time_integrator_parameters(ud);
    
	/* Grid and space discretization */
	ud->inx =   3+1; /*  */
	ud->iny =  96+1; /*  */
	ud->inz =  48+1;

    /* explicit predictor step */
	/* Recovery */
	ud->recovery_order        = SECOND;
	ud->limiter_type_scalars  = NONE; 
	ud->limiter_type_velocity = NONE; 
    /*  RUPE; NONE; MONOTONIZED_CENTRAL; MINMOD; VANLEER; SWEBY_MUNZ; SUPERBEE; */
        
    /* parameters for SWEBY_MUNZ limiter family */
    ud->kp = 0.0; /* 1.4; */
    ud->kz = 0.0; /* 1.4; */
    ud->km = 0.0; /* 1.4; */
    ud->kY = 0.0; /* 1.4; */
    ud->kZ = 0.0; /* 1.4; */
	
    /* all explicit predictor operations are done on ncache-size data worms to save memory */ 
	ud->ncache =  201; /* (ud->inx+3); */
	
	/* linear solver-stuff */
    double tol = 1.e-10;
    ud->flux_correction_precision         = tol;
    ud->flux_correction_local_precision   = tol;    /* 1.e-05 should be enough */
    ud->second_projection_precision       = tol;
    ud->second_projection_local_precision = tol;  /* 1.e-05 should be enough */
    ud->flux_correction_max_iterations    = 6000;
    ud->second_projection_max_iterations  = 6000;
            
    ud->initial_projection                = CORRECT;   /* to be tested: WRONG;  CORRECT; */
    ud->initial_impl_Euler                = WRONG;   /* to be tested: WRONG;  CORRECT; */
    
    ud->column_preconditioner             = WRONG; /* WRONG; CORRECT; */

	/* numerics parameters */
	ud->eps_Machine = sqrt(DBL_EPSILON);
		
	/* ================================================================================== */
    /* =====  CODE FLOW CONTROL  ======================================================== */
	/* ================================================================================== */
    
    ud->tout[0] =  1.0;      
    //ud->tout[1] =  2.0;      
    //ud->tout[2] =  3.0;      
    ud->tout[1] = -1.0;
    
    /*
     ud->tout[0] =  0.5;      
     ud->tout[1] =  1.0;      
     ud->tout[2] =  1.5;      
     ud->tout[3] =  2.0;      
     ud->tout[4] =  2.5;      
     ud->tout[5] =  3.0;      
     ud->tout[6] = -1.0;
     */
    
    ud->stepmax = 20000;

	ud->write_stdout = ON;
	ud->write_stdout_period = 1;
	ud->write_file = ON;
	ud->write_file_period = 20;
	ud->file_format = HDF;

    ud->n_time_series = 500; /* n_t_s > 0 => store_time_series_entry() called each timestep */

    {
#ifdef TOMMASO
        char *OutputBaseFolder      = "/home/tommaso/work/repos/RKLM_Reference/";
#else
        char *OutputBaseFolder      = "/Users/rupert/Documents/Computation/RKLM_Reference/";
#endif
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
    extern double* diss_midpnt;
    
	const double u0    = 0.0;
	const double v0    = 1.0*ud.wind_speed;
	const double w0    = 1.0*ud.wind_speed;
    
    const double rotdir = 1.0;
    
    const double p0      = 1.0;
    const double rho0    = 0.5;  /* 0.5 standard;  1.0 stable configuration; */
    const double del_rho = 0.5;  /* 0.5 standard; -0.5 stable configuration; 0.0; for homentropic */
    const double R0      = 0.4;
    const double fac     = 1*1024.0; /* 4*1024.0 */
    const double yc      = 0.0;
    const double zc      = 0.0;
    
    /*periodic setting: */
    const double ycm     = yc-(ud.ymax-ud.ymin);
    const double ycp     = yc+(ud.ymax-ud.ymin);
    const double zcm     = zc-(ud.zmax-ud.zmin);
    const double zcp     = zc+(ud.zmax-ud.zmin);
			
	const int icx  = elem->icx;
	const int icy  = elem->icy;
	const int icz  = elem->icz;
    const int igx  = elem->igx;
	const int igy  = elem->igy;
	const int igz  = elem->igz;
    
    const int icxn  = node->icx;
    const int icyn  = node->icy;
    const int iczn  = node->icz;
    const int igxn  = node->igx;
    const int igyn  = node->igy;
    const int igzn  = node->igz;

	int i, j, k, l, m, n;
	double x, y, z;
	double rho, u, v, w, rhoY, theta, T, p_hydro;
    double r, uth;
    double ycc, zcc;
    
    /* data needed for the pressure distribution in the Kadioglu et al. paper */ 
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

    Hydrostatics_State(mpv, elem, node);
    
    /* Initial data and hydro-states in the flow domain */
	for(k = igz; k < icz - igz; k++) {
        l = k * icx * icy; 
		z = elem->z[k];
        
        zcc = (fabs(z-zc) < fabs(z-zcm) ? (fabs(z-zc) < fabs(z-zcm) ? zc : zcp) : zcm);

        for(j = igy; j < icy - igy; j++) {
            m = l + j * icx;
            y = elem->y[j];

            ycc = (fabs(y-yc) < fabs(y-ycm) ? (fabs(y-yc) < fabs(y-ycm) ? yc : ycp) : ycm);

            for(i = 0; i < icx; i++) {
                
                double p2c = 0.0;
                double dp2c = 0.0;

                n = m + i;                
                x       = elem->x[i];

                r       = sqrt((y-ycc)*(y-ycc) + (z-zcc)*(z-zcc));
                uth     = rotdir * (r < R0 ? fac * pow( 1.0-r/R0, 6) * pow( r/R0, 6) : 0.0);
                
                u       = u0;
                v       = v0 + uth * (+(z-zcc)/r);
                w       = w0 + uth * (-(y-ycc)/r);
                p_hydro = mpv->HydroState->p0[j];
                rhoY    = mpv->HydroState->rhoY0[j];
                theta   = stratification(y);
                rho     =  (r < R0 ? (rho0 + del_rho*pow( 1-(r/R0)*(r/R0) , 6)) : rho0);
                T       = T_from_p_rho(p_hydro,rho);
                
                dp2c = 0.0;
                if ( r/R0 < 1.0 ) {
                    for (int ip = 0; ip < 25; ip++)
                    {
                        dp2c += (rho0+del_rho) * coe[ip] * (pow(r/R0 ,12+ip) - 1.0) * rotdir * rotdir;
                    }
                }
                
                p2c += dp2c;
                Sol->rho[n]  = rho;
                Sol->rhou[n] = rho * u;
                Sol->rhov[n] = rho * v;
                Sol->rhow[n] = rho * w;

                if (ud.is_compressible) {
                    double p     = p0 + ud.Msq*fac*fac*dp2c;
                    Sol->rhoY[n] += pow(p,th.gamminv);
                    Sol->rhoe[n] += rhoe(rho, u, v, w, p);
                } else {                    
                    Sol->rhoe[n] += rhoe(rho, u, v, w, p_hydro);
                    Sol->rhoY[n] += rhoY;
                }

                /* Exner pressure */
                mpv->p2_cells[n]  = th.Gamma*fac*fac*p2c/mpv->HydroState->rhoY0[j];
            }            
		}                
	}  
    set_ghostcells_p2(mpv->p2_cells, elem, igx);

    /* nodal pressure */
    for(k = igzn; k < iczn - igzn; k++) {
        l = k * icxn * icyn; 
        z = node->z[k];
        
        zcc = (fabs(z-zc) < fabs(z-zcm) ? (fabs(z-zc) < fabs(z-zcm) ? zc : zcp) : zcm);

        for(j = igyn; j < icyn - igyn; j++) {
            m = l + j * icxn;
            y = node->y[j];
                        
            ycc = (fabs(y-yc) < fabs(y-ycm) ? (fabs(y-yc) < fabs(y-ycm) ? yc : ycp) : ycm);

            for(i = igxn; i < icxn - igxn; i++) {
                n = m + i;                
                x       = node->x[i];
                
                r       = sqrt((y-ycc)*(y-ycc) + (z-zcc)*(z-zcc));
                
                if ( r/R0 < 1.0 ) {
                    mpv->p2_nodes[n] = 0.0;
                    for (int ip = 0; ip < 25; ip++)
                    {
                        mpv->p2_nodes[n] += coe[ip] * (pow(r/R0 ,12+ip) - 1.0) * rotdir * rotdir;
                    }
                }
                else {
                    mpv->p2_nodes[n] = 0.0;
                }
                
                /* Exner pressure */
                mpv->p2_nodes[n]  = th.Gamma*fac*fac*mpv->p2_nodes[n]/mpv->HydroState->rhoY0[j];
            }            
        }                
    }
    
    ud.nonhydrostasy   = nonhydrostasy(0);
    ud.compressibility = compressibility(0);
    
    set_wall_rhoYflux(bdry, Sol, mpv, elem);
    Set_Explicit_Boundary_Data(Sol, elem, OUTPUT_SUBSTEPS);
    
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
            Sol->rhow[nc] -= w0*Sol->rho[nc];
        }
        
        //euler_backward_non_advective_expl_part(Sol, mpv, elem, ud.dtfixed);
        euler_backward_non_advective_impl_part(Sol, mpv, diss_midpnt, (const ConsVars*)Sol0, elem, node, 0.0, ud.dtfixed, 1.0);
        for (int nn=0; nn<node->nc; nn++) {
            mpv->p2_nodes[nn] = p2aux[nn];
            mpv->dp2_nodes[nn] = 0.0;
        }
        free(p2aux);
        
        for (int nc=0; nc<elem->nc; nc++) {
            Sol->rhou[nc] += u0*Sol->rho[nc];
            Sol->rhov[nc] += v0*Sol->rho[nc];
            Sol->rhow[nc] += w0*Sol->rho[nc];
        }
        
        ud.is_compressible = is_compressible;
        ud.compressibility = compressibility;
    }
    

}

/* ================================================================================== */

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
