/*******************************************************************************
 File:   userdata.c
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
    
    /* references for non-dimensionalization */
    double h_ref    = 1;                 /* [m]               */
    double t_ref    = 1;                   /* [s]               */
    double T_ref    = 300.00;                /* [K]               */
    double p_ref    = 10e+5;                 /* [Pa]              */
    double u_ref    = h_ref/t_ref;           /* [m/s]; Sr = 1     */
    double rho_ref  = p_ref / (R_gas*T_ref); /* [kg/m^3]          */

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

	/* Low Mach */
    ud->is_nonhydrostatic = 1;
    ud->is_compressible   = 0;
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
	ud->xmin = - 0.5/ud->h_ref;  
	ud->xmax =   3.5/ud->h_ref;  
	ud->ymin = - 0.5/ud->h_ref;
	ud->ymax =   0.5/ud->h_ref; 
	ud->zmin = - 0.5/ud->h_ref/8.0;
	ud->zmax =   0.5/ud->h_ref/8.0;

	/* boundary/initial conditions */
	ud->wind_speed        = 1.0*1.0/ud->u_ref;              /* velocity in [u_ref] */
	ud->wind_shear        = -0.0;              /* velocity in [u_ref/h_ref] */             
	ud->hill_height       =  0.0;              /* height   in [h_ref]   */ 
	ud->hill_length_scale =  99999.9;          /* width    in [h_ref]   */   
	
	ud->bdrytype_min[0] = PERIODIC; /* DIRICHLET; */
	ud->bdrytype_min[1] = WALL; /* SLANTED_WALL; */
	ud->bdrytype_min[2] = WALL;
	ud->bdrytype_max[0] = PERIODIC; /* DIRICHLET; */  
	ud->bdrytype_max[1] = WALL;  
	ud->bdrytype_max[2] = WALL;
	
	ud->absorber = WRONG; /* CORRECT; */ 
	
	/* ======================================================================== */
	/* =====  NUMERICS  ======================================================= */
	/* ======================================================================== */
	
    /* time discretization */
    ud->time_integrator       = SI_MIDPT;
    ud->advec_time_integrator = STRANG; /* HEUN; EXPL_MIDPT;   default: STRANG;  */
	ud->CFL                   = 0.96;       
    ud->dtfixed0              = 10000.999;
    ud->dtfixed               = 10000.999;   
    
    set_time_integrator_parameters(ud);
    
	/* Grid and space discretization */
	ud->inx = 80+1; /*  */
	ud->iny = 20+1; /*  */
	ud->inz =   1;

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
	
    /* al explicit predictor operations are done on ncache-size data worms to save memory */ 
	ud->ncache =  201; /* (ud->inx+3); */
	
	/* linear solver-stuff */
    double tol = 1.e-10;
    ud->flux_correction_precision         = tol;
    ud->flux_correction_local_precision   = tol;    /* 1.e-05 should be enough */
    ud->second_projection_precision       = tol;
    ud->second_projection_local_precision = tol;  /* 1.e-05 should be enough */
    ud->flux_correction_max_iterations    = 6000;
    ud->second_projection_max_iterations  = 6000;
    
    ud->initial_projection                = WRONG;   /* to be tested: WRONG;  CORRECT; */
    ud->initial_impl_Euler                = WRONG;   /* to be tested: WRONG;  CORRECT; */
    
    ud->column_preconditioner             = WRONG; /* WRONG; CORRECT; */
    ud->synchronize_nodal_pressure        = WRONG;   /* WRONG; CORRECT; */
    ud->synchronize_weight                = 0.0;    /* relevant only when prev. option is "CORRECT"
                                                     Should ultimately be a function of dt . */  
            
	/* numerics parameters */
	ud->eps_Machine = sqrt(DBL_EPSILON);
		
	/* ================================================================================== */
    /* =====  CODE FLOW CONTROL  ======================================================== */
	/* ================================================================================== */
    ud->tout[0] =  0.5;      
    ud->tout[1] =  1.0;      
    ud->tout[2] =  1.5;      
    ud->tout[3] =  2.0;      
    ud->tout[4] =  2.5;      
    ud->tout[5] =  3.0;      
    ud->tout[6] = -1.0;
    /*
    ud->tout[0] =  3.00 * (ud->xmax-ud->xmin)/(10.0/ud->u_ref);      
    ud->tout[1] = -1.0;
     */

    ud->stepmax = 10000;

	ud->write_stdout = ON;
	ud->write_stdout_period = 1;
	ud->write_file = ON;
	ud->write_file_period = 100000;
	ud->file_format = HDF;

    ud->n_time_series = 500; /* n_t_s > 0 => store_time_series_entry() called each timestep */

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
    
	const double u0    = 1.0*ud.wind_speed;
	const double v0    = 0.0*ud.wind_speed;
	const double w0    = 0.0;
    
    const double urot = -1.0;  /* the origin of the March 24 - trouble ... ;^) */
    
    const double p0      = 1.0;
    const double rho0    = 1.0;  /* 0.5 standard;  1.0 stable configuration; */
    const double del_rho = 0.0;  /* 0.5 standard; -0.5 stable configuration; 0.0; for homentropic */
    const double R0      = 0.2;
    const double xc      = 0.0;
    const double yc      = 0.0;
    
    const int nhires   = 1;
    const int nhiressq = nhires*nhires;
    
    /*periodic setting: */
    const double xcm     = xc-(ud.xmax-ud.xmin);
    const double xcp     = xc+(ud.xmax-ud.xmin);
    const double ycm     = yc-(ud.ymax-ud.ymin);
    const double ycp     = yc+(ud.ymax-ud.ymin);
    
    const double dx = elem->dx;
    const double dy = elem->dy;
			
    const double dxx = dx/nhires;
    const double dyy = dy/nhires;
    
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
    double xcc, ycc;
	
    Hydrostatics_State(mpv, elem, node);
    
   /* Initial data and hydro-states in the flow domain */
	for(k = igz; k < icz - igz; k++) {
        l = k * icx * icy; 
		z = elem->z[k];
        
        for(j = igy; j < icy - igy; j++) {
            m = l + j * icx;
            y = elem->y[j];
			            
            for(i = 0; i < icx; i++) {
                
                double p2c = 0.0;

                n = m + i;                
                x = elem->x[i];

                Sol->rho[n]  = 0.0;
                Sol->rhou[n] = 0.0;
                Sol->rhov[n] = 0.0;
                Sol->rhow[n] = 0.0;
                Sol->rhoY[n] = 0.0;

                /* more accurate initialization by discrete integral averaging */
                for (int jj=0; jj<nhires; jj++) {
                    
                    double yy = (y-0.5*dy+0.5*dyy) + jj*dyy;

                    ycc = (fabs(yy-yc) < fabs(yy-ycm) ? (fabs(yy-yc) < fabs(yy-ycm) ? yc : ycp) : ycm);

                    for (int ii=0; ii<nhires; ii++) {
                                                
                        double xx = (x-0.5*dx+0.5*dxx) + ii*dxx;
                        
                        xcc = (fabs(xx-xc) < fabs(xx-xcm) ? (fabs(xx-xc) < fabs(xx-xcm) ? xc : xcp) : xcm);
                        
                        r       = sqrt((xx-xcc)*(xx-xcc) + (yy-ycc)*(yy-ycc));
                        uth     = urot * (r < R0 ? r/R0 : (r < 2.0*R0 ? 2.0 - r/R0 : 0.0));
                        
                        u       = u0 + uth * (-(yy-ycc)/r);
                        v       = v0 + uth * (+(xx-xcc)/r);
                        w       = w0;
                        p_hydro = mpv->HydroState->p0[j];
                        rhoY    = mpv->HydroState->rhoY0[j];
                        theta   = stratification(yy);
                        rho     =  (r < R0 ? (rho0 + del_rho*pow( 1-(r/R0)*(r/R0) , 6)) : rho0);
                        T       = T_from_p_rho(p_hydro,rho);
                        
                        if ( r/R0 < 1.0 ) {
                            p2c += rho*urot*urot * 0.5*(r/R0)*(r/R0);
                        }
                        else if ( r/R0 < 2.0 ) {
                            p2c += rho*urot*urot * (0.5 + (4.0*log(r/R0) - 4.0*(r/R0-1) + 0.5*((r/R0)*(r/R0) - 1.0)));
                        }
                        else {
                            p2c += rho*urot*urot * (0.5 + (4.0*log(2.0) - 4.0 + 1.5));;
                        }
                                                
                        Sol->rho[n]  += rho;
                        Sol->rhou[n] += rho * u;
                        Sol->rhov[n] += rho * v;
                        Sol->rhow[n] += rho * w;
                        
                        if (ud.is_compressible) {
                            double p     = p0 + ud.Msq*mpv->p2_cells[n];
                            Sol->rhoY[n] += pow(p,th.gamminv);
                            Sol->rhoe[n] += rhoe(rho, u, v, w, p);
                        } else {                    
                            Sol->rhoe[n] += rhoe(rho, u, v, w, p_hydro);
                            Sol->rhoY[n] += rhoY;
                        }
                    }
                }
                
                Sol->rho[n]  /= nhiressq;
                Sol->rhou[n] /= nhiressq;
                Sol->rhov[n] /= nhiressq;
                Sol->rhow[n] /= nhiressq;
                Sol->rhoY[n] /= nhiressq;
                
                mpv->p2_cells[n] = th.Gamma*(p2c/nhiressq)/mpv->HydroState->rhoY0[j];

            }            
		}                
	}  
    set_ghostcells_p2(mpv->p2_cells, elem, igx);

    /* nodal pressure */
    for(k = igzn; k < iczn - igzn; k++) {
        l = k * icxn * icyn; 
        z = node->z[k];
        
        for(j = igyn; j < icyn - igyn; j++) {
            m = l + j * icxn;
            y = node->y[j];
            
            ycc = (fabs(y-yc) < fabs(y-ycm) ? (fabs(y-yc) < fabs(y-ycm) ? yc : ycp) : ycm);
            
            for(i = igxn; i < icxn - igxn; i++) {
                n = m + i;                
                x       = node->x[i];
                xcc = (fabs(x-xc) < fabs(x-xcm) ? (fabs(x-xc) < fabs(x-xcm) ? xc : xcp) : xcm);
                
                r = sqrt((x-xcc)*(x-xcc) + (y-ycc)*(y-ycc));
                
                if ( r/R0 < 1.0 ) {
                    mpv->p2_nodes[n] = rho0*u0*u0 * 0.5*(r/R0)*(r/R0);
                }
                else if ( r/R0 < 2.0 ) {
                    mpv->p2_nodes[n] =  rho0*urot*urot * (0.5 + (4.0*log(r/R0) - 4.0*(r/R0-1) + 0.5*((r/R0)*(r/R0) - 1.0)));
                }
                else {
                    mpv->p2_nodes[n] = rho0*urot*urot * (0.5 + (4.0*log(2.0) - 4.0 + 1.5));
                }
                
                /* Exner pressure */
                mpv->p2_nodes[n]  = th.Gamma*mpv->p2_nodes[n]/mpv->HydroState->rhoY0[j];
            }            
        }                
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
