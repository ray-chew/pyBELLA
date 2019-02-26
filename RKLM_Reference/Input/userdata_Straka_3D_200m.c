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
	
	/* ========================================================================= */
	/* =====  PROBLEM SET UP  ================================================== */
	/* ========================================================================= */
	
    /* Earth */
	double grav  = 9.81;                               /* [m/s^2]                */
	double omega = 0.0;//2*PI*sin(0.25*PI)/(24.0*3600.0);    /*  [s^-1]                */
    
    /* thermodynamics and chemistry */
    double R_gas = 287.4;            /* [J/kg/K]                        */
    double R_vap = 461.00;           /* [J/kg/K]                        */
    double Q_vap = 2.53e+06;         /* [J]                             */
    double gamma = 1.4;              /* dimensionless                   */

    double viscm  = 75.0;            /* [m^2/s]                         */
    double viscbm = 0.0;             /* [m^2/s]                         */
    double visct  = 0.0;             /* [m^2/s]                         */
    double viscbt = 0.0;             /* [m^2/s]                         */
    double cond   = 75.0;            /* [m^2/s]                         */

    /* references for non-dimensionalization */
	double h_ref = 10000;            /* [m]                             */
	double t_ref = 100;              /* [s]                             */
	double T_ref = 300;              /* [K]                             */
	double p_ref = 100000;           /* [Pa]                            */
	double u_ref = h_ref/t_ref;      /* Strouhal No == 1 always assumed */
    double rho_ref  = p_ref / (R_gas*T_ref); /* [kg/m^3]          */

    double Nsq_ref  = grav*1.3e-05;     /* [] */
    
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
    ud->mol_trans   = STRAKA_DIFFUSION_MODEL; 
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
        
	/* flow domain */
	ud->xmin = - 1*2.56;
	ud->xmax =   1*2.56;
	ud->ymin =   0.0;
	ud->ymax =   0.64;
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
	
	/* ================================================================================== */
	/* =====  NUMERICS  ================================================================= */
	/* ================================================================================== */
	
    /* time discretization */
    ud->time_integrator       = SI_MIDPT;  /* this code version has only one option */
    ud->advec_time_integrator = STRANG; /* HEUN; EXPL_MIDPT;   best tested: STRANG; */
	ud->CFL                   = 0.96; /* 0.45; 0.9; 0.8; */
    ud->dtfixed0              = 0.0225;
	ud->dtfixed               = 0.0225; /* 0.0052; */ /*  0.004; */
    
    set_time_integrator_parameters(ud);
    
	/* Grid and space discretization */
	ud->inx = 257+1; /* 641; 321; 161; 129; 81; */
	ud->iny = 32+1; /* 321; 161;  81;  65; 41;  */
	ud->inz =  1;
	
	/* explicit predictor step */
	/* Recovery */
	ud->recovery_order        = SECOND;
	ud->limiter_type_scalars  = NONE; 
	ud->limiter_type_velocity = NONE; 
    /*  RUPE; NONE; MONOTONIZED_CENTRAL; MINMOD; VANLEER; SWEBY_MUNZ; SUPERBEE; */
	    
    /* parameters for SWEBY_MUNZ limiter family */
    ud->kp = 1.4;
	ud->kz = 1.4; /* Entro abused for velocity in split-step-aligned velocity ! */
	ud->km = 1.4;
	ud->kY = 1.4;
	ud->kZ = 1.4; /* 2.0 */
	
	ud->ncache =  333; /* (ud->inx+3); */
	
    /* linear solver-stuff */
    double tol                            = 1.e-6;
    ud->flux_correction_precision         = tol;
    ud->flux_correction_local_precision   = tol;    /* 1.e-05 should be enough */
    ud->second_projection_precision       = tol;
    ud->second_projection_local_precision = tol;  /* 1.e-05 should be enough */
    ud->flux_correction_max_iterations    = 6000;
    ud->second_projection_max_iterations  = 6000;
    ud->initial_projection                = WRONG;   /* WRONG;  CORRECT; */
    ud->initial_impl_Euler                = CORRECT; /* WRONG;  CORRECT; */
    
    ud->initial_projection                = WRONG;   /* If this is changed, implement div-control in Sol_initial() */
    ud->initial_impl_Euler                = WRONG;   /* to be tested: WRONG;  CORRECT; */

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
    ud->tout[0] =  3.0;
    ud->tout[1] =  6.0;
    ud->tout[2] =  9.0;
    ud->tout[3] = -1.0;
        
    ud->stepmax = 5000000;
    
	ud->write_stdout = ON;
	ud->write_stdout_period = 1;
	ud->write_file = ON;
	ud->write_file_period = 100000;
	ud->file_format = HDF;
    
    ud->n_time_series = 500; /* n_t_s > 0 => store_time_series_entry() called each timestep */
    
    {
        char *OutputBaseFolder      = "/home/tommaso/work/repos/RKLM_Reference/";
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
    
    const int compressible = (ud.is_compressible == 0 ? 0 : 1);

	const double u0 = ud.wind_speed;
	const double v0 = 0.0;
	const double w0 = 0.0;
	const double delT = 15.0/300.0; /* temp perturbation in K */
	
    const double x0 = 0.0;
    const double y0 = 0.3;
    const double rx = 0.4;
    const double ry = 0.2;
    		
	const int icx = elem->icx;
	const int icy = elem->icy;
	const int icz = elem->icz;
	const int igy = elem->igy;
    
    const int icxn = node->icx;
    const int icyn = node->icy;
	
	double x, y, xn, yn, u, v, w, p, rho, rhoY;
	
    States *HySt, *HyStn;
    double *Y, *Yn;

    double g;
	
	g = ud.gravity_strength[1];
	
    /* store background hydrostatic state in mpv auxiliary struct */
    Hydrostatics_State(mpv, elem, node);

    /* Aux space for cell- and node-based column-wise hydrostatics */
    HySt  = States_new(node->icy);
    Y     = (double*)malloc(node->icy*sizeof(double));
    HyStn = States_new(node->icy);
    Yn    = (double*)malloc(node->icy*sizeof(double));
    
    /* TODO: In contrast to Straka, I am using  delT  here as a 
     potential temperature perturbation, not as one of temperature
     itself
     */
    
    /* computations for the vertical slice at  k=0 */
    for(int i = 0; i < icx; i++) {
        
        /* set potential temperature stratification in the column */
        for(int j = 0; j < elem->icy; j++) {
            double dtemp, r;
            x     = elem->x[i];
            y     = elem->y[j];
            r     = sqrt((x-x0)*(x-x0)/(rx*rx) + (y-y0)*(y-y0)/(ry*ry));
            dtemp = delT * (r < 1.0 ? 0.5 * (1.0 + cos(PI*r)) : 0.0);
            Y[j]  = stratification(y) - dtemp/(ud.Msq*mpv->HydroState->p20[j]);
        }  
        
        for(int j = 0; j < node->icy; j++) {
            double dtemp, rn;
            xn    = node->x[i];
            yn    = node->y[j];
            rn    = sqrt((xn-x0)*(xn-x0)/(rx*rx) + (yn-y0)*(yn-y0)/(ry*ry));
            dtemp = delT * (rn < 1.0 ? 0.5 * (1.0 + cos(PI*rn)) : 0.0);
            Yn[j] = stratification(yn) - dtemp/(ud.Msq*mpv->HydroState_n->p20[j]);
        }    
        
        /* determine hydrostatic pressure distributions column-wise (lateral relation still neglected) */
        Hydrostatics_Column(HySt, HyStn, Y, Yn, elem, node);

        /* initialization of field variables */
        for(int j = igy; j < icy - igy + 1; j++) {
            
            int n  = j*icx + i;
            int nn   = j*icxn+i;
            
            u   = u0;
            v   = v0;
            w   = w0;
            
            p    = (compressible ? HySt->p0[j] : mpv->HydroState->p0[j]);            
            rhoY = (compressible ? HySt->rhoY0[j] : mpv->HydroState->rhoY0[j]);
            rho  = rhoY/Y[j];
            
            Sol->rho[n]    = rho;
            Sol->rhou[n]   = rho * u;
            Sol->rhov[n]   = rho * v;
            Sol->rhow[n]   = rho * w;
            Sol->rhoe[n]   = rhoe(rho, u, v, w, p);
            Sol->rhoY[n]   = rhoY;
            
            mpv->p2_cells[n]   = HySt->p20[j];
            Sol->rhoX[BUOY][n] = Sol->rho[n] * (Sol->rho[n]/Sol->rhoY[n] - mpv->HydroState->S0[j]);
            
            /* nodal pressure */
            mpv->p2_nodes[nn] = HyStn->p20[j];
        }
    }
    
    States_free(HySt);
    States_free(HyStn);
    free(Y);
    free(Yn);
    
    /* Find hydrostatic surface pressure and readjust pressure in the columns */
    Hydrostatic_Initial_Pressure(Sol, mpv, elem, node);  
    
    /* Copy the data just generated into all the parallel vertical slices */
    for (int j=0; j<icy; j++) {
        int mc = j*icx;
        int mn = j*icxn;
        for (int i=0; i<icx; i++) {
            int nc = mc + i;
            int nn = mn + i;
            for (int k=1; k<icz; k++) {
                int lc = nc + k*icx*icy;
                int ln = nn + k*icxn*icyn;
                
                Sol->rho[lc]    = Sol->rho[nc] ;
                Sol->rhou[lc]   = Sol->rhou[nc];
                Sol->rhov[lc]   = Sol->rhov[nc];
                Sol->rhow[lc]   = Sol->rhow[nc];
                Sol->rhoe[lc]   = Sol->rhoe[nc];
                Sol->rhoY[lc]   = Sol->rhoY[nc];
                
                Sol->rhoX[BUOY][lc] = Sol->rhoX[BUOY][nc];
                mpv->p2_cells[lc]   = mpv->p2_cells[nc];
                mpv->p2_nodes[ln]   = mpv->p2_nodes[nn];
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
