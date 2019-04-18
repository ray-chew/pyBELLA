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

/* mollification function used in the Skamarock-Klemp 1994 tests for theta perturbations */
double molly(double x);

/* horizontal stretch for S&K94 IGWs: planetary -> 160.0;  long-wave -> 20.0;  standard -> 1.0; */

void User_Data_init(User_Data* ud) {
    
    int i; 
    
    /* ================================================================================== */
    /* =====  PROBLEM SET UP  =========================================================== */
    /* ================================================================================== */
    
    /* Earth */
    double grav     = 9.81;             /* gravitational acceleration [m/s^2]    */
    double omega    = 0.0*0.0001;       /* Coriolis parameter [1/s]              */
                                        /* sin(0.5*PI) * 2.0 * 0.00007272205217; */
    
    /* thermodynamics and chemistry */
    double R_gas    = 287.4;            /* [J/kg/K]               */
    double R_vap    = 461.00;           /* [J/kg/K]               */
    double Q_vap    = 2.53e+06;         /* [J]                    */
    double gamma    = 1.4;              /* dimensionless; 5.0/3.0       */

    /* references for non-dimensionalization */
    double h_ref    = 10000;                 /* [m]               */
    double t_ref    = 100;                   /* [s]               */
    double T_ref    = 222.65;                /* [K]               */
    double p_ref    = 0.237e+5;                  /* [Pa]              */
    double u_ref    = h_ref/t_ref;           /* [m/s]; Sr = 1     */
    double rho_ref  = p_ref / (R_gas*T_ref); /* [kg/m^3]          */
    
    /* reference stratification as (buoyancy frequency)^2 */
    double Nsq_ref  = 1.0e-4;                /* [1/s^2]           */
    
    ud->h_ref       = h_ref;
    ud->t_ref       = t_ref;
    ud->T_ref       = T_ref;
    ud->p_ref       = p_ref;
    ud->u_ref       = u_ref;
    ud->rho_ref     = rho_ref;
    ud->Nsq_ref     = Nsq_ref;
    ud->g_ref       = grav;
    ud->gamm        = gamma;
    ud->Rg_over_Rv  = R_gas/R_vap;
    ud->Q           = Q_vap/(R_gas*T_ref);
    
    /* number of advected species */
    ud->nspec       = NSPEC;

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
        
    /* flow domain; all lengths in units of  href  */
    ud->xmin = -20.0e+6/ud->h_ref;
    ud->xmax =  20.0e+6/ud->h_ref;
    ud->ymin =      0.0/ud->h_ref;
    ud->ymax =   3.0e+5/ud->h_ref;
    ud->zmin = - 1.0;
    ud->zmax =   1.0;
    
    /* boundary/initial conditions */
    ud->wind_speed  =  1.0 * 20.0/u_ref;
    ud->wind_shear  = -0.0;              /* velocity in [u_ref/h_ref] */
    ud->hill_shape  = AGNESI;            /* AGNESI, SCHLUTOW */
    ud->hill_height = 0.0 * 0.096447; 
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
    ud->CFL                   = 0.9; /* 0.45; 0.9; 0.8; */
    ud->dtfixed0              = 2000 / ud->t_ref;
    ud->dtfixed               = 2000 / ud->t_ref;
    
    set_time_integrator_parameters(ud);
    
    /* Grid and space discretization */
    ud->inx =  300+1; /* 641; 321; 161; 129; 81; */
    ud->iny =   30+1; /* 321; 161;  81;  65; 41;  */
    ud->inz =      1;
    
    /* explicit predictor step */
    /* Recovery */
    ud->recovery_order = SECOND; /* FIRST, SECOND */ 
    ud->limiter_type_scalars  = NONE; 
    ud->limiter_type_velocity = NONE; 
    /*  RUPE; NONE; MONOTONIZED_CENTRAL; MINMOD; VANLEER; SWEBY_MUNZ; SUPERBEE; */
    
    /* parameters for SWEBY_MUNZ limiter family */
    ud->kp = 1.4;
    ud->kz = 1.4; /* Entro abused for velocity in split-step-aligned velocity ! */
    ud->km = 1.4;
    ud->kY = 1.4;
    ud->kZ = 1.4; /* 2.0 */
        
    /* al explicit predictor operations are done on ncache-size data worms to save memory */ 
    ud->ncache = 75; /* 71+4; 304*44; 604*44; (ud->inx+3); (ud->inx+3)*(ud->iny+3);*/
    
    /* linear solver-stuff */
    double tol                            = 1.e-10;
    ud->flux_correction_precision         = tol;
    ud->flux_correction_local_precision   = tol;    /* 1.e-05 should be enough */
    ud->second_projection_precision       = tol;
    ud->second_projection_local_precision = tol;  /* 1.e-05 should be enough */
    ud->flux_correction_max_iterations    = 6000;
    ud->second_projection_max_iterations  = 6000;
    ud->initial_projection                = WRONG;   /* WRONG;  CORRECT; */
    ud->initial_impl_Euler                = CORRECT; /* WRONG;  CORRECT; */
    
    ud->column_preconditioner             = CORRECT; /* WRONG; CORRECT; */
    ud->synchronize_nodal_pressure        = WRONG;   /* WRONG; CORRECT; */
    ud->synchronize_weight                = 0.0;    /* relevant only when prev. option is "CORRECT"
                                                     Should ultimately be a function of dt . */  

    /* numerics parameters */
    ud->eps_Machine = sqrt(DBL_EPSILON);
    
    /* ================================================================================== */
    /* =====  CODE FLOW CONTROL  ======================================================== */
    /* ================================================================================== */
    
    ud->tout[0] = 45000.0 / ud->t_ref; /* 3000 */
    ud->tout[1] = -1.0;

    ud->stepmax = 10000;
    
    ud->write_stdout = ON;
    ud->write_stdout_period = 1;
    ud->write_file = ON;
    ud->write_file_period = 40;
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

void Sol_initial(ConsVars* Sol, const ElemSpaceDiscr* elem, const NodeSpaceDiscr* node) {
    
    extern Thermodynamic th;
    extern User_Data ud;
    extern MPV* mpv;
    
    const int compressible = (ud.is_compressible == 0 ? 0 : 1);
    
    
    /* Skamarock-Klemp 1994 - type test */
    const double u0    = ud.wind_speed;
    const double v0    = 0.0;
    const double w0    = 0.0;
    const double delth = 0.01 / ud.T_ref;                    /* pot. temp. perturbation amplitude; standard:  0.01 / ud.T_ref */
    const double xc    = -15.0e+6/ud.h_ref; /* initial position of center of pot temp perturbation */
    const double a     =   6.0e+5/ud.h_ref;       /* characteristic width of the witch of Agnesi type mollifier */
    
    const int icx = elem->icx;
    const int icy = elem->icy;
    const int icz = elem->icz;
    const int igy = elem->igy;
    
    const int icxn = node->icx;
    const int icyn = node->icy;
    
    int i, j, n, nn;
    double x, y, xn, yn;
    double rho, u, v, w, p, rhoY;
        
    States *HySt, *HyStn;
    double *Y, *Yn;
    
    /* This will become a 3Dified 2D test case. So I first compute
     the 2D vertical slice data and then simply copy them into all
     the vertical slices for different z values.
     (y is the vertical coordinate!)
     */
        
    /* store background hydrostatic state in mpv auxiliary struct */
    Hydrostatics_State(mpv, elem, node);
    
    /* Aux space for cell- and node-based column-wise hydrostatics */
    HySt  = States_new(node->icy);
    Y     = (double*)malloc(node->icy*sizeof(double));
    HyStn = States_new(node->icy);
    Yn    = (double*)malloc(node->icy*sizeof(double));
     
    /* computations for the vertical slice at  k=0 */
    for(i = 0; i < icx; i++) {
        
        /* set potential temperature stratification in the column */
        for(j = 0; j < elem->icy; j++) {
            x     = elem->x[i];
            y     = elem->y[j];
            Y[j]  = stratification(y)  + delth * molly(x) * sin(PI*y/ud.ymax)  / (1.0 + (x-xc)*(x-xc) / (a*a));
        }  
        
        for(j = 0; j < node->icy; j++) {
            xn    = node->x[i];
            yn    = node->y[j];
#if 0
            Yn[j] = stratification(yn)  + delth * sin(PI*yn/ud.ymax);
#else
            Yn[j] = stratification(yn)  + delth * molly(xn) * sin(PI*yn/ud.ymax)  / (1.0 + (xn-xc)*(xn-xc) / (a*a));
#endif
        }        
        /* determine hydrostatic pressure distributions column-wise (lateral relation still neglected) */
        Hydrostatics_Column(HySt, HyStn, Y, Yn, elem, node);
        
        /* initialization of field variables */
        for(j = igy; j < icy - igy + 1; j++) {
            
            n  = j*icx + i;
                        
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
            nn   = j*icxn+i;
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
}

/* ================================================================================== */

double stratification(
                      double y) {
    
    extern User_Data ud;
    
    double Hinv   = 4.08163e-5 * ud.h_ref;
    double Theta0 = 222.65 / ud.T_ref;
    return( Theta0*exp(Hinv*y) );
}


/* ================================================================================== */

double molly(
             double x) {
    
    extern User_Data ud;
    
    const double del  = 0.25;
    const double L    = ud.xmax-ud.xmin;
    const double xi_l = MIN_own(1.0, (x-ud.xmin)/(del*L));
    const double xi_r = MIN_own(1.0, (ud.xmax-x)/(del*L));
    
    return(0.5* MIN_own(1.0-cos(PI*xi_l), 1.0-cos(PI*xi_r)));
}



/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: userdata.c,v $
 Revision 1.2  1998/03/07 09:56:47  nicola
 Added flux computation and multiple pressure variables.
 
 Revision 1.1  1998/03/01 18:50:44  nicola
 Added csr.c and sod1d.c (user data for 3d code)
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
