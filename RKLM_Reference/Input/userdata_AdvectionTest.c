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

/* mollification function to set a localized potential temperature perturbation for advection */
double molly(double r, const double a);

void User_Data_init(User_Data* ud) {
    
    int i; 
    
    /* ================================================================================== */
    /* =====  PROBLEM SET UP  =========================================================== */
    /* ================================================================================== */
    
    /* Earth */
    double grav     = 0.0;              /* gravitational acceleration [m/s^2]    */
    double omega    = 0.0;              /* Coriolis parameter [1/s]              */
                                        /* sin(0.5*PI) * 2.0 * 0.00007272205217; */
    
    /* thermodynamics and chemistry */
    double R_gas    = 287.4;            /* [J/kg/K]               */
    double R_vap    = 461.00;           /* [J/kg/K]               */
    double Q_vap    = 2.53e+06;         /* [J]                    */
    double gamma    = 1.4;              /* dimensionless          */

    /* references for non-dimensionalization */
    double h_ref    = 10000;                 /* [m]               */
    double t_ref    = 100;                   /* [s]               */
    double T_ref    = 300.00;                /* [K]               */
    double p_ref    = 10e+5;                 /* [Pa]              */
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

    /* Low Mach */
    ud->is_compressible   = 1;    /* 0: psinc;  1: comp;  -1: transition (see compressibility()) */
    ud->acoustic_timestep = 0;    /* advective time step -> 0;  acoustic time step -> 1; */
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
    ud->xmin = - 1.0;
    ud->xmax =   1.0;
    ud->ymin = - 1.0;
    ud->ymax =   1.0;
    ud->zmin = - 1.0;
    ud->zmax =   1.0;
    
    /* boundary/initial conditions */
    ud->wind_speed  =  1.0 * 20.0/u_ref;
    ud->wind_shear  = -0.0;              /* velocity in [u_ref/h_ref] */
    ud->hill_height = 0.0 * 0.096447; 
    ud->hill_length_scale = 0.1535;   /* hill_length * l_ref = 1.0 km */
    
    ud->bdrytype_min[0] = PERIODIC; /* DIRICHLET; PERIODIC; WALL; */
    ud->bdrytype_min[1] = PERIODIC; /* SLANTED_WALL; */
    ud->bdrytype_min[2] = PERIODIC;
    ud->bdrytype_max[0] = PERIODIC; /* DIRICHLET; PERIODIC; WALL; */
    ud->bdrytype_max[1] = PERIODIC;
    ud->bdrytype_max[2] = PERIODIC;
    
    ud->absorber = WRONG; /* CORRECT;  WRONG; */ /*  BREAKING WAVE CHANGE */
    
    /* ================================================================================== */
    /* =====  NUMERICS  ================================================================= */
    /* ================================================================================== */
    
    /* time discretization */
    ud->time_integrator       = SI_MIDPT; /* this code version has only one option */
    ud->advec_time_integrator = STRANG; /* HEUN; EXPL_MIDPT;   best tested: STRANG; */
    ud->CFL                   = 0.96; /* 0.45; 0.9; 0.8; */
    ud->dtfixed0              = 9999.9;
    ud->dtfixed               = 9999.9;
    
    set_time_integrator_parameters(ud);
    
    /* Grid and space discretization */
    ud->inx =  31+1; /* 641; 321; 161; 129; 81; */
    ud->iny =  31+1; /* 321; 161;  81;  65; 41;  */
    ud->inz =  1;
    
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
    
    /* numerics parameters */
    ud->eps_Machine = sqrt(DBL_EPSILON);
    
    /* ================================================================================== */
    /* =====  CODE FLOW CONTROL  ======================================================== */
    /* ================================================================================== */
    
    ud->tout[0] = 3000.0 / ud->t_ref;
    ud->tout[1] = -1.0;

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

void Sol_initial(ConsVars* Sol, 
                 const ElemSpaceDiscr* elem, 
                 const NodeSpaceDiscr* node) {
    
    extern Thermodynamic th;
    extern User_Data ud;
    extern MPV* mpv;
    
    const int compressible = (ud.is_compressible == 0 ? 0 : 1);
    
    
    /* advection test */
    const double phi   = 2.0*PI/4.0;
    const double theta = 2.0*PI/4.0;
    const double u0    = cos(phi)*cos(theta)*ud.wind_speed;
    const double v0    = sin(phi)*cos(theta)*ud.wind_speed;
    const double w0    = sin(theta)*ud.wind_speed;
    const double delth = 0.25; /* pot. temp. perturbation amplitude; standard:  0.01 / ud.T_ref */
    const double xc    = 0.0;  /* initial position of center of pot temp perturbation */
    const double yc    = 0.0;  /* initial position of center of pot temp perturbation */
    const double zc    = 0.0;  /* initial position of center of pot temp perturbation */
    const double a     = 0.25; /* characteristic width of the witch of Agnesi type mollifier */
    
    const int icx = elem->icx;
    const int icy = elem->icy;
    const int icz = elem->icz;
    const int igx = elem->igx;
    const int igy = elem->igy;
    
    const int icxn = node->icx;
    const int icyn = node->icy;
    
    double rho, u, v, w, p, rhoY;
    
    double g;
        
    /* This will become a 3Dified 2D test case. So I first compute
     the 2D vertical slice data and then simply copy them into all
     the vertical slices for different z values.
     (y is the vertical coordinate!)
     */
    
    g = ud.gravity_strength[1];

    /* store background hydrostatic state in mpv auxiliary struct */
    Hydrostatics_State(mpv, elem, node);

    /* store background hydrostatic state in mpv auxiliary struct */
    Hydrostatics_State(mpv, elem, node);
        
    /* initialization of field variables */
    for(int k = 0; k < icz; k++) {
        int lc = k*icx*icy;
        for(int j = 0; j < icy; j++) {
            int mc = lc + j*icx;
            for(int i = 0; i < icx; i++) {
                int nc = mc + i;
                
                double x = elem->x[i];
                double y = elem->y[j];
                double z = elem->z[k];
                double Y = 1.0 + delth*molly(sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc)+(z-zc)*(z-zc)), a);
                
                double p    = mpv->HydroState->p0[j];            
                double rhoY = mpv->HydroState->rhoY0[j];
                double rho  = rhoY/Y;
                
                Sol->rho[nc]    = rho;
                Sol->rhou[nc]   = rho * u0;
                Sol->rhov[nc]   = rho * v0;
                Sol->rhow[nc]   = rho * w0;
                Sol->rhoe[nc]   = rhoe(rho, u0, v0, w0, p);
                Sol->rhoY[nc]   = rhoY;
                
                mpv->p2_cells[nc]   = 0.0;
                Sol->rhoX[BUOY][nc] = Sol->rho[nc] * ( Sol->rho[nc]/Sol->rhoY[nc] - mpv->HydroState->S0[j]);
            }
        }
    }
    for (int nn=0; nn < node->nc; nn++) {
        mpv->p2_nodes[nn] = 0.0;
    }
}

/* ================================================================================== */

double stratification(
                      double y) {
    
    extern User_Data ud;
    
    double Nsq = ud.Nsq_ref * ud.t_ref * ud.t_ref;
    double g   = ud.gravity_strength[1] / ud.Msq;
    
    return( exp(Nsq*y/MAX_own(ud.eps_Machine, g)) );
}


/* ================================================================================== */

double molly(
             double r,
             const double a) {
    
    extern User_Data ud;
        
    return(0.5*(1.0+cos(MIN_own(0.5*PI*r/a,PI))));
}



/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: userdata.c,v $
 Revision 1.2  1998/03/07 09:56:47  nicola
 Added flux computation and multiple pressure variables.
 
 Revision 1.1  1998/03/01 18:50:44  nicola
 Added csr.c and sod1d.c (user data for 3d code)
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
