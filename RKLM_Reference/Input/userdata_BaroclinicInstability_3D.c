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

/* ====================================================================================== */

double molly(double x);
double F(double z, double dth);
double fH(double z, double dth);
double HT(double z);
double Zeta(double y, double z);
double Thetat(double y, double z);
double Thetas(double y, double z);
double Thetae(double y, double z);

/* ====================================================================================== */

void User_Data_init(User_Data* ud) {
    
    int i; 
    
    /* ================================================================================== */
    /* =====  PROBLEM SET UP  =========================================================== */
    /* ================================================================================== */
    
    /* Earth */
    double grav     = 9.81;             /* gravitational acceleration [m/s^2]    */
    double omega    = 1.0*0.0001;       /* Coriolis parameter [1/s]              */
    /* sin(0.5*PI) * 2.0 * 0.00007272205217; */
    
    /* thermodynamics and chemistry */
    double R_gas    = 287.4;            /* [J/kg/K]               */
    double R_vap    = 461.00;           /* [J/kg/K]               */
    double Q_vap    = 2.53e+06;         /* [J]                    */
    double gamma    = 1.4;              /* dimensionless; 5.0/3.0       */
    
    double viscm  = 0.0;            /* [m^2/s]                         */
    double viscbm = 0.0;             /* [m^2/s]                         */
    double visct  = 0.0;             /* [m^2/s]                         */
    double viscbt = 0.0;             /* [m^2/s]                         */
    double cond   = 0.0;            /* [m^2/s]                         */
    
    /* references for non-dimensionalization */
    double h_ref    = 1000;                  /* [m]               */
    double t_ref    = 100;                   /* [s]               */
    double T_ref    = 300.00;                /* [K]               */
    double p_ref    = 1e+5;                  /* [Pa]              */
    double u_ref    = h_ref/t_ref;           /* [m/s]; Sr = 1     */
    double rho_ref  = p_ref / (R_gas*T_ref); /* [kg/m^3]          */
    
    /* reference stratification as (buoyancy frequency)^2 */
    double Nsq_ref  = 0.00011772;                /* [1/s^2]           */
    
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
    ud->is_ArakawaKonor   =  0;    /* 1: AK09;  0: model determined by is_nonhydrostatic and is_compressible    */
    ud->is_nonhydrostatic =  1;    /* 0: hydrostatic;  1: nonhydrostatic;  -1: transition (see nonhydrostasy()) */
    ud->is_compressible   =  0;    /* 0: psinc;  1: comp;  -1: psinc-comp-transition (see compressibility()) */
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
    ud->xmin =   0.0e+6/ud->h_ref;
    ud->xmax =  10.0e+6/ud->h_ref;
    ud->ymin =   0.0e+0/ud->h_ref;
    ud->ymax =  18.0e+3/ud->h_ref;
    ud->zmin =  -4.0e+6/ud->h_ref;
    ud->zmax =   4.0e+6/ud->h_ref;
    
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
    ud->bottom_theta_bc = BOTTOM_BC_DEFAULT; /* ZERO_ORDER_EXTRAPOL, BOTTOM_BC_DEFAULT */
    
    /* ================================================================================== */
    /* =====  NUMERICS  ================================================================= */
    /* ================================================================================== */
    
    /* time discretization */
    ud->time_integrator       = SI_MIDPT;  /* this code version has only one option */
    ud->advec_time_integrator = STRANG; /* HEUN; EXPL_MIDPT; NO_ADVECTION;  best tested: STRANG; */
    ud->CFL                   = 0.45; /* 0.45; 0.9; 0.8; */
    /* large time step test variant  (N*dt = 20.0, or  dt = 2000 s in the planetary IGW test) */
    ud->dtfixed0              = 10000.0;
    ud->dtfixed               = 10000.0;     
    
    set_time_integrator_parameters(ud);
    
    /* Grid and space discretization */
    ud->inx =  32+1; /* 641; 321; 161; 129; 81; */
    ud->iny = 4*16+1; /* 321; 161;  81;  65; 41;  */
    ud->inz = 4*16+1;
    
    /* explicit predictor step */
    /* Recovery */
    ud->recovery_order = SECOND; /* FIRST, SECOND */ 
    ud->limiter_type_scalars  = NONE; 
    ud->limiter_type_velocity = NONE; 
    /* RUPE; NONE; MONOTONIZED_CENTRAL; MINMOD; VANLEER; SWEBY_MUNZ; SUPERBEE; */
    
    /* parameters for SWEBY_MUNZ limiter family */
    ud->kp = 1.4;
    ud->kz = 1.4; /* Entro abused for velocity in split-step-aligned velocity ! */
    ud->km = 1.4;
    ud->kY = 1.4;
    ud->kZ = 1.4; /* 2.0 */
    
    /* all explicit predictor operations are done on ncache-size data worms to save memory */ 
    ud->ncache = 154; /* 71+4; 304*44; 604*44; (ud->inx+3); (ud->inx+3)*(ud->iny+3);*/
    
    /* linear solver-stuff */
    double tol                            = 1.e-8; // 1.e-8 * (ud->is_compressible == 1 ? 0.01 : 1.0);
    ud->flux_correction_precision         = tol;
    ud->flux_correction_local_precision   = tol;    /* 1.e-05 should be enough */
    ud->second_projection_precision       = tol;
    ud->second_projection_local_precision = tol;  /* 1.e-05 should be enough */
    ud->flux_correction_max_iterations    = 6000;
    ud->second_projection_max_iterations  = 6000;
    ud->initial_projection                = WRONG; /* WRONG;  CORRECT; */
    
    ud->column_preconditioner             = CORRECT; /* WRONG; CORRECT; */
    
    /* numerics parameters */
    ud->eps_Machine = sqrt(DBL_EPSILON);
    
    /* ================================================================================== */
    /* =====  CODE FLOW CONTROL  ======================================================== */
    /* ================================================================================== */
    double day = 24*3600/ud->t_ref;
    
    ud->tout[0] = 0.25*day; /* 3000 */
    ud->tout[1] = 0.50*day; /* 3000 */
    ud->tout[2] = 0.75*day; /* 3000 */
    ud->tout[3] = 1.00*day; /* 3000 */
    ud->tout[4] = -1.0;
    
    ud->stepmax = 10000;
    
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
    
    const int compressible = (ud.is_compressible == 0 ? 0 : 1);
    
    
    /* Skamarock-Klemp 1994 - type test */
    const double xc    =   0.5*(ud.xmax+ud.xmin);
    const double yc    =   9.0e+3/ud.h_ref;
    const double wxz   = 500.0e+3/ud.h_ref;
    const double wy    =   2.0e+3/ud.h_ref;
    const double u0    = 0.0;
    const double v0    = 0.0;
    const double w0    = 0.0;
    const double delth = 0.0 / ud.T_ref;                    /* pot. temp. perturbation amplitude; standard:  0.01 / ud.T_ref */
    
    const int icx = elem->icx;
    const int icy = elem->icy;
    const int icz = elem->icz;
    const int igy = elem->igy;
    
    const int icxn = node->icx;
    const int icyn = node->icy;
    
    int i, j, k, n, nn;
    double x, y, z, xn, yn, zn, r, the, thp;
    double rho, u, v, w, p, rhoY;
    
    States *HySt, *HyStn;
    double *Y, *Yn;
    
    /* store background hydrostatic state in mpv auxiliary struct */
    Hydrostatics_State(mpv, elem, node);
    
    /* Aux space for cell- and node-based column-wise hydrostatics */
    HySt  = States_new(node->icy);
    Y     = (double*)malloc(node->icy*sizeof(double));
    HyStn = States_new(node->icy);
    Yn    = (double*)malloc(node->icy*sizeof(double));
    
    
    for(k = 0; k < icz; k++) {
        
        for(i = 0; i < icx; i++) {
            
            /* set potential temperature stratification in the column */
            for(j = 0; j < icy; j++) {
                x     = elem->x[i];
                y     = elem->y[j];
                z     = elem->z[k];
                r     = sqrt(((x-xc)*(x-xc) + z*z)/(wxz*wxz) + (y-yc)*(y-yc)/(wy*wy));
                r     = MIN_own(r,1.0);
                thp   = (z >= 0.0 ? 1.0 : -1.0) * delth * pow(cos(0.5*PI*r),2.0);
                the   = Thetae(y,z);
                Y[j]  = the + thp;
            }  
            
            for(j = 0; j < icyn; j++) {
                xn    = node->x[i];
                yn    = node->y[j];
                zn    = node->z[k];
                r     = sqrt(((xn-xc)*(xn-xc) + zn*zn)/(wxz*wxz) + (yn-yc)*(yn-yc)/(wy*wy));
                r     = MIN_own(r,1.0);
                thp   = (zn >= 0.0 ? 1.0 : -1.0) * delth * pow(cos(0.5*PI*r),2.0);
                the   = Thetae(yn,zn);
                Yn[j] = the + thp;
            }        

            /* determine hydrostatic pressure distributions column-wise (lateral relation still neglected) */
            Hydrostatics_Column(HySt, HyStn, Y, Yn, elem, node);
            
            /* initialization of field variables */
            for(j = 0; j < icy; j++) {
                
                n  = k*icy*icx + j*icx + i;
                
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
#ifdef FULL_VARIABLES
                Sol->rhoX[BUOY][n] = Sol->rho[n]/Y[j];
#else
                Sol->rhoX[BUOY][n] = Sol->rho[n] * (1.0/Y[j] - mpv->HydroState->S0[j]);
#endif
                
                /* nodal pressure */
                nn   = k*icyn*icxn + j*icxn+i;
                mpv->p2_nodes[nn] = HyStn->p20[j];
            }
        }
    }
    
    States_free(HySt);
    States_free(HyStn);
    free(Y);
    free(Yn);
    
    /* Find geostrophically balanced state corresponding to the Y-distribution */
    // Geostrophic_Initial_State(Sol, mpv, elem, node);  
    
    ud.nonhydrostasy   = nonhydrostasy(0);
    ud.compressibility = compressibility(0);
    
    set_wall_rhoYflux(bdry, Sol, mpv, elem);
    Set_Explicit_Boundary_Data(Sol, elem, OUTPUT_SUBSTEPS);
    
    ConsVars_set(Sol0, Sol, elem->nc);
}

/* ================================================================================== */

double stratification(
                      double y) {
    
    extern User_Data ud;
    
    double Nsq = ud.Nsq_ref * ud.t_ref * ud.t_ref;
    double g   = ud.gravity_strength[1] / ud.Msq;
    
    return( exp(Nsq*y/g) );
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

/* ================================================================================== */

double F(double z,
         double dth) {
    
    const double eta = z/dth;
    
    // return (eta < - 1.0 ? -1 : (eta > 1 ? 1 : sin(PI*eta/2.0) ) );
    return ( tanh(PI*eta/2.0) );
}

/* ================================================================================== */

double fH(double z,
          double dth) {
    
    extern User_Data ud;
    
    // double zmin = ud.zmin;
    // double zmax = ud.zmax;
    // double z0    = 0.25 * ( z > 0.0 ? zmin+3.0*zmax : 3.0*zmin+zmax );
    double z0    = 0.0;
    
    // return (z > 0.0 ? F(z-z0,dth) : -F(z-z0,dth) );
    return ( F(z-z0,dth) );
}

/* ================================================================================== */

double HT(double z) {
    
    extern User_Data ud;
    
    double dth  = 525.0e+3/ud.h_ref; /* [km] */
    double HT0  =   8.0e+3/ud.h_ref; /* [km] */
    double DTh  =  30.0/ud.T_ref;  /* [K] */
    double Th0  = 273.0/ud.T_ref;  /* [K] */
    double Nsqs = 0.0245*0.0245*ud.t_ref*ud.t_ref; /* [s^{-2}] */
    double Nsqt =   0.01*  0.01*ud.t_ref*ud.t_ref; /* [s^{-2}] */
    double dNsq = Nsqs - Nsqt;
    double g    = ud.gravity_strength[1] / ud.Msq;

    double DHT  = 0.5 * g * DTh / Th0 / dNsq;
    
    return ( HT0 - DHT * fH(z,dth) );
}


/* ================================================================================== */

double Zeta(double y,
            double z) {
    
    extern User_Data ud;
    
    double H = 24.0e+3/ud.h_ref;  
    
    return ( sin(0.5*PI*(H - y)/(H-HT(z))) );
}

/* ================================================================================== */

double Thetat(double y,
              double z) {
    
    extern User_Data ud;
    
    double dth  = 525.0e+3/ud.h_ref; /* [km] */
    double DTh  =  30.0/ud.T_ref;  /* [K] */
    double Th0  = 273.0/ud.T_ref;  /* [K] */
    
    double Nsqt = 0.01*0.01*ud.t_ref*ud.t_ref;   /* [s^{-2}] */
    double g    = ud.gravity_strength[1] / ud.Msq;
    
    double kappa = 70.0;    
    double fHv   = fH(z-kappa*y,dth);
    
    return ( Th0 * (1.0 + Nsqt * y / g) - 0.5 * DTh * fHv );
}

/* ================================================================================== */

double Thetas(double y,
              double z) {
    
    extern User_Data ud;
    
    double dth  = 525.0e+3/ud.h_ref; /* [km] */
    double DTh  =  30.0/ud.T_ref;  /* [K] */
    double Th0  = 273.0/ud.T_ref;  /* [K] */
    
    double Nsqs = 0.0245*0.0245*ud.t_ref*ud.t_ref; /* [s^{-2}] */
    double Nsqt =   0.01*  0.01*ud.t_ref*ud.t_ref; /* [s^{-2}] */
    double g    = ud.gravity_strength[1] / ud.Msq;
    
    // double zmin  = ud.zmin;
    // double zmax  = ud.zmax;
    // double z0    = 0.25 * ( z > 0.0 ? zmin+3.0*zmax : 3.0*zmin+zmax );
    double z0    = 0.0;
    double kappa = 70.0;
    
    double zeta = Zeta(y,z);
    // double zeta = 1.0;
    double HTz  = HT(z);
    double H    = 3.0e+3/ud.h_ref;  
    double fHv  = fH(z-kappa*HTz,dth);

    return ( Th0 * (1.0 + Nsqs * y / g) - Th0 * (Nsqs-Nsqt) * HTz * zeta / g - zeta * 0.5*DTh*fH(z-kappa*HTz,dth) * (1.0+(HTz-y)/H)  );
    // best so far: return ( Th0 * (1.0 + Nsqs * y / g - (Nsqs-Nsqt) * HTz / g ) -  0.0 * zeta * 0.5*DTh*fH(z-kappa*HTz,dth) * (1.0+0.5*(HTz-y)/H)  );
    // return ( Th0 * (1.0 + Nsqs * y / g - (Nsqs-Nsqt) * HTz / g ) -  0.0 * zeta * 0.5*DTh*fHv * (1.0+0.5*(HTz-y)/H)  );
    // return ( Th0 * (1.0 + Nsqt * HTz * zeta / g + Nsqs * (y - HTz * zeta) / g ) -  zeta * 0.5*DTh*fHv * (1.0+1.0*(HTz-y)/H)  );
    // return ( Thetat(y,z) + Th0 * Nsqs * (y - HTz) / g  +  1.0 * zeta * 0.5*DTh*fHv * (1.0+0.0*(HTz-y)/H)  );
}

/* ================================================================================== */

double Thetae(double y,
              double z) {
    
    double HTz  = HT(z);
    
    return ( y < HTz ? Thetat(y,z) : Thetas(y,z) );
    // return ( Thetat(y,z) );
}
/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: userdata.c,v $
 Revision 1.2  1998/03/07 09:56:47  nicola
 Added flux computation and multiple pressure variables.
 
 Revision 1.1  1998/03/01 18:50:44  nicola
 Added csr.c and sod1d.c (user data for 3d code)
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
