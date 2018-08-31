/* ========================================================================== */

void euler_forward_non_advective(ConsVars* Sol,
                                 MPV* mpv,
                                 const ConsVars* Sol0,
                                 const ElemSpaceDiscr* elem,
                                 const NodeSpaceDiscr* node,
                                 const double dt,
                                 const enum EXPLICIT_PRESSURE wp)
{
    /* 
     re-evaluates the update computed in the previous time step
     in the second projection in the sense of the beginning 
     Euler forward half time step.
     Straightforwardly, this routine only works with  dt = const.
     Whether it survives for variable time steps (at least with 
     moderate changes from one loop to the next) remains to be 
     analyzed and tested. 
     */        
    double** hplus   = mpv->wplus;
    double*  hcenter = mpv->wcenter;
    
    // euler_backward_gravity(Sol, mpv, elem, dt);
    operator_coefficients_nodes(hplus, hcenter, elem, node, Sol, Sol0, mpv, dt);
    correction_nodes(Sol, elem, node, (const double**)hplus, mpv->p2_nodes, 2.0*dt);
    /* factor 2.0 of dt due to usage of correction_nodes() in second_projection() */
    
    for(int nn=0; nn<node->nc; nn++) {
        mpv->p2_nodes[nn] += mpv->dp2_nodes[nn];
    }
    
    set_ghostnodes_p2(mpv->p2_nodes, node, 2);       
    Set_Explicit_Boundary_Data(Sol, elem);
}
#else

/* ========================================================================== */

void euler_forward_non_advective(ConsVars* Sol,
                                 MPV* mpv,
                                 const ConsVars* Sol0,
                                 const ElemSpaceDiscr* elem,
                                 const NodeSpaceDiscr* node,
                                 const double dt,
                                 const enum EXPLICIT_PRESSURE wp)
{
    /* 
     evaluates Euler forward for the pressure gradient, gravity, 
     and background stratification advection terms based on cell-
     centered data. 
     */
    extern User_Data ud;
    extern Thermodynamic th;
    extern BDRY* bdry;
    
    extern double *W0;
    extern enum Boolean W0_in_use;
    
    double nonhydro = ud.nonhydrostasy;
    
    /*
     double* dp2n = mpv->dp2_nodes;
     */
    assert(W0_in_use == WRONG);
    W0_in_use = CORRECT;
    double* dp2n = W0;
    
    double* p2n  = mpv->p2_nodes;
    
    const double g        = ud.gravity_strength[1];
    const double Msq      = ud.Msq;
    const double Ginv     = th.Gammainv; 
    const double coriolis = ud.coriolis_strength[0];
    const double u0       = ud.wind_speed;
    
    double *div = mpv->rhs;
    double div_max;
    
    memset(dp2n, 0.0, node->nc*sizeof(double));
    memset(div, 0.0, node->nc*sizeof(double));
    
    /* TODO: call to divergence_nodes() might be unnecessary at least after the first
     step, because mpv->dp2_nodes[] should contain the relevant information already
     from the last time step 
     */
    /* last two arguments to  divergence_nodes() tuned to pure divergence calculation */
    int x_periodic, y_periodic, z_periodic;    
    x_periodic = 0;
    y_periodic = 0;
    z_periodic = 0;
    if(ud.bdrytype_min[0] == PERIODIC) x_periodic = 1;
    if(ud.bdrytype_min[1] == PERIODIC) y_periodic = 1;
    if(ud.bdrytype_min[2] == PERIODIC) z_periodic = 1;
    
    double weight = dt * pow(2.0, -(elem->ndim-1));
    div_max = divergence_nodes(div, elem, node, (const ConsVars*)Sol, mpv, bdry, dt, weight);
    catch_periodic_directions(div, node, elem, x_periodic, y_periodic, z_periodic);
    
    switch (elem->ndim) {
        case 1:
        {
            const int icx   = elem->icx;
            const int igx   = elem->igx;
            const double dx = node->dx;
            
            for (int i=igx; i<icx-igx+1; i++) {
                int nc = i;
                int nn0 = i;
                int nn1 = i+1;
                
                double dpdx    = wp*(p2n[nn1]-p2n[nn0])/dx;
                double rhoYovG = Ginv*Sol->rhoY[nc];
                double drhou   = Sol->rhou[nc] - u0*Sol->rho[nc];
                double dpidP   = (th.gm1 / ud.Msq) * \
                0.5 * (pow(Sol->rhoY[nc], th.gamm - 2.0) + pow(Sol->rhoY[nc-1], th.gamm - 2.0));
                /* alternative without need to call pow():  
                 double dpidP   = th.gm1 * mpv->p2_cells[nc] / Sol->rhoY[nc]; 
                 */
                
                Sol->rhou[nc]  = u0*Sol->rho[nc] + dt * (- rhoYovG * dpdx + coriolis * Sol->rhow[nc]);
                Sol->rhow[nc]  = Sol->rhow[nc] - dt * coriolis * drhou;
                Sol->rhoY[nc]  = Sol->rhoY[nc] - dt * div[nc];
                
                dp2n[nn0] -= dt * dpidP * div[nn0];
            }
            ERROR("boundary fix in  euler_forward_non_advective()  not implemented in 1D yet\n");           
        }
            break;
        case 2:
        {
            const int igxe = elem->igx;
            const int icxe = elem->icx;
            const int igye = elem->igy;
            const int icye = elem->icy;
            
            const int icxn = node->icx;
            
            const double dx = node->dx;
            const double dy = node->dy;
            
            for (int j=igye; j<icye-igye+1; j++) {
                int mc = j*icxe;
                int mn = j*icxn;
                double S0p = mpv->HydroState_n->S0[j+1];
                double S0m = mpv->HydroState_n->S0[j];
                
                for (int i=igxe; i<icxe-igxe+1; i++) {
                    int nc  = mc + i;
                    
                    int nn00 = mn  + i;
                    int nn10 = nn00 + icxn;
                    int nn01 = nn00 + 1;
                    int nn11 = nn00 + 1 + icxn;
                    
                    double dpdx    = wp*0.5*(p2n[nn01]-p2n[nn00]+p2n[nn11]-p2n[nn10])/dx;
                    double dpdy    = wp*0.5*(p2n[nn10]-p2n[nn00]+p2n[nn11]-p2n[nn01])/dy;
                    double dSdy    = (S0p-S0m) / dy;
                    
                    double rhoYovG = Ginv*Sol->rhoY[nc];
                    double v       = Sol->rhov[nc]/Sol->rho[nc];
                    double dchi    = Sol->rhoX[BUOY][nc]/Sol->rho[nc];
                    double chi     = Sol->rho[nc]/Sol->rhoY[nc];
                    double dbuoy   = -Sol->rho[nc]*dchi/chi;  /* -dchi/chibar; */
                    double drhou   = Sol->rhou[nc] - u0*Sol->rho[nc];
                    double dpidP   = (th.gm1 / ud.Msq) * \
                    0.25 * (pow(Sol->rhoY[nc], th.gamm - 2.0)      + pow(Sol->rhoY[nc-1], th.gamm - 2.0) + \
                            pow(Sol->rhoY[nc-icxe], th.gamm - 2.0) + pow(Sol->rhoY[nc-icxe-1], th.gamm - 2.0));
                    /* alternative without need to call pow():  
                     double dpidP   = th.gm1 * mpv->p2_cells[nc] / Sol->rhoY[nc]; 
                     */
                    
                    Sol->rhou[nc]  = Sol->rhou[nc] + dt * ( - rhoYovG * dpdx + coriolis * Sol->rhow[nc]);
                    Sol->rhov[nc]  = Sol->rhov[nc] + dt * ( - rhoYovG * dpdy + (g/Msq) * dbuoy) * nonhydro; 
                    Sol->rhow[nc]  = Sol->rhow[nc] - dt * coriolis * drhou;
                    Sol->rhoX[BUOY][nc] += dt * ( - v * dSdy) * Sol->rho[nc];
                    
                    dp2n[nn00] -= dt * dpidP * div[nn00];
                }
            }
        }
            break;
        case 3: 
        {
            const int icx = elem->icx;
            const int icy = elem->icy;
            const int icz = elem->icz;
            const int icxy = icx*icy;
            
            const int igx = elem->igx;
            const int igy = elem->igy;
            const int igz = elem->igz;
            
            const int inx = node->icx;
            const int iny = node->icy;
            
            const double dx = node->dx;
            const double dy = node->dy;
            const double dz = node->dz;
            
            for (int k=igz; k<icz-igz+1; k++) {
                int lc = k*icy*icx;
                int ln = k*iny*inx;
                for (int j=igy; j<icy-igy+1; j++) {
                    int mc = lc + j*icx;
                    int mn = ln + j*inx;
                    double S0p    = mpv->HydroState_n->S0[j+1];
                    double S0m    = mpv->HydroState_n->S0[j];
                    
                    for (int i=igx; i<icx-igx+1; i++) {
                        int nc        = mc + i;
                        
                        int nn000 = mn   + i;
                        int nn010 = nn000 + inx;
                        int nn001 = nn000 + 1;
                        int nn011 = nn000 + 1 + inx;
                        int nn100 = mn   + i + inx*iny;
                        int nn110 = nn100 + inx;
                        int nn101 = nn100 + 1;
                        int nn111 = nn100 + 1 + inx;
                        
                        double dpdx   = wp*0.25*(p2n[nn001]-p2n[nn000]+p2n[nn011]-p2n[nn010]+p2n[nn101]-p2n[nn100]+p2n[nn111]-p2n[nn110])/dx;
                        double dpdy   = wp*0.25*(p2n[nn010]-p2n[nn000]+p2n[nn011]-p2n[nn001]+p2n[nn110]-p2n[nn100]+p2n[nn111]-p2n[nn101])/dy;
                        double dpdz   = wp*0.25*(p2n[nn100]-p2n[nn000]+p2n[nn110]-p2n[nn010]+p2n[nn101]-p2n[nn001]+p2n[nn111]-p2n[nn011])/dz;
                        double dSdy   = (S0p-S0m) / dy;
                        
                        double rhoYovG = Ginv*Sol->rhoY[nc];
                        double v       = Sol->rhov[nc]/Sol->rho[nc];
                        double dchi    = Sol->rhoX[BUOY][nc]/Sol->rho[nc];
                        double chi     = Sol->rho[nc]/Sol->rhoY[nc];
                        double dbuoy   = -Sol->rho[nc]*dchi/chi;  /* -dchi/chibar; */
                        double drhou   = Sol->rhou[nc] - u0*Sol->rho[nc];
                        double dpidP   = (th.gm1 / ud.Msq) * \
                        0.125 * (pow(Sol->rhoY[nc], th.gamm - 2.0)          + pow(Sol->rhoY[nc-1], th.gamm - 2.0) + \
                                 pow(Sol->rhoY[nc-icx], th.gamm - 2.0)      + pow(Sol->rhoY[nc-icx-1], th.gamm - 2.0) + \
                                 pow(Sol->rhoY[nc-icxy], th.gamm - 2.0)     + pow(Sol->rhoY[nc-1-icxy], th.gamm - 2.0) + \
                                 pow(Sol->rhoY[nc-icx-icxy], th.gamm - 2.0) + pow(Sol->rhoY[nc-icx-1-icxy], th.gamm - 2.0));
                        /* alternative without need to call pow():  
                         double dpidP   = th.gm1 * mpv->p2_cells[nc] / Sol->rhoY[nc]; 
                         */
                        
                        Sol->rhou[nc]  = Sol->rhou[nc] + dt * ( - rhoYovG * dpdx + coriolis * Sol->rhow[nc]);
                        Sol->rhov[nc]  = Sol->rhov[nc] + dt * ( - rhoYovG * dpdy + (g/Msq) * dbuoy) * nonhydro; 
                        Sol->rhow[nc]  = Sol->rhow[nc] + dt * ( - rhoYovG * dpdz - coriolis * drhou);
                        Sol->rhoY[nc]  = Sol->rhoY[nc] - dt * div[nc];
                        Sol->rhoX[BUOY][nc] += dt * ( - v * dSdy) * Sol->rho[nc];
                        
                        dp2n[nn000] -= dt * dpidP * div[nn000];
                    }
                }
            }
            ERROR("boundary fix in  euler_forward_non_advective()  not implemented in 3D yet\n");
        }
            break;
            
        default:
            break;
    }
    
    /* last half Euler backward step equals first half Euler forward step */
    if (ud.is_compressible) {
        for (int nn=0; nn<node->nc; nn++) {
#if 1
            mpv->p2_nodes[nn] += dp2n[nn];
#else        
            mpv->p2_nodes[nn] += mpv->dp2_nodes[nn];
#endif
        }
    }
    
    W0_in_use = WRONG;
    
    set_ghostnodes_p2(mpv->p2_nodes, node, 2);       
    Set_Explicit_Boundary_Data(Sol, elem);
    
}

#endif
