/*******************************************************************************
 File:   grid.c
 Author: Nicola                         Rupert 
 Date:   Thu Feb 26 10:25:23 CET 1998   Feb. 2004
 *******************************************************************************/
#include "Common.h"
#include "kgrid.h"
#include "math_own.h"
#include "error.h"
#include <stdlib.h>
#include <assert.h>
#include <math.h>

Grid* Grid_new(
			   const int inx, 
			   const int iny, 
			   const int inz,
			   const double x0,
			   const double x1,
			   const double y0,
			   const double y1,
			   const double z0,
			   const double z1,
			   const enum BdryType left,
			   const enum BdryType right,
			   const enum BdryType bottom,
			   const enum BdryType top,
			   const enum BdryType back,
			   const enum BdryType front) {
	
	Grid* grid = (Grid*)malloc(sizeof(Grid));
	
	double dx, dy, dz;
	int i;
	
	assert(inx >  1);
	assert(iny >= 1);
	assert(inz >= 1);
	
	grid->ndim = inz > 1 ? 3 : (iny > 1 ? 2 : 1);
	
	grid->inx = inx;
	grid->iny = iny;
	grid->inz = inz;
	
	dx =           (x1 - x0) / (inx - 1)      ;
	dy = iny > 1 ? (y1 - y0) / (iny - 1) : 0.0;
	dz = inz > 1 ? (z1 - z0) / (inz - 1) : 0.0;
	/* GRIDDING
	 dx =           (x1 - x0) / (inx )      ;
	 dy = iny > 1 ? (y1 - y0) / (iny ) : 0.0;
	 dz = inz > 1 ? (z1 - z0) / (inz ) : 0.0;
	 */
	assert(dx >  0.0);
	assert(dy >= 0.0);
	assert(dz >= 0.0); 
	
	grid->dx = dx;
	grid->dy = dy;
	grid->dz = dz;
	
	grid->x0 = x0;
	grid->x1 = x1;
	grid->y0 = y0;
	grid->y1 = y1;
	grid->z0 = z0;
	grid->z1 = z1;
	
	grid->x = (double*)malloc(inx * sizeof(double));
	grid->y = (double*)malloc(iny * sizeof(double));
	grid->z = (double*)malloc(inz * sizeof(double));
	
	for(i = 0; i < inx; i++)
		grid->x[i] = x0 + dx * i;
	
	for(i = 0; i < iny; i++)
		grid->y[i] = y0 + dy * i;
	
	for(i = 0; i < inz; i++)
		grid->z[i] = z0 + dz * i; 
	
	grid->left = left;
	grid->right = right;
	grid->bottom = bottom;
	grid->top = top;
	grid->back = back;
	grid->front = front;
	
	if(iny == 1) {
		grid->bottom = TUNIX;
		grid->top = TUNIX; 
	}
	
	if(inz == 1) {
		grid->back = TUNIX;
		grid->front = TUNIX;
	}
	
	return grid; 
}


void Grid_free(Grid* grid) {   
	free(grid->x);
	free(grid->y);
	free(grid->z);
	free(grid);
}


/* =================================================== */

/* static const double big = 1000000.0;
 */
static const double big = 1.0;



ElemSpaceDiscr* ElemSpaceDiscr_new(Grid* g) {
	
	double x0, y0, z0;
	int i, j, k;
	
	ElemSpaceDiscr* d = (ElemSpaceDiscr*)malloc(sizeof(ElemSpaceDiscr));
	
	assert(g->inx >  1);
	assert(g->iny >= 1);
	assert(g->inz >= 1);
	
	d->ndim   = g->ndim;
    d->normal = big;  /* default value */
		
    d->igx = d->ig[0] = 2;
	d->igy = d->ig[1] = g->iny > 1 ? 2 : 0;
	d->igz = d->ig[2] = g->inz > 1 ? 2 : 0;
    
	d->icx = d->ic[0] =              g->inx - 1 + 2 * d->igx    ;
	d->icy = d->ic[1] = g->iny > 1 ? g->iny - 1 + 2 * d->igy : 1;
	d->icz = d->ic[2] = g->inz > 1 ? g->inz - 1 + 2 * d->igz : 1;
	
	d->nc = d->icx * d->icy * d->icz;
    
    d->stride[0] = 1;
    d->stride[1] = g->iny > 1 ? d->icx : 0;  
    d->stride[2] = g->inz > 1 ? d->icx * d->icy : 0;
	
	d->ifx =              d->icx + 1    ;
	d->ify = g->iny > 1 ? d->icy + 1 : 0;
	d->ifz = g->inz > 1 ? d->icz + 1 : 0;
	
	d->nfx =  d->ifx * d->icy * d->icz;
	d->nfy =  d->icx * d->ify * d->icz;
	d->nfz =  d->icx * d->icy * d->ifz;
	
	d->nf = d->nfx + d->nfy + d->nfz;
	
	d->dx =              g->dx      ;
	d->dy = d->icy > 1 ? g->dy : big;
	d->dz = d->icz > 1 ? g->dz : big;
	d->dxyz[0] = d->dx;
	d->dxyz[1] = d->dy;
	d->dxyz[2] = d->dz;
	
	assert(d->dx > 0.0);
	assert(d->dy > 0.0);
	assert(d->dz > 0.0); 
	
	d->dxmin = MIN_own(d->dx, d->dy);	
	d->dxmin = MIN_own(d->dxmin, d->dz);	
	
	d->x = (double*)malloc(d->icx * sizeof(double));
	d->y = (double*)malloc(d->icy * sizeof(double));
	d->z = (double*)malloc(d->icz * sizeof(double));
	
	x0 =              g->x0 - d->igx * d->dx + 0.5 * d->dx       ;
	y0 = d->icy > 1 ? g->y0 - d->igy * d->dy + 0.5 * d->dy : g->y0;
	z0 = d->icz > 1 ? g->z0 - d->igz * d->dz + 0.5 * d->dz : g->z0;
	
	for(i = 0; i < d->icx; i++) d->x[i] = x0 + d->dx * i; 
	for(j = 0; j < d->icy; j++) d->y[j] = y0 + d->dy * j;
	for(k = 0; k < d->icz; k++) d->z[k] = z0 + d->dz * k;  
	
	d->left = g->left;
	d->right = g->right;
	d->bottom = g->bottom;
	d->top = g->top;
	d->back = g->back;
	d->front = g->front;
	
	d->scale_factor = 1.0;  /* default value */
	
	return d;
}

void ElemSpaceDiscr_free(ElemSpaceDiscr* d) {   
	free(d->x);
	free(d->y);
	free(d->z);
	free(d);
}



NodeSpaceDiscr* NodeSpaceDiscr_new(Grid* g) {
	
	double x0, y0, z0;
	int i, j, k;
	
	NodeSpaceDiscr* d = (NodeSpaceDiscr*)malloc(sizeof(NodeSpaceDiscr));
	
	assert(g->inx >  1);
	assert(g->iny >= 1);
	assert(g->inz >= 1);
	
	d->ndim = g->ndim;
		
    d->igx = d->ig[0] = 2;
	d->igy = d->ig[1] = g->iny > 1 ? 2 : 0;
	d->igz = d->ig[2] = g->inz > 1 ? 2 : 0;
    
	d->icx = d->ic[0] =              g->inx + 2 * d->igx    ;
	d->icy = d->ic[1] = g->iny > 1 ? g->iny + 2 * d->igy : 1;
	d->icz = d->ic[2] = g->inz > 1 ? g->inz + 2 * d->igz : 1;
	
	d->nc = d->icx * d->icy * d->icz;
	
    d->stride[0] = 1;
    d->stride[1] = g->iny > 1 ? d->icx : 0;  
    d->stride[2] = g->inz > 1 ? d->icx * d->icy : 0;
    
	d->ifx = d->icx + 1;
	d->ify = g->iny > 1 ? d->icy + 1 : 0;
	d->ifz = g->inz > 1 ? d->icz + 1 : 0;
	
	d->nfx =  d->ifx * d->icy * d->icz;
	d->nfy =  d->icx * d->ify * d->icz;
	d->nfz =  d->icx * d->icy * d->ifz;
	
	d->nf = d->nfx + d->nfy + d->nfz;
	
	d->dx =              g->dx      ;
	d->dy = d->icy > 1 ? g->dy : big;
	d->dz = d->icz > 1 ? g->dz : big;
    d->dxyz[0] = d->dx;
    d->dxyz[1] = d->dy;
    d->dxyz[2] = d->dz;
	
	assert(d->dx > 0.0);
	assert(d->dy > 0.0);
	assert(d->dz > 0.0); 
    
    d->dxmin = MIN_own(d->dx, d->dy);
    d->dxmin = MIN_own(d->dxmin, d->dz);

	d->x = (double*)malloc(d->icx * sizeof(double));
	d->y = (double*)malloc(d->icy * sizeof(double));
	d->z = (double*)malloc(d->icz * sizeof(double));
	
	x0 =              g->x0 - d->igx * d->dx;
	y0 = d->icy > 1 ? g->y0 - d->igy * d->dy : g->y0;
	z0 = d->icz > 1 ? g->z0 - d->igz * d->dz : g->z0;
	
	for(i = 0; i < d->icx; i++) d->x[i] = x0 + d->dx * i; 
	for(j = 0; j < d->icy; j++) d->y[j] = y0 + d->dy * j;
	for(k = 0; k < d->icz; k++) d->z[k] = z0 + d->dz * k; 
	
	d->left = g->left;
	d->right = g->right;
	d->bottom = g->bottom;
	d->top = g->top;
	d->back = g->back;
	d->front = g->front; 
	
	d->scale_factor = 1.0;  /* default value */
	
	return d;
}

void NodeSpaceDiscr_free(NodeSpaceDiscr* d) {   
	free(d->x);
	free(d->y);
	free(d->z);
	free(d);
}

/*------------------------------------------------------------------------------
 derive surface elements from volume elements 
 ------------------------------------------------------------------------------*/
ElemSpaceDiscr* surface_elems(const ElemSpaceDiscr* elem)
{
    int i, j;
    
    ElemSpaceDiscr* elem_s = (ElemSpaceDiscr*)malloc(sizeof(ElemSpaceDiscr));
        
    elem_s->ndim = elem->ndim-1;
    elem_s->normal = 1;   /* should have normal direction as parameter to fct. */
    
    elem_s->igx = elem_s->ig[0] = elem->igx;
    elem_s->igy = elem_s->ig[1] = elem->igz;
    elem_s->igz = elem_s->ig[2] = 0;
    
    elem_s->icx = elem_s->ic[0] = elem->icx;
    elem_s->icy = elem_s->ic[1] = elem->icz;
    elem_s->icz = elem_s->ic[2] = 1;
    
    elem_s->nc = elem_s->icx * elem_s->icy;
    
    elem_s->stride[0] = 1;
    elem_s->stride[1] = (elem_s->ndim > 1 ? elem_s->icx : 0);  
    elem_s->stride[2] = 0;
        
    elem_s->ifx = elem_s->icx + 1;
    elem_s->ify = (elem_s->ndim > 1 ? elem_s->icy + 1 : 0);
    elem_s->ifz = 0;
    
    elem_s->nfx =  elem_s->ifx * elem_s->icy * elem_s->icz;
    elem_s->nfy =  elem_s->icx * elem_s->ify * elem_s->icz;
    elem_s->nfz =  elem_s->icx * elem_s->icy * elem_s->ifz;
    
    elem_s->nf = elem_s->nfx + elem_s->nfy + elem_s->nfz;
    
    elem_s->dx = elem->dx;
    elem_s->dy = elem->dz;
    elem_s->dz = big;
    elem_s->dxyz[0] = elem_s->dx;
    elem_s->dxyz[1] = elem_s->dy;
    elem_s->dxyz[2] = elem_s->dz;
    
    assert(elem_s->dx > 0.0);
    assert(elem_s->dy > 0.0);
    assert(elem_s->dz > 0.0); 
    
    elem_s->dxmin = MIN_own(elem_s->dx, elem_s->dy);    
    elem_s->dxmin = MIN_own(elem_s->dxmin, elem_s->dz);    
    
    elem_s->x = (double*)malloc(elem_s->icx * sizeof(double));
    elem_s->y = (double*)malloc(elem_s->icy * sizeof(double));
    elem_s->z = (double*)malloc(elem_s->icz * sizeof(double));
        
    for(i = 0; i < elem_s->icx; i++) elem_s->x[i] = elem->x[i]; 
    for(j = 0; j < elem_s->icy; j++) elem_s->y[j] = elem->z[j];
    elem_s->z[0] = elem->y[0];  
    
    elem_s->left = elem->left;
    elem_s->right = elem->right;
    elem_s->bottom = elem->back;
    elem_s->top = elem->front;
    elem_s->back = elem->bottom;
    elem_s->front = elem->top;
    
    elem_s->scale_factor = 1.0;  /* default value */
        
    return elem_s;
}

/*------------------------------------------------------------------------------
 derive surface nodes from volume nodes 
 ------------------------------------------------------------------------------*/
NodeSpaceDiscr* surface_nodes(const NodeSpaceDiscr* node)
{
    int i, j;
    
    NodeSpaceDiscr* node_s = (NodeSpaceDiscr*)malloc(sizeof(NodeSpaceDiscr));
    
    node_s->ndim = node->ndim-1;
    
    node_s->igx = node_s->ig[0] = node->igx;
    node_s->igy = node_s->ig[1] = node->igz;
    node_s->igz = node_s->ig[2] = 0;
    
    node_s->icx = node_s->ic[0] = node->icx;
    node_s->icy = node_s->ic[1] = node->icz;
    node_s->icz = node_s->ic[2] = 1;
    
    node_s->nc = node_s->icx * node_s->icy * node_s->icz;
    
    node_s->stride[0] = 1;
    node_s->stride[1] = (node_s->ndim > 1 ? node_s->icx : 0);  
    node_s->stride[2] = 0;
    
    node_s->ifx = node->ifx;
    node_s->ify = node->ifz;
    node_s->ifz = 0;
    
    node_s->nfx = node_s->ifx * node_s->icy * node_s->icz;
    node_s->nfy = node_s->icx * node_s->ify * node_s->icz;
    node_s->nfz = node_s->icx * node_s->icy * node_s->ifz;
    
    node_s->nf = node_s->nfx + node_s->nfy + node_s->nfz;
    
    node_s->dx = node->dx;
    node_s->dy = node->dz;
    node_s->dz = big;
    node_s->dxyz[0] = node_s->dx;
    node_s->dxyz[1] = node_s->dy;
    node_s->dxyz[2] = node_s->dz;
    
    assert(node_s->dx > 0.0);
    assert(node_s->dy > 0.0);
    assert(node_s->dz > 0.0); 
    
    node_s->dxmin = MIN_own(node_s->dx, node_s->dy);
    node_s->dxmin = MIN_own(node_s->dxmin, node_s->dz);
    
    node_s->x = (double*)malloc(node_s->icx * sizeof(double));
    node_s->y = (double*)malloc(node_s->icy * sizeof(double));
    node_s->z = (double*)malloc(node_s->icz * sizeof(double));
        
    for(i = 0; i < node_s->icx; i++) node_s->x[i] = node->x[i]; 
    for(j = 0; j < node_s->icy; j++) node_s->y[j] = node->z[j];
    node_s->z[0] = node->y[0]; 
    
    node_s->left   = node->left;
    node_s->right  = node->right;
    node_s->bottom = node->back;
    node_s->top    = node->front;
    node_s->back   = node->bottom;
    node_s->front  = node->top; 
    
    node_s->scale_factor = 1.0;  /* default value */
    
    return node_s;
}

/*------------------------------------------------------------------------------
 extrude surface date to volume data / cell-based
 ------------------------------------------------------------------------------*/

void extrude_cells(double* q_vol, 
                   const double* q_surf, 
                   const ElemSpaceDiscr* elem, 
                   const ElemSpaceDiscr* elem_surf)
{
    /* lift the surface quantity  q_surf  to the volume quantity  q_vol
     by zeroeth order extrapolation
     */
    switch (elem_surf->normal) {
            
        case 0:
            assert(elem->icy == elem_surf->icx);
            assert(elem->icz == elem_surf->icy);
            for (int j=0; j<elem->icy; j++) {
                const int lcv = j*elem->stride[1];
                const int lcs = j*elem_surf->stride[0];
                for (int k=0; k<elem->icz; k++) {
                    const int mcv = lcv + k*elem->stride[2];
                    const int mcs = lcs + k*elem_surf->stride[1];
                    for (int i=0; i<elem->icx; i++) {
                        const int ncv = mcv + i*elem->stride[0];
                        q_vol[ncv] = q_surf[mcs];
                    }
                }
            }
            break;
            
        case 1:
            assert(elem->icx == elem_surf->icx);
            assert(elem->icz == elem_surf->icy);
            for (int k=0; k<elem->icz; k++) {
                const int lcv = k*elem->stride[2];
                const int lcs = k*elem_surf->stride[1];
                for (int i=0; i<elem->icx; i++) {
                    const int mcv = lcv + i*elem->stride[0];
                    const int mcs = lcs + i*elem_surf->stride[0];
                    for (int j=0; j<elem->icy; j++) {
                        const int ncv = mcv + j*elem->stride[1];
                        q_vol[ncv] = q_surf[mcs];
                    }
                }
            }
            break;
            
        case 2:
            assert(elem->icx == elem_surf->icx);
            assert(elem->icy == elem_surf->icy);
            for (int j=0; j<elem->icy; j++) {
                const int lcv = j*elem->stride[1];
                const int lcs = j*elem_surf->stride[1];
                for (int i=0; i<elem->icx; i++) {
                    const int mcv = lcv + i*elem->stride[0];
                    const int mcs = lcs + i*elem_surf->stride[0];
                    for (int k=0; k<elem->icz; k++) {
                        const int ncv = mcv + k*elem->stride[2];
                        q_vol[ncv] = q_surf[mcs];
                    }
                }
            }
            break;
            
        default:
            ERROR("\nNo plausible data extrusion direction provided.\n");
            break;
    }
}

/*------------------------------------------------------------------------------
 extrude surface date to volume data / node-based
 ------------------------------------------------------------------------------*/

void extrude_nodes(double* q_vol, 
                   const double* q_surf, 
                   const NodeSpaceDiscr* node, 
                   const NodeSpaceDiscr* node_surf)
{
    /* lift the surface quantity  q_surf  to the volume quantity  q_vol
     by zeroeth order extrapolation
     */
    switch (node_surf->normal) {
            
        case 0:
            assert(node->icy == node_surf->icx);
            assert(node->icz == node_surf->icy);
            for (int j=0; j<node->icy; j++) {
                const int lcv = j*node->stride[1];
                const int lcs = j*node_surf->stride[0];
                for (int k=0; k<node->icz; k++) {
                    const int mcv = lcv + k*node->stride[2];
                    const int mcs = lcs + k*node_surf->stride[1];
                    for (int i=0; i<node->icx; i++) {
                        const int ncv = mcv + i*node->stride[0];
                        q_vol[ncv] = q_surf[mcs];
                    }
                }
            }
            break;
            
        case 1:
            assert(node->icx == node_surf->icx);
            assert(node->icz == node_surf->icy);
            for (int k=0; k<node->icz; k++) {
                const int lcv = k*node->stride[2];
                const int lcs = k*node_surf->stride[1];
                for (int i=0; i<node->icx; i++) {
                    const int mcv = lcv + i*node->stride[0];
                    const int mcs = lcs + i*node_surf->stride[0];
                    for (int j=0; j<node->icy; j++) {
                        const int ncv = mcv + j*node->stride[1];
                        q_vol[ncv] = q_surf[mcs];
                    }
                }
            }
            break;
            
        case 2:
            assert(node->icx == node_surf->icx);
            assert(node->icy == node_surf->icy);
            for (int j=0; j<node->icy; j++) {
                const int lcv = j*node->stride[1];
                const int lcs = j*node_surf->stride[1];
                for (int i=0; i<node->icx; i++) {
                    const int mcv = lcv + i*node->stride[0];
                    const int mcs = lcs + i*node_surf->stride[0];
                    for (int k=0; k<node->icz; k++) {
                        const int ncv = mcv + k*node->stride[2];
                        q_vol[ncv] = q_surf[mcs];
                    }
                }
            }
            break;
            
        default:
            ERROR("\nNo plausible data extrusion direction provided.\n");
            break;
    }
}

