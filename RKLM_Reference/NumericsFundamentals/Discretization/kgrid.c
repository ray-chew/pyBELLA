/*******************************************************************************
 File:   grid.c
 Author: Nicola                         Rupert 
 Date:   Thu Feb 26 10:25:23 CET 1998   Feb. 2004
 *******************************************************************************/
#include "Common.h"
#include "kgrid.h"
#include "math_own.h"
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
	
	d->ndim = g->ndim;
		
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
	
	d->nci = (g->inx - 1) 
	* (g->iny > 1 ? g->iny - 1 : 1) 
	* (g->inz > 1 ? g->inz - 1 : 1);
	
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


