/*******************************************************************************
 File:   kgrid.h
 Author: Nicola                        Rupert
 Date:   Thu Feb 26 09:57:21 CET 1998  Feb. 2004
 *******************************************************************************/
#ifndef KGRID_H
#define KGRID_H

#include "enum_bdry.h"


/*------------------------------------------------------------------------------
 Grid (regular grid)
 ------------------------------------------------------------------------------*/
typedef struct {
	
	int ndim;    /* space dimension */
	
	int inx;     /* number of nodes in x */
	int iny;
	int inz;
	
	double dx;   /* spacing */
	double dy;
	double dz;
	double dxmin;
	
	double x0;
	double x1;  
	double y0;
	double y1;
	double z0;
	double z1;
	
	double* x;   /* coordinates of the nodes */
	double* y;
	double* z;
	
	enum BdryType left;
	enum BdryType right;
	enum BdryType bottom;
	enum BdryType top;
	enum BdryType back;
	enum BdryType front;
	
} Grid;


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
Grid* Grid_new(
			   const int nnx, 
			   const int nny, 
			   const int nnz,
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
			   const enum BdryType front);


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void Grid_free(Grid* grid);

/*------------------------------------------------------------------------------
 ElemSpaceDiscr (element centered discretization of a regular grid)
 ------------------------------------------------------------------------------*/
typedef struct {
	
	int ndim;    /* number of space dimensions */
    int normal;  /* big: volume grid; i: surface grid normal to direction i */
	
	int nc;      /* total number of cells (internal + ghost) */
	
	int nf;      /* total number of faces (internal + ghost) */
	
	int nfx;     /* total number of faces (internal + ghost) normal to x */
	int nfy;
	int nfz;
	
	int igx;     /* number of ghosts cells and faces in x per each side */
	int igy; 
	int igz; 
    int ig[3];
	
	int icx;     /* number of cells (internal + ghost) in x */ 
	int icy; 
	int icz;
	int ic[3]; 
	
	int ifx;     /* number of faces (internal + ghost) normal to x in x */
	int ify;
	int ifz;
	
    int stride[3];  /* strides on the grid in the three directions */
	
    double dx;   /* spacing in x */
	double dy;
	double dz;
	double dxmin;
	double dxyz[3];
	
	double scale_factor;
	
	double* x;   /* x coordinate of the cell centers */
	double* y;
	double* z;
	
	enum BdryType left;
	enum BdryType right;
	enum BdryType bottom;
	enum BdryType top;
	enum BdryType back;
	enum BdryType front;
	
} ElemSpaceDiscr;

ElemSpaceDiscr* ElemSpaceDiscr_new(Grid* g);

void ElemSpaceDiscr_free(ElemSpaceDiscr* d);


/*------------------------------------------------------------------------------
 NodeSpaceDiscr (node centered discretization of a regular grid)
 ------------------------------------------------------------------------------*/
typedef struct {
	
	int ndim;    /* space dimension */
    int normal;  /* 0: volume grid; i: surface grid normal to direction i */
	
	int nc;      /* total number of cells (internal + ghost) */
	
	int nf;      /* total number of faces (internal + ghost) */
	
	int nfx;     /* total number of faces (internal + ghost) normal to x */
	int nfy;
	int nfz;
	
	int igx;     /* number of ghosts cells and faces in x per each side */
	int igy; 
	int igz; 
    int ig[3];

	int icx;     /* number of cells (internal + ghost) in x */ 
	int icy; 
	int icz; 
    int ic[3];

	int ifx;     /* number of faces (internal + ghost) normal to x in x */
	int ify;
	int ifz;
	
    int stride[3];  /* strides on the grid in the three directions */

	double dx;   /* spacing in x */
	double dy;
	double dz;
    double dxmin;
    double dxyz[3];
    
	double scale_factor;
	
	double* x;   /* x coordinate of the cell centers */
	double* y;
	double* z;
	
	enum BdryType left;
	enum BdryType right;
	enum BdryType bottom;
	enum BdryType top;
	enum BdryType back;
	enum BdryType front;
	
} NodeSpaceDiscr;

NodeSpaceDiscr* NodeSpaceDiscr_new(Grid* g);

void NodeSpaceDiscr_free(NodeSpaceDiscr* d);


/*------------------------------------------------------------------------------
 derive surface grids from volume grids
 ------------------------------------------------------------------------------*/
ElemSpaceDiscr* surface_elems(const ElemSpaceDiscr* elem);
NodeSpaceDiscr* surface_nodes(const NodeSpaceDiscr* node);

/*------------------------------------------------------------------------------
 extrude surface data to volume data / cell-based
 ------------------------------------------------------------------------------*/
void extrude_cells(double* p2_aux, 
                   const double* p2_surf, 
                   const ElemSpaceDiscr* elem, 
                   const ElemSpaceDiscr* elem_surf);

/*------------------------------------------------------------------------------
 extrude surface data to volume data / node-based
 ------------------------------------------------------------------------------*/
void extrude_nodes(double* p2_aux, 
                   const double* p2_surf, 
                   const NodeSpaceDiscr* node, 
                   const NodeSpaceDiscr* node_surf);

#endif /* KGRID_H */



