/*******************************************************************************
 File:   memory.h
 Author:       
 Date:   
 *******************************************************************************/
#ifndef MEMORY_H
#define MEMORY_H

#include "ProjectionType.h"
#include "enumerator.h"

void update(
			ConsVars* sol, 
			const ConsVars* flux[3], 
			const VectorField* buoy,
			const ElemSpaceDiscr* elem, 
			const double dt);


/*------------------------------------------------------------------------------
 global array memory managment: allocates (if allocate) or releases 
 (if !allocate) memory for global arrays
 ------------------------------------------------------------------------------*/
void arrays(int allocate);


/*------------------------------------------------------------------------------
 here a permutation of the x,y,z directions is carried out. the permutation goes
 forward when "direction" = FORWARD, otherwise ...
 ------------------------------------------------------------------------------*/
void rotate2D(ConsVars* Sol, double *rhs, double *Yinvbg, const enum Direction direction);


/*------------------------------------------------------------------------------
 here a permutation of the x,y,z directions is carried out. the permutation goes
 forward when "direction" = FORWARD, otherwise ...
 ------------------------------------------------------------------------------*/
void rotate3D(ConsVars* Sol, double *rhs, double *Yinvbg, const enum Direction direction);


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void flip3D_f( double *f, int ix, int iy, int iz, int n, double *W );


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void flip3D_b( double *f, int ix, int iy, int iz, int n, double *W );


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void flip2D(double *f, int ix, int iy, int n, double *W );


#ifdef SOLVER_1_HYPRE

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void map2D_to_michaels_memory_cells(
                                    double *U,
                                    const ElemSpaceDiscr *elem,
                                    double *W0);

void map2D_to_michaels_memory_faces_x(
                                      double *U,
                                      const ElemSpaceDiscr *elem,
                                      double *W0,
                                      double *W1);

void map2D_to_michaels_memory_faces_y(
                                      double *U,
                                      const ElemSpaceDiscr *elem,
                                      double *W0,
                                      double *W1);

void map2D_to_ruperts_memory_cells(
                                   double *U,
                                   const ElemSpaceDiscr *elem,
                                   double *W0);

void map2D_to_ruperts_memory_faces_x(
                                     double *U,
                                     const ElemSpaceDiscr *elem,
                                     double *W0);

void map2D_to_ruperts_memory_interface_coeffs_x(
                                                double *U,
                                                const ElemSpaceDiscr *elem,
                                                double *W0);
void map2D_to_ruperts_memory_faces_y(
                                     double *U,
                                     const ElemSpaceDiscr *elem,
                                     double *W0);

void map2D_to_ruperts_memory_interface_coeffs_y(
                                                double *U,
                                                const ElemSpaceDiscr *elem,
                                                double *W0);
void map2D_to_michaels_memory_nodes(
                                    double *rhs,
                                    const NodeSpaceDiscr *node,
                                    double *W0);

void map2D_to_ruperts_memory_nodes(
                                   double *U,
                                   const NodeSpaceDiscr *node,
                                   double *W0);

void map2D_to_michaels_memory_nodefaces_x(
                                          double *U,
                                          const NodeSpaceDiscr *node,
                                          double *W0,
                                          double *W1);

void map2D_to_michaels_memory_nodefaces_y(
                                          double *U,
                                          const NodeSpaceDiscr *node,
                                          double *W0,
                                          double *W1);

void map2D_to_michaels_memory_opcoeffs_cells(
                                       double *h[3],
                                       double *hplus[3],
                                       const ElemSpaceDiscr *elem);

void map2D_to_michaels_memory_nodefaces(
                                        double *h[3],
                                        double *hplus[3],
                                        const ElemSpaceDiscr *elem,
                                        const NodeSpaceDiscr *node);


#endif /* SOLVER_1_HYPRE */


#endif /* MEMORY_H */




/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: memory.h,v $
 Revision 1.1  1998/03/01 18:43:35  nicola
 This is the initial revision of 3d. It comes out after two weeks of work on
 Matthias' version of Rupert's F3D code. It is the starting point for imple
 menting a low Mach number extension.
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/


