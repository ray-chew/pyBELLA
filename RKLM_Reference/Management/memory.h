/*******************************************************************************
 File:   memory.h
 Author:       
 Date:   
 *******************************************************************************/
#ifndef MEMORY_H
#define MEMORY_H

#include "enumerator.h"

/*------------------------------------------------------------------------------
 global array memory managment: allocates (if allocate) or releases 
 (if !allocate) memory for global arrays
 ------------------------------------------------------------------------------*/
void arrays(int allocate);


/*------------------------------------------------------------------------------
 here a permutation of the x,y,z directions is carried out. the permutation goes
 forward when "direction" = FORWARD, otherwise ...
 ------------------------------------------------------------------------------*/
void rotate2D(ConsVars* Sol, const enum Direction direction);


/*------------------------------------------------------------------------------
 here a permutation of the x,y,z directions is carried out. the permutation goes
 forward when "direction" = FORWARD, otherwise ...
 ------------------------------------------------------------------------------*/
void rotate3D(ConsVars* Sol, const enum Direction direction);


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void flip3D_f( double *f, int ix, int iy, int iz, int n, double *W );


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void flip3D_b( double *f, int ix, int iy, int iz, int n, double *W );


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void flip2D(double *f, int ix, int iy, int n, double *W );

#endif /* MEMORY_H */

/* LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL */

