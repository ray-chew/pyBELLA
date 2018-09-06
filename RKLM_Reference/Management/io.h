/*******************************************************************************
 File:   io.h
 Author: Nicola
 Date:   Wed Feb 18 08:42:53 CET 1998
 
 Notice: we should distinguish between read/write functions for the grid and
 for the solution! 
 *******************************************************************************/
#ifndef IO_H
#define IO_H

#include <stdio.h> 
/* #include "space_discretization.h" */ 
#include "kgrid.h"
#include "enumerator.h"
#include "variable.h"


/*------------------------------------------------------------------------------
 putout: write out grid and solution
 ------------------------------------------------------------------------------*/
void putout(ConsVars* Sol, 
            char* dir_name, 
            char* field_name,
            const ElemSpaceDiscr *elem,
            const NodeSpaceDiscr *node,
            const int writeout) ;

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void WriteHDF(
			  FILE* pfile, 
			  int rows, 
			  int cols, 
			  int layers, 
			  int ndim,
			  double* Data, 
			  char* file_name, 
			  char* var_name);


#ifdef ONE_POINT_TIME_SERIES

void initialize_time_series(void);
void store_time_series_entry(const ConsVars *Sol,
                             const ElemSpaceDiscr *elem,
                             const int step);
void close_time_series(void);

#endif /* ONE_POINT_TIME_SERIES */


#endif /* IO_H */




/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: io.h,v $
 Revision 1.2  1998/03/07 09:56:46  nicola
 Added flux computation and multiple pressure variables.
 
 Revision 1.1  1998/03/01 18:43:34  nicola
 This is the initial revision of 3d. It comes out after two weeks of work on
 Matthias' version of Rupert's F3D code. It is the starting point for imple
 menting a low Mach number extension.
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
