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
void putout(
			ConsVars* Sol, 
			const double t, 
			const double tout, 
			const int step, 
			const int SplitStep,
			char* dir_name, 
			char* field_name, 
            const int writeout);


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void WriteTimeHistories(const ConsVars* Sol,
                               const ElemSpaceDiscr* elem,
                               const double time,
                               const int step,
                               const int first_running_last);

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


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void ElemSpaceDiscrWriteASCII(
							  double* var, 
							  const ElemSpaceDiscr* elem,
							  const char* filename,
							  const char* varname);


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void NodeSpaceDiscrWriteASCII(
							  double* var, 
							  const NodeSpaceDiscr* node,
							  const char* filename,
							  const char* varname);




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
