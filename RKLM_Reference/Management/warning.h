/*******************************************************************************
 File:   warning.h
 Author: Nicola         
 Date:   Tue Feb 17 16:04:45 WET 1998
 *******************************************************************************/
#ifndef WARNING_H
#define WARNING_H


#include <stdlib.h>


#define WARNING(message)\
printf("file %s, line %d: %s\n", __FILE__, __LINE__, message)




#endif /* WARNING_H */




/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: warning.h,v $
 Revision 1.1  1998/03/01 18:43:37  nicola
 This is the initial revision of 3d. It comes out after two weeks of work on
 Matthias' version of Rupert's F3D code. It is the starting point for imple
 menting a low Mach number extension.
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/



