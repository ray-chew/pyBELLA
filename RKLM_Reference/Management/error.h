/*******************************************************************************
 File:   error.h
 Author: Nicola         
 Date:   Tue Dec 30 16:48:42 GMT 1997
 *******************************************************************************/
#ifndef ERROR_H
#define ERROR_H

#include <stdlib.h>
#include <stdio.h>


#define ERROR(message) {\
printf("file %s, line %d: %s\n", __FILE__, __LINE__, message);\
exit(EXIT_FAILURE);\
}




#endif /* ERROR_H */




/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: error.h,v $
 Revision 1.1  1998/03/01 18:43:34  nicola
 This is the initial revision of 3d. It comes out after two weeks of work on
 Matthias' version of Rupert's F3D code. It is the starting point for imple
 menting a low Mach number extension.
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/



