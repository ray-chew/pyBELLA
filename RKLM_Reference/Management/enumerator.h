/*******************************************************************************
 File:   enumerator.h
 Author: Nicola
 Date:   Tue Feb 24 09:45:43 CET 1998
 *******************************************************************************/
#ifndef ENUMERATOR_H
#define ENUMERATOR_H


/*------------------------------------------------------------------------------
 enumerators
 ------------------------------------------------------------------------------*/
enum Boolean {
    WRONG,
    CORRECT};

enum Constraint {
    VIOLATED,
    SATISFIED};

enum Switch {
    OFF,
    ON};

/* only HDF available in March 2018 code */
enum FileFormat {
    ASCII,
    HDF,
    SILO};

enum No_of_Strang_Sweeps {
    SINGLE_STRANG_SWEEP,
    DOUBLE_STRANG_SWEEP};

enum RecoveryOrder {
    FIRST,
    SECOND};

/* Only the SI_MIDPT option is implemented in the March 2018 code;
   other options were available in oder versions of the code
 */
enum TimeIntegrator {
    OP_SPLIT,
    OP_SPLIT_MD_UPDATE,
    HEUN,
    EXPL_MIDPT,
    RK3_SKAMA,
    RK3_TEST,
    SI_MIDPT};

enum FluxesFrom {
    FLUX_EXTERNAL,
    FLUX_INTERNAL
};

enum MUSCL_ON_OFF {
    WITHOUT_MUSCL,
    WITH_MUSCL
};

enum FORCES_ON_OFF {
    WITHOUT_FORCES,
    WITH_FORCES
};

enum EXPLICIT_PRESSURE {
    WITHOUT_PRESSURE,
    WITH_PRESSURE
};

enum Direction {
    BACKWARD,
    FORWARD};

enum LimiterType {
    NONE,
	MINMOD, 
	VANLEER, 
	VANLEERSmooth,
	SUPERBEE,
	MONOTONIZED_CENTRAL, 
	SWEBY_MUNZ,
	RUPE,
    NO_SLOPE,
    NUMBER_OF_LIMITER};

#endif /* ENUMERATOR_H */




/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: enumerator.h,v $
 Revision 1.1  1998/03/01 18:43:33  nicola
 This is the initial revision of 3d. It comes out after two weeks of work on
 Matthias' version of Rupert's F3D code. It is the starting point for imple
 menting a low Mach number extension.
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
