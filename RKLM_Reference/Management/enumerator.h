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

enum PoissonCentering {
    CELLPOISSON,
    NODEPOISSON};

enum Constraint {
    VIOLATED,
    SATISFIED};

enum Switch {
    OFF,
    ON};

enum FileFormat {
    ASCII,
    HDF,
    SILO};

enum RecoveryOrder {
    FIRST,
    SECOND};

enum InitialCondition {
    ZERO,
    SOD1D,
    ONE_FLAME_BALL};

enum TimeIntegrator {
    OP_SPLIT,
    OP_SPLIT_MD_UPDATE,
    HEUN,
    EXPL_MIDPT,
    RK3_SKAMA,
    RK3_TEST};

enum TimeLevel {
    OLD,
    NEW};

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

enum SolverType {
    JACOBI,
	BICGSTAB, 
    BICGSTAB_PRECON};


enum ToWhichEdge {
    TO_THE_RIGHT,
    TO_THE_LEFT};

#endif /* ENUMERATOR_H */




/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: enumerator.h,v $
 Revision 1.1  1998/03/01 18:43:33  nicola
 This is the initial revision of 3d. It comes out after two weeks of work on
 Matthias' version of Rupert's F3D code. It is the starting point for imple
 menting a low Mach number extension.
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
