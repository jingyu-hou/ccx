/*     WeICME - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998 Guido Dhondt                          */

#ifndef __CCX_SPOOLES_H
#define __CCX_SPOOLES_H

/*
 * seperated from WeICME.h: otherwise everyone would have to include
 * the spooles header files
 */

#include <pthread.h>
#include <misc.h>
#include <FrontMtx.h>
#include <SymbFac.h>
#if USE_MT
#include <MT/spoolesMT.h>
#endif

/* increase this for debugging */
#define DEBUG_LVL	0

struct factorinfo 
{
	ITG size;
	double cpus[11];
	IV *newToOldIV, *oldToNewIV;
	SolveMap *solvemap;
	FrontMtx *frontmtx;
	SubMtxManager *mtxmanager;
	ETree *frontETree;
	ITG nthread;
	FILE *msgFile;

};

void spooles_factor(double *ad, double *au, double *adb, double *aub, 
                    double *sigma, ITG *icol, ITG *irow,
                    ITG *neq, ITG *nzs, ITG *symmetryflag,
                    ITG *inputformat, ITG *nzs3,ITG *singularflag);

void spooles_solve(double *b, ITG *neq);

void spooles_cleanup();

void spooles_factor_rad(double *ad, double *au, double *adb, double *aub, 
                    double *sigma, ITG *icol, ITG *irow,
                    ITG *neq, ITG *nzs, ITG *symmetryflag,
                    ITG *inputformat,ITG *singularflag);

void spooles_solve_rad(double *b, ITG *neq);

void spooles_cleanup_rad();

#endif
