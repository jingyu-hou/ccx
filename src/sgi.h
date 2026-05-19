/*     WeICME - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998 Guido Dhondt                          */

#include <scsl_sparse.h>

void sgi_main(double *ad, double *au, double *adb, double *aub, double *sigma,
         double *b, ITG *icol, ITG *irow, 
         ITG *neq, ITG *nzs, ITG token);

void sgi_factor(double *ad, double *au, double *adb, double *aub, 
                double *sigma,ITG *icol, ITG *irow, 
                ITG *neq, ITG *nzs, ITG token);

void sgi_solve(double *b,ITG token);

void sgi_cleanup(ITG token);



