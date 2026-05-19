/*     WeICME - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998 Guido Dhondt                          */

void tau(double *ad, double **aup, double *adb, double *aubp, double *sigma,
         double *b, ITG *icol, ITG **irowp, 
         ITG *neq, ITG *nzs);

void tau_factor(double *ad, double **aup, double *adb, double *aub, 
                double *sigma,ITG *icol, ITG **irowp, 
                ITG *neq, ITG *nzs);

void tau_solve(double *b,ITG *neq);

void tau_cleanup();



