/*     WeICME - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998 Guido Dhondt                          */

void pardiso_main(double *ad, double *au, double *adb, double *aub, 
         double *sigma,double *b, ITG *icol, ITG *irow, 
	 ITG *neq, ITG *nzs,ITG *symmetryflag,ITG *inputformat,ITG *jq,
	 ITG *nzs3,ITG *nrhs);

void pardiso_factor(double *ad, double *au, double *adb, double *aub, 
                double *sigma,ITG *icol, ITG *irow, 
		ITG *neq, ITG *nzs,ITG *symmetryflag,ITG *inputformat,
		ITG *jq,ITG *nzs3);

void pardiso_solve(double *b,ITG *neq,ITG *symmetryflag,ITG *nrhs);

void pardiso_cleanup(ITG *neq,ITG *symmetryflag);

void FORTRAN(pardiso,(long long *pt,ITG *maxfct,ITG *mnum,ITG *mtype,ITG *phase,
                   ITG *neq,double *aupardiso,ITG *pointers,ITG *irowpardiso,
                   ITG *perm,ITG *nrhs,ITG *iparm,ITG *msglvl,double *b,
                   double *x,ITG *error));

char envMKL[32];
