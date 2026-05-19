/*     WeICME - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998 Guido Dhondt                          */

   
void matrixstorage(double *ad,double **aup,double *adb,double *aub,
                   double *sigma,ITG *icol,ITG **irowp,ITG *neq,ITG *nzs,
                   ITG *ntrans,ITG *inotr,double *trab,double *co,ITG *nk,
                   ITG *nactdof,char *jobnamec,ITG *mi,ITG *ipkon,
                   char *lakon,ITG *kon,ITG *ne,ITG *mei,ITG *nboun,
                   ITG *nmpc,double *cs,ITG *mcs);
