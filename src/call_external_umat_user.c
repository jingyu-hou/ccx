/*                                                                              */
/*      WeICME (Wedge Integrated Computational Materials Engineering)           */
/*                 - A 3-dimensional finite element program.                    */   
/*     Developed and maintained by Shenzhen Wedge Central                       */
/*    South Research Institute co., Ltd., Shenzhen, China                       */   
/*     Copy Right 2019-2023.                                                      */

#include<stdio.h>
#include<stdlib.h>

#include"WeICME.h"

typedef double ccxreal;
typedef int    ccxint;

typedef void (*WeICMEptr)(const char * const,
			    const ccxint* const,
			    const ccxint* const,
			    const ccxint* const,
			    const ccxreal* const,
			    const ccxreal* const,
			    const ccxreal* const,
			    const ccxreal* const,
			    const ccxreal* const,
			    const ccxreal* const,
			    const ccxreal* const,
			    const ccxreal* const,
			    const ccxint* const,
			    const ccxreal* const,
			    const ccxreal* const,
			    const ccxreal* const,
			    const ccxreal* const,
			    const ccxint* const,
			    const ccxint* const,
			    const ccxint* const,
			    const ccxint* const,
			    const ccxreal* const,
			    ccxreal* const,
			    ccxreal* const,
			    ccxreal* const,
			    const ccxint* const,
			    const ccxreal* const,
			    const ccxreal* const,
			    ccxreal* const,
			    const ccxint* const,
			    const int);

void call_external_umat_user_(const char * const amat,
			      const ccxint* const iel,
			      const ccxint* const iint,
			      const ccxint* const NPROPS,
			      const ccxreal* const MPROPS,
			      const ccxreal* const STRAN1,
			      const ccxreal* const STRAN0,
			      const ccxreal* const beta,
			      const ccxreal* const F0,
			      const ccxreal* const voj,
			      const ccxreal* const F1,
			      const ccxreal* const vj,
			      const ccxint* const ithermal,
			      const ccxreal* const TEMP1,
			      const ccxreal* const DTIME,
			      const ccxreal* const time,
			      const ccxreal* const ttime,
			      const ccxint* const icmd,
			      const ccxint* const ielas,
			      const ccxint* const mi,
			      const ccxint* const NSTATV,
			      const ccxreal* const STATEV0,
			      ccxreal* const STATEV1,
			      ccxreal* const STRESS,
			      ccxreal* const DDSDDE,
			      const ccxint* const iorien,
			      const ccxreal* const pgauss,
			      const ccxreal* const orab,
			      ccxreal* const PNEWDT,
			      const ccxint* const ipkon,
			      const int size){
#ifdef WeICME_EXTERNAL_BEHAVIOURS_SUPPORT
  const WeICMEExternalBehaviour* uf = WeICME_searchExternalBehaviour(amat);
  if(uf==NULL){
    printf("*ERROR: invalid material\n");
    exit(-1);
  }
  WeICMEptr f = (WeICMEptr) uf->ptr;
  f(amat,iel,iint,NPROPS,MPROPS,STRAN1,STRAN0,beta,
    F0, voj,F1,vj,ithermal,TEMP1,DTIME,time,ttime,
    icmd, ielas,mi,NSTATV,STATEV0,STATEV1,
    STRESS,DDSDDE,iorien,pgauss,orab,PNEWDT,
    ipkon,size);
#endif /* WeICME_EXTERNAL_BEHAVIOURS_SUPPORT */
}
