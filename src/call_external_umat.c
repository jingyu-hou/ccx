/*                                                                              */
/*      WeICME (Wedge Integrated Computational Materials Engineering)           */
/*                 - A 3-dimensional finite element program.                    */   
/*     Developed and maintained by Shenzhen Wedge Central                       */
/*    South Research Institute co., Ltd., Shenzhen, China                       */   
/*     Copy Right 2019-2023.                                                      */

#include<stdio.h>
#include<stdlib.h>

#include"WeICME.h"

typedef double abaqus_real;
typedef int    abaqus_int;

typedef void (*umatptr)(abaqus_real *const,
			abaqus_real *const,
			abaqus_real *const,
			abaqus_real *const,
			abaqus_real *const,
			abaqus_real *const,
			abaqus_real *const,
			abaqus_real *const,
			abaqus_real *const,
			abaqus_real *const,
			const abaqus_real *const,
			const abaqus_real *const,
			const abaqus_real *const,
			const abaqus_real *const,
			const abaqus_real *const,
			const abaqus_real *const,
			const abaqus_real *const,
			const abaqus_real *const,
			const char        *const,
			const abaqus_int  *const,
			const abaqus_int  *const,
			const abaqus_int  *const,
			const abaqus_int  *const,
			const abaqus_real *const,
			const abaqus_int  *const,
			const abaqus_real *const,
			const abaqus_real *const,
			abaqus_real *const,
			const abaqus_real *const,
			const abaqus_real *const,
			const abaqus_real *const,
			const abaqus_int  *const,
			const abaqus_int  *const,
			const abaqus_int  *const,
			const abaqus_int  *const,
			const abaqus_int  *const,
			abaqus_int  *const,
			const int);


void call_external_umat_(abaqus_real *const STRESS,
			 abaqus_real *const STATEV,
			 abaqus_real *const DDSDDE,
			 abaqus_real *const SSE,
			 abaqus_real *const SPD,
			 abaqus_real *const SCD,
			 abaqus_real *const RPL,
			 abaqus_real *const DDSDDT,
			 abaqus_real *const DRPLDE,
			 abaqus_real *const DRPLDT,
			 const abaqus_real *const STRAN,
			 const abaqus_real *const DSTRAN,
			 const abaqus_real *const TIME,
			 const abaqus_real *const DTIME,
			 const abaqus_real *const TEMP,
			 const abaqus_real *const DTEMP,
			 const abaqus_real *const PREDEF,
			 const abaqus_real *const DPRED,
			 const char           *const CMNAME,
			 const abaqus_int  *const NDI,
			 const abaqus_int  *const NSHR,
			 const abaqus_int  *const NTENS,
			 const abaqus_int  *const NSTATV,
			 const abaqus_real *const PROPS,
			 const abaqus_int  *const NPROPS,
			 const abaqus_real *const COORDS,
			 const abaqus_real *const DROT,
			 abaqus_real *const PNEWDT,
			 const abaqus_real *const CELENT,
			 const abaqus_real *const DFGRD0,
			 const abaqus_real *const DFGRD1,
			 const abaqus_int  *const NOEL,
			 const abaqus_int  *const NPT,
			 const abaqus_int  *const LAYER,
			 const abaqus_int  *const KSPT,
			 const abaqus_int  *const KSTEP,
			 abaqus_int  *const KINC,
			 const int size){

#ifdef WeICME_EXTERNAL_BEHAVIOURS_SUPPORT

  const WeICMEExternalBehaviour* uf = WeICME_searchExternalBehaviour(CMNAME);
  if(uf==NULL){
    printf("*ERROR: invalid material\n");
    exit(-1);
  }
  umatptr f = (umatptr) uf->ptr;
  f(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,RPL,
    DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,
    DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
    NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,
    COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,
    NOEL,NPT,LAYER,KSPT,KSTEP,KINC,size);

#endif /* WeICME_EXTERNAL_BEHAVIOURS_SUPPORT */

}
