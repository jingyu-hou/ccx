/*                                                                              */
/*      WeICME (Wedge Integrated Computational Materials Engineering)           */
/*                 - A 3-dimensional finite element program.                    */
/*     Developed and maintained by Shenzhen Wedge Central                       */
/*    South Research Institute co., Ltd., Shenzhen, China                       */
/*     Copy Right 2019-2023.                                                    */

#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <pthread.h>
#include "WeICME.h"

static char *lakon1,*matname1,*sideload1;

static ITG *kon1,*ipkon1,*ne1,*nelcon1,*nrhcon1,*nalcon1,*ielmat1,*ielorien1,
    *norien1,*ntmat1_,*ithermal1,*iprestr1,*iperturb1,*iout1,*nmethod1,
    *nplicon1,*nplkcon1,*npmat1_,*mi1,*ielas1,*icmd1,*ncmat1_,*nstate1_,
    *istep1,*iinc1,calcul_fn1,calcul_qa1,calcul_cauchy1,*nener1,ikin1,
    *nal=NULL,*ipompc1,*nodempc1,*nmpc1,*ncocon1,*ikmpc1,*ilmpc1,
    num_cpus,mt1,*nk1,*ne01,*nshcon1,*nelemload1,*nload1,*mortar1,
    *ielprop1,*kscale1,*iponoel1,*inoel1,*network1,*ipobody1,*ibody1,
    *neapar=NULL,*nebpar=NULL,*nmpcon1,*nmpmat1_,*nphase1,*phase_inf1;

static double *co1,*v1,*stx1,*elcon1,*rhcon1,*alcon1,*alzero1,*orab1,*t01,*t11,
    *prestr1,*eme1,*fn1=NULL,*qa1=NULL,*vold1,*veold1,*dtime1,*time1,
    *ttime1,*plicon1,*plkcon1,*xstateini1,*xstiff1,*xstate1,*stiini1,
    *vini1,*ener1,*eei1,*enerini1,*springarea1,*reltime1,*coefmpc1,
    *cocon1,*qfx1,*thicke1,*emeini1,*shcon1,*xload1,*prop1,
    *xloadold1,*pslavsurf1,*pmastsurf1,*clearini1,*xbody1,*mpcon1,
    *pphase1, *cphase1,*phaseother1,*rdpcon1,*gscon1,*co2 ;


void results(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,ITG *ne,
       double *v,double *stn,ITG *inum,double *stx,double *elcon,ITG *nelcon,
       double *rhcon,ITG *nrhcon,double *alcon,ITG *nalcon,double *alzero,
       ITG *ielmat,ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
       double *t0,
       double *t1,ITG *ithermal,double *prestr,ITG *iprestr,char *filab,
       double *eme,double *emn,
       double *een,ITG *iperturb,double *f,double *fn,ITG *nactdof,ITG *iout,
       double *qa,double *vold,double *b,ITG *nodeboun,ITG *ndirboun,
       double *xboun,ITG *nboun,ITG *ipompc,ITG *nodempc,double *coefmpc,
       char *labmpc,ITG *nmpc,ITG *nmethod,double *cam,ITG *neq,double *veold,
       double *accold,double *bet,double *gam,double *dtime,double *time,
       double *ttime,double *plicon,ITG *nplicon,double *plkcon,
       ITG *nplkcon,double *xstateini,double *xstiff,double *xstate,ITG *npmat_,
       double *epn,char *matname,ITG *mi,ITG *ielas,ITG *icmd,ITG *ncmat_,
       ITG *nstate_,
       double *stiini,double *vini,ITG *ikboun,ITG *ilboun,double *ener,
       double *enern,double *emeini,double *xstaten,double *eei,double *enerini,
       double *cocon,ITG *ncocon,char *set,ITG *nset,ITG *istartset,
       ITG *iendset,
       ITG *ialset,ITG *nprint,char *prlab,char *prset,double *qfx,double *qfn,
       double *trab,
       ITG *inotr,ITG *ntrans,double *fmpc,ITG *nelemload,ITG *nload,
       ITG *ikmpc,ITG *ilmpc,
       ITG *istep,ITG *iinc,double *springarea,double *reltime, ITG *ne0,
       double *thicke,
       double *shcon,ITG *nshcon,char *sideload,double *xload,
       double *xloadold,ITG *icfd,ITG *inomat,double *pslavsurf,
       double *pmastsurf,ITG *mortar,ITG *islavact,double *cdn,
       ITG *islavnode,ITG *nslavnode,ITG *ntie,double *clearini,
       ITG *islavsurf,ITG *ielprop,double *prop,double *energyini,
       double *energy,ITG *kscale,ITG *iponoel,ITG *inoel,ITG *nener,
       char *orname,ITG *network,ITG *ipobody,double *xbody,ITG *ibody,
       double *mpcon,ITG *nmpcon,ITG *nmpmat_,double *pphase,
       double *cphase,double *phaseother,ITG *nphase,ITG *phase_inf,double *rdpcon,
       double *gscon){

    ITG intpointvarm,calcul_fn,calcul_f,calcul_qa,calcul_cauchy,ikin,
        intpointvart,mt=mi[1]+1,i,j,k,list1,*ilist1=NULL,nea,neb;

    /*

     calculating integration point values (strains, stresses,
     heat fluxes, material tangent matrices and nodal forces)

     storing the nodal and integration point results in the
     .dat file

     iout=-2: v is assumed to be known and is used to
              calculate strains, stresses..., no result output
              corresponds to iout=-1 with in addition the
              calculation of the internal energy density
     iout=-1: v is assumed to be known and is used to
              calculate strains, stresses..., no result output;
              is used to take changes in SPC's and MPC's at the
              start of a new increment or iteration into account
     iout=0: v is calculated from the system solution
             and strains, stresses.. are calculated, no result output
     iout=1:  v is calculated from the system solution and strains,
              stresses.. are calculated, requested results output
     iout=2: v is assumed to be known and is used to
             calculate strains, stresses..., requested results output */

    /* variables for multithreading procedure */

    ITG sys_cpus,*ithread=NULL,force,ndim,nfield,iorienloc;
    char *env,*envloc,*envsys,*cflag;

    double *pwork,press,pressx[mi[0]**ne],pressn[*nk];
    FILE *pFile;

    num_cpus = 0;
    sys_cpus=0;

    /* explicit user declaration prevails */

    envsys=getenv("NUMBER_OF_CPUS");
    if(envsys){
	sys_cpus=atoi(envsys);
	if(sys_cpus<0) sys_cpus=0;
    }

    /* automatic detection of available number of processors */

    if(sys_cpus==0){
	sys_cpus = getSystemCPUs();
	if(sys_cpus<1) sys_cpus=1;
    }

    /* local declaration prevails, if strictly positive */

    envloc = getenv("CCX_NPROC_RESULTS");
    if(envloc){
	num_cpus=atoi(envloc);
	if(num_cpus<0){
	    num_cpus=0;
	}else if(num_cpus>sys_cpus){
	    num_cpus=sys_cpus;
	}

    }

    /* else global declaration, if any, applies */

    env = getenv("OMP_NUM_THREADS");
    if(num_cpus==0){
	if (env)
	    num_cpus = atoi(env);
	if (num_cpus < 1) {
	    num_cpus=1;
	}else if(num_cpus>sys_cpus){
	    num_cpus=sys_cpus;
	}
    }

// next line is to be inserted in a similar way for all other paralell parts

    if(*ne<num_cpus) num_cpus=*ne;

    pthread_t tid[num_cpus];

    /* 1. nodewise storage of the primary variables
       2. determination which derived variables have to be calculated */
//  intpointvart=222;printf("intpointvart=%d, mi(2)=%d\n",intpointvart,mi[1]);
//  for(i=0;i<*nk;i++){printf("check1 %d %d %d %d %d %d\n",i,nactdof[5*i],nactdof[5*i+1],nactdof[5*i+2],nactdof[5*i+3],nactdof[5*i+4]);}
    FORTRAN(resultsini,(nk,v,ithermal,filab,iperturb,f,fn,
       nactdof,iout,qa,vold,b,nodeboun,ndirboun,
       xboun,nboun,ipompc,nodempc,coefmpc,labmpc,nmpc,nmethod,cam,neq,
       veold,accold,bet,gam,dtime,mi,vini,nprint,prlab,
       &intpointvarm,&calcul_fn,&calcul_f,&calcul_qa,&calcul_cauchy,
       &ikin,&intpointvart));
//for(i=0;i<*nk;i++){printf("check2 %d %d %d %d %d %d\n",i,nactdof[5*i],nactdof[5*i+1],nactdof[5*i+2],nactdof[5*i+3],nactdof[5*i+4]);}
//SFREE(fn);printf("fn freed\n");
   /* next statement allows for storing the displacements in each
      iteration: for debugging purposes */

    if((strcmp1(&filab[3],"I")==0)&&(*iout==0)){
	FORTRAN(frditeration,(co,nk,kon,ipkon,lakon,ne,v,
		ttime,ielmat,matname,mi,istep,iinc,ithermal));
    }

    /* calculating the stresses and material tangent at the
       integration points; calculating the internal forces */

    if(((ithermal[0]<=1)||(ithermal[0]>=3))&&(intpointvarm==1)){

        /* determining the element bounds in each thread */

	NNEW(neapar,ITG,num_cpus);
	NNEW(nebpar,ITG,num_cpus);
	elementcpuload(neapar,nebpar,ne,ipkon,&num_cpus);

	NNEW(fn1,double,num_cpus*mt**nk);
	NNEW(qa1,double,num_cpus*4);
	NNEW(nal,ITG,num_cpus);

	co1=co;kon1=kon;ipkon1=ipkon;lakon1=lakon;ne1=ne;v1=v;
        stx1=stx;elcon1=elcon;nelcon1=nelcon;rhcon1=rhcon;
        nrhcon1=nrhcon;alcon1=alcon;nalcon1=nalcon;alzero1=alzero;
        ielmat1=ielmat;ielorien1=ielorien;norien1=norien;orab1=orab;
        ntmat1_=ntmat_;t01=t0;t11=t1;ithermal1=ithermal;prestr1=prestr;
        iprestr1=iprestr;eme1=eme;iperturb1=iperturb;iout1=iout;
        vold1=vold;nmethod1=nmethod;veold1=veold;dtime1=dtime;
        time1=time;ttime1=ttime;plicon1=plicon;nplicon1=nplicon;
        plkcon1=plkcon;nplkcon1=nplkcon;xstateini1=xstateini;
        xstiff1=xstiff;xstate1=xstate;npmat1_=npmat_;matname1=matname;
        mi1=mi;ielas1=ielas;icmd1=icmd;ncmat1_=ncmat_;nstate1_=nstate_;
        stiini1=stiini;vini1=vini;ener1=ener;eei1=eei;enerini1=enerini;
        istep1=istep;iinc1=iinc;springarea1=springarea;reltime1=reltime;
        calcul_fn1=calcul_fn;calcul_qa1=calcul_qa;calcul_cauchy1=calcul_cauchy;
        nener1=nener;ikin1=ikin;mt1=mt;nk1=nk;ne01=ne0;thicke1=thicke;
        emeini1=emeini;pslavsurf1=pslavsurf;clearini1=clearini;
        pmastsurf1=pmastsurf;mortar1=mortar;ielprop1=ielprop;prop1=prop;
        kscale1=kscale;mpcon1=mpcon;nmpcon1=nmpcon;nmpmat1_=nmpmat_;
        pphase1=pphase;cphase1=cphase;phaseother1=phaseother;nphase1=nphase;
	phase_inf1=phase_inf;rdpcon1=rdpcon;gscon1=gscon;

        if(vini1 == NULL){
           vini1=vold;
	}

	/* calculating the stresses */

	if(((*nmethod!=4)&&(*nmethod!=5))||(iperturb[0]>1)){
		printf(" Using up to %" ITGFORMAT " cpu(s) for the stress calculation.\n\n", num_cpus);
	}

	/* create threads and wait */

	NNEW(ithread,ITG,num_cpus);
	for(i=0; i<num_cpus; i++)  {
	    ithread[i]=i;
	    pthread_create(&tid[i], NULL, (void *)resultsmechmt, (void *)&ithread[i]);
	}
	for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);

	for(i=0;i<mt**nk;i++){
	    fn[i]=fn1[i];
	}
	for(i=0;i<mt**nk;i++){
	    for(j=1;j<num_cpus;j++){
		fn[i]+=fn1[i+j*mt**nk];
	    }
	}
	SFREE(fn1);SFREE(ithread);SFREE(neapar);SFREE(nebpar);

        /* determine the internal force */

	qa[0]=qa1[0];
	for(j=1;j<num_cpus;j++){
	    qa[0]+=qa1[j*4];
	}

        /* determine the decrease of the time increment in case
           the material routine diverged */

	qa[2]=qa1[2];
        for(j=1;j<num_cpus;j++){
	    if(qa1[2+j*4]>0.){
		if(qa[2]<0.){
		    qa[2]=qa1[2+j*4];
		}else{
		    if(qa1[2+j*4]<qa[2]){qa[2]=qa1[2+j*4];}
		}
	    }
	}

        /* maximum change in creep strain increment in the
           present time increment */

	qa[3]=qa1[3];
        for(j=1;j<num_cpus;j++){
	    if(qa1[3+j*4]>0.){
		if(qa[3]<0.){
		    qa[3]=qa1[3+j*4];
		}else{
		    if(qa1[3+j*4]>qa[3]){qa[3]=qa1[3+j*4];}
		}
	    }
	}

	SFREE(qa1);

	for(j=1;j<num_cpus;j++){
	    nal[0]+=nal[j];
	}

	if(calcul_qa==1){
	    if(nal[0]>0){
		qa[0]/=nal[0];
	    }
	}
	SFREE(nal);
    }

    /* calculating the thermal flux and material tangent at the
       integration points; calculating the internal point flux */

    if((ithermal[0]>=2)&&(intpointvart==1)){

        /* determining the element bounds in each thread */

	NNEW(neapar,ITG,num_cpus);
	NNEW(nebpar,ITG,num_cpus);
	elementcpuload(neapar,nebpar,ne,ipkon,&num_cpus);

	NNEW(fn1,double,num_cpus*mt**nk);
	NNEW(qa1,double,num_cpus*4);
	NNEW(nal,ITG,num_cpus);

	co1=co;kon1=kon;ipkon1=ipkon;lakon1=lakon;v1=v;
        elcon1=elcon;nelcon1=nelcon;rhcon1=rhcon;nrhcon1=nrhcon;
	ielmat1=ielmat;ielorien1=ielorien;norien1=norien;orab1=orab;
        ntmat1_=ntmat_;t01=t0;iperturb1=iperturb;iout1=iout;vold1=vold;
        ipompc1=ipompc;nodempc1=nodempc;coefmpc1=coefmpc;nmpc1=nmpc;
        dtime1=dtime;time1=time;ttime1=ttime;plkcon1=plkcon;
        nplkcon1=nplkcon;xstateini1=xstateini;xstiff1=xstiff;
        xstate1=xstate;npmat1_=npmat_;matname1=matname;mi1=mi;
        ncmat1_=ncmat_;nstate1_=nstate_;cocon1=cocon;ncocon1=ncocon;
        qfx1=qfx;ikmpc1=ikmpc;ilmpc1=ilmpc;istep1=istep;iinc1=iinc;
        springarea1=springarea;calcul_fn1=calcul_fn;calcul_qa1=calcul_qa;
        mt1=mt;nk1=nk;shcon1=shcon;nshcon1=nshcon;ithermal1=ithermal;
        nelemload1=nelemload;nload1=nload;nmethod1=nmethod;reltime1=reltime;
        sideload1=sideload;xload1=xload;xloadold1=xloadold;
        pslavsurf1=pslavsurf;pmastsurf1=pmastsurf;mortar1=mortar;
        clearini1=clearini;plicon1=plicon;nplicon1=nplicon;ne1=ne;
        ielprop1=ielprop,prop1=prop;iponoel1=iponoel;inoel1=inoel;
    	network1=network;ipobody1=ipobody;ibody1=ibody;xbody1=xbody;
        mpcon1=mpcon,nmpcon1=nmpcon,nmpmat1_=nmpmat_;vini1=vini;
        pphase1=pphase;cphase1=cphase;phaseother1=phaseother;nphase1=nphase;
	phase_inf1=phase_inf,rdpcon1=rdpcon,gscon1=gscon;


	// updated coordinates for thermal calculation:
	//double * co2 ;
//	co2 = malloc ( sizeof(double) * (3) * (*nk) ) ;
//
//	for ( i = 0; i != 3*(*nk); i++ ) {
//	  co2[i] = co[i] ;
//	}
//
//	for ( i = 0; i != (*nk); i++ ) {
//	  co2[3*i]   = co2[3*i]   + vini[mt*i+1] ;
//	  co2[3*i+1] = co2[3*i+1] + vini[mt*i+2] ;
//	  co2[3*i+2] = co2[3*i+2] + vini[mt*i+3] ;
//	}
	// end





	/* calculating the heat flux */

	printf(" Using up to %" ITGFORMAT " cpu(s) for the heat flux calculation.\n\n", num_cpus);

	/* create threads and wait */

	NNEW(ithread,ITG,num_cpus);
	for(i=0; i<num_cpus; i++)  {
	    ithread[i]=i;
	    pthread_create(&tid[i], NULL, (void *)resultsthermmt, (void *)&ithread[i]);
	}
	for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);

	for(i=0;i<*nk;i++){
		fn[mt*i]=fn1[mt*i];
	}
	for(i=0;i<*nk;i++){
	    for(j=1;j<num_cpus;j++){
		fn[mt*i]+=fn1[mt*i+j*mt**nk];
	    }
	}
	SFREE(fn1);SFREE(ithread);SFREE(neapar);SFREE(nebpar);

        /* determine the internal concentrated heat flux */

	qa[1]=qa1[1];
	for(j=1;j<num_cpus;j++){
	    qa[1]+=qa1[1+j*4];
	}

	SFREE(qa1);

	for(j=1;j<num_cpus;j++){
	    nal[0]+=nal[j];
	}

	if(calcul_qa==1){
	    if(nal[0]>0){
		qa[1]/=nal[0];
	    }
	}
	SFREE(nal);

//	free(co2);

    }

    /* calculating the matrix system internal force vector */

    FORTRAN(resultsforc,(nk,f,fn,nactdof,ipompc,nodempc,
	    coefmpc,labmpc,nmpc,mi,fmpc,&calcul_fn,&calcul_f));

    /* storing results in the .dat file
       extrapolation of integration point values to the nodes
       interpolation of 3d results for 1d/2d elements */

    
    
    pwork = malloc ( sizeof(double) * (*mi) * (*ne) ) ;



    plasticwork ( matname,
		  xstate, nstate_, 
		  stx, ener, 
		  pwork,
		  shcon, nshcon,
		  rhcon, nrhcon, t1, 
		  v, vold, vini,
		  mi, ne, lakon, kon, co, nk,
		  elcon, nelcon, ipkon, rdpcon1, 
		  ielmat1, ntmat_, iout,
		  ttime, time ) ;

    free ( pwork ) ;
    
    

    /* use external function to calculate
    drxc (  ) ; // to be continued
    */


//    if (*nstate_==20){
//      list1=0;
//      nea=1;
//      neb=*ne;
//      FORTRAN(damagecal,(&nea,&neb,lakon,ipkon,mi,&list1,ilist1,
//          pressx,pressn,stx,co,kon,xstate,nstate_,dtime,inum,nk,ne,
//          orab,ielorien,vold,ielmat,thicke,ielprop,prop));
//    }

// Phase laten heat
    FORTRAN(phaselatentheat,(shcon,nshcon,rhcon,nrhcon,v,vold,vini,
       nk,ne,lakon,xstateini,xstate,nstate_,mi,ipkon,kon,ntmat_,cphase,
        phase_inf,nphase,iout));

/*    if(*iout>0){
           xstateini12=xstateini;
        }
*/

    FORTRAN(resultsprint,(co,nk,kon,ipkon,lakon,ne,v,stn,inum,
       stx,ielorien,norien,orab,t1,ithermal,filab,een,iperturb,fn,
       nactdof,iout,vold,nodeboun,ndirboun,nboun,nmethod,ttime,xstate,
       epn,mi,
       nstate_,ener,enern,xstaten,eei,set,nset,istartset,iendset,
       ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
       nelemload,nload,&ikin,ielmat,thicke,eme,emn,rhcon,nrhcon,shcon,
       nshcon,cocon,ncocon,ntmat_,sideload,icfd,inomat,pslavsurf,islavact,
       cdn,mortar,islavnode,nslavnode,ntie,islavsurf,time,ielprop,prop,
       veold,ne0,nmpc,ipompc,nodempc,labmpc,energyini,energy,orname,
       xload));

  return;

}

/* subroutine for multithreading of resultsmech */

void *resultsmechmt(ITG *i){

    ITG indexfn,indexqa,indexnal,nea,neb,list1,*ilist1=NULL;

    indexfn=*i*mt1**nk1;
    indexqa=*i*4;
    indexnal=*i;

    nea=neapar[*i]+1;
    neb=nebpar[*i]+1;

    list1=0;
    FORTRAN(resultsmech,(co1,kon1,ipkon1,lakon1,ne1,v1,
          stx1,elcon1,nelcon1,rhcon1,nrhcon1,alcon1,nalcon1,alzero1,
          ielmat1,ielorien1,norien1,orab1,ntmat1_,t01,t11,ithermal1,prestr1,
          iprestr1,eme1,iperturb1,&fn1[indexfn],iout1,&qa1[indexqa],vold1,
          nmethod1,
          veold1,dtime1,time1,ttime1,plicon1,nplicon1,plkcon1,nplkcon1,
          xstateini1,xstiff1,xstate1,npmat1_,matname1,mi1,ielas1,icmd1,
          ncmat1_,nstate1_,stiini1,vini1,ener1,eei1,enerini1,istep1,iinc1,
          springarea1,reltime1,&calcul_fn1,&calcul_qa1,&calcul_cauchy1,nener1,
	  &ikin1,&nal[indexnal],ne01,thicke1,emeini1,
	  pslavsurf1,pmastsurf1,mortar1,clearini1,&nea,&neb,ielprop1,prop1,
	  kscale1,&list1,ilist1,mpcon1,nmpcon1,nmpmat1_,rdpcon1,gscon1));

    return NULL;
}

/* subroutine for multithreading of resultsmech */

void *resultsthermmt(ITG *i){

    ITG indexfn,indexqa,indexnal,nea,neb;

    indexfn=*i*mt1**nk1;
    indexqa=*i*4;
    indexnal=*i;

    nea=neapar[*i]+1;
    neb=nebpar[*i]+1;


    FORTRAN(resultstherm,(co1,kon1,ipkon1,lakon1,v1, // co2 is changed from co1
	   elcon1,nelcon1,rhcon1,nrhcon1,ielmat1,ielorien1,norien1,orab1,
	   ntmat1_,t01,iperturb1,&fn1[indexfn],shcon1,nshcon1,
	   iout1,&qa1[indexqa],vold1,ipompc1,nodempc1,coefmpc1,nmpc1,
           dtime1,time1,ttime1,plkcon1,nplkcon1,xstateini1,xstiff1,xstate1,
           npmat1_,matname1,mi1,ncmat1_,nstate1_,cocon1,ncocon1,
           qfx1,ikmpc1,ilmpc1,istep1,iinc1,springarea1,
	   &calcul_fn1,&calcul_qa1,&nal[indexnal],&nea,&neb,ithermal1,
	   nelemload1,nload1,nmethod1,reltime1,sideload1,xload1,xloadold1,
	   pslavsurf1,pmastsurf1,mortar1,clearini1,plicon1,nplicon1,ielprop1,
	   prop1,iponoel1,inoel1,network1,ipobody1,xbody1,ibody1,
	   mpcon1,nmpcon1,nmpmat1_,pphase1,cphase1,phaseother1,nphase1,phase_inf1,
	   vini1,gscon1));

    return NULL;
}
