/*                                                                              */
/*      WeICME (Wedge Integrated Computational Materials Engineering)           */
/*                 - A 3-dimensional finite element program.                    */   
/*     Developed and maintained by Shenzhen Wedge Central                       */
/*    South Research Institute co., Ltd., Shenzhen, China                       */   
/*     Copy Right 2019-2023.                                                      */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "WeICME.h"
#ifdef SPOOLES
   #include "spooles.h"
#endif
#ifdef SGI
   #include "sgi.h"
#endif
#ifdef TAUCS
   #include "tau.h"
#endif
#ifdef PARDISO
#include "pardiso.h"
#endif

#define max(a,b) ((a) >= (b) ? (a) : (b))

void castfvm(double **cop, ITG *nk, ITG **konp, ITG **ipkonp, char **lakonp,
	     ITG *ne, 
	     ITG *nodeboun, ITG *ndirboun, double *xboun, ITG *nboun, 
	     ITG **ipompcp, ITG **nodempcp, double **coefmpcp, char **labmpcp,
             ITG *nmpc, 
	     ITG *nodeforc, ITG *ndirforc,double *xforc, ITG *nforc, 
	     ITG **nelemloadp, char **sideloadp, double *xload,ITG *nload, 
	     ITG *nactdof, 
	     ITG **icolp, ITG *jq, ITG **irowp, ITG *neq, ITG *nzl, 
	     ITG *nmethod, ITG **ikmpcp, ITG **ilmpcp, ITG *ikboun, 
	     ITG *ilboun,
             double *elcon, ITG *nelcon, double *rhcon, ITG *nrhcon,
	     double *alcon, ITG *nalcon, double *alzero, ITG **ielmatp,
	     ITG **ielorienp, ITG *norien, double *orab, ITG *ntmat_,
	     double *t0, double *t1, double *t1old, 
	     ITG *ithermal,double *prestr, ITG *iprestr, 
	     double **voldp,ITG *iperturb, double *sti, ITG *nzs,  
	     ITG *kode, char *filab, 
             ITG *idrct, ITG *jmax, ITG *jout, double *timepar,
             double *eme,
	     double *xbounold, double *xforcold, double *xloadold,
             double *veold, double *accold,
	     char *amname, double *amta, ITG *namta, ITG *nam,
             ITG *iamforc, ITG **iamloadp,
             ITG *iamt1, double *alpha, ITG *iexpl,
	     ITG *iamboun, double *plicon, ITG *nplicon, double *plkcon,
             ITG *nplkcon,
             double **xstatep, ITG *npmat_, ITG *istep, double *ttime,
             char *matname, double *qaold, ITG *mi,
             ITG *isolver, ITG *ncmat_, ITG *nstate_, ITG *iumat,
             double *cs, ITG *mcs, ITG *nkon, double **enerp, ITG *mpcinfo,
             char *output,
             double *shcon, ITG *nshcon, double *cocon, ITG *ncocon,
             double *physcon, ITG *nflow, double *ctrl, 
             char *set, ITG *nset, ITG *istartset,
             ITG *iendset, ITG *ialset, ITG *nprint, char *prlab,
             char *prset, ITG *nener,ITG *ikforc,ITG *ilforc, double *trab, 
             ITG *inotr, ITG *ntrans, double **fmpcp, char *cbody,
             ITG *ibody, double *xbody, ITG *nbody, double *xbodyold,
             ITG *ielprop, double *prop, ITG *ntie, char *tieset,
	     ITG *itpamp, ITG *iviewfile, char *jobnamec, double *tietol,
	     ITG *nslavs, double *thicke, ITG *ics, 
	     ITG *nintpoint,ITG *mortar,ITG *ifacecount,char *typeboun,
	     ITG **islavsurfp,double **pslavsurfp,double **clearinip,
	     ITG *nmat,double *xmodal,ITG *iaxial,ITG *inext,ITG *nprop,
	     ITG *network,char *orname,double *mpcon, ITG *nmpcon,ITG *nmpmat_,
             double*pphase,double *cphase,double *phaseother, ITG *nphase,
	     ITG *phase_inf,double *flcon,ITG *nflcon,double *lh,double *rdpcon,
	     double *gscon){

  char description[13]="            ",*lakon=NULL,jobnamef[396]="",
      *sideface=NULL,*labmpc=NULL,fnffrd[132]="",*lakonf=NULL,
      *sideloadref=NULL,*sideload=NULL,stiffmatrix[132]=""; 
 
  ITG *inum=NULL,k,iout=0,icntrl,iinc=0,jprint=0,iit=-1,jnz=0,
      icutb=0,istab=0,ifreebody,uncoupled,n1,n2,itruecontact,
      iperturb_sav[2],ilin,*icol=NULL,*irow=NULL,ielas=0,icmd=0,
      memmpc_,mpcfree,icascade,maxlenmpc,*nodempc=NULL,*iaux=NULL,
      *nodempcref=NULL,memmpcref_,mpcfreeref,*itg=NULL,*ineighe=NULL,
      *ieg=NULL,ntg=0,ntr,*kontri=NULL,*nloadtr=NULL,idamping=0,
      *ipiv=NULL,ntri,newstep,mode=-1,noddiam=-1,nasym=0,im,
      ntrit,*inocs=NULL,inewton=0,*ipobody=NULL,*nacteq=NULL,
      *nactdog=NULL,nteq,*itietri=NULL,*koncont=NULL,istrainfree=0,
      ncont,ne0,nkon0,*ipkon=NULL,*kon=NULL,*ielorien=NULL,
      *ielmat=NULL,itp=0,symmetryflag=0,inputformat=0,kscale=1,
      *iruc=NULL,iitterm=0,iturbulent,ngraph=1,ismallsliding=0,
      *ipompc=NULL,*ikmpc=NULL,*ilmpc=NULL,i0ref,irref,icref,
      *itiefac=NULL,*islavsurf=NULL,*islavnode=NULL,*imastnode=NULL,
      *nslavnode=NULL,*nmastnode=NULL,*imastop=NULL,imat,
      *iponoels=NULL,*inoels=NULL,*islavsurfold=NULL,maxlenmpcref,
      *islavact=NULL,mt=mi[1]+1,*nactdofinv=NULL,*ipe=NULL, 
      *ime=NULL,*ikactmech=NULL,nactmech,inode,idir,neold,neini,
      iemchange=0,nzsrad,*mast1rad=NULL,*irowrad=NULL,*icolrad=NULL,
      *jqrad=NULL,*ipointerrad=NULL,*integerglob=NULL,negpres=0,
      mass[2]={0,0}, stiffness=1, buckling=0, rhsi=1, intscheme=0,idiscon=0,
      coriolis=0,*ipneigh=NULL,*neigh=NULL,maxprevcontel,nslavs_prev_step,
      *nelemface=NULL,*ipoface=NULL,*nodface=NULL,*ifreestream=NULL,iex,
      *isolidsurf=NULL,*neighsolidsurf=NULL,*iponoel=NULL,*inoel=NULL,
      nef,nface,nfreestream,nsolidsurf,i,icfd=0,id,*neij=NULL,
      node,networknode,iflagact=0,*nodorig=NULL,*ipivr=NULL,iglob=0,
      *inomat=NULL,*ipnei=NULL,ntrimax,*nx=NULL,*ny=NULL,*nz=NULL,
      *neifa=NULL,*neiel=NULL,*ielfa=NULL,*ifaext=NULL,nflnei,nfaext,
      idampingwithoutcontact=0,*nactdoh=NULL,*nactdohinv=NULL,*ipkonf=NULL,
      *ielmatf=NULL,*ielorienf=NULL,ialeatoric=0,nloadref,isym,
      *nelemloadref=NULL,*iamloadref=NULL,*idefload=NULL,nload_,
      *nelemload=NULL,*iamload=NULL,ncontacts=0,inccontact=0,nrhs=1,
      j=0,*ifatie=NULL,n,inoelsize=0,isensitivity=0,*istartblk=NULL,
      *iendblk=NULL,*nblket=NULL,*nblkze=NULL,nblk,*konf=NULL,*ielblk=NULL,
      *iwork=NULL,nelt,lrgw,*igwk=NULL,itol,itmax,iter,ierr,iunit,ligw,inoelfree,
      singularflag,initial;

  double *stn=NULL,*v=NULL,*een=NULL,cam[5],*epn=NULL,*cg=NULL,
         *cdn=NULL,*vel=NULL,*velo=NULL,*veloo=NULL,*vfa=NULL,*pslavsurfold=NULL,
         *f=NULL,*fn=NULL,qa[4]={0.,0.,-1.,0.},qam[2]={0.,0.},dtheta,theta,
	 err,ram[8]={0.,0.,0.,0.,0.,0.,0.,0.},*areaslav=NULL,
         *springarea=NULL,ram1[8]={0.,0.,0.,0.,0.,0.,0.,0.},
	 ram2[8]={0.,0.,0.,0.,0.,0.,0.,0.},deltmx,ptime,smaxls,sminls,
         uam[2]={0.,0.},*vini=NULL,*ac=NULL,qa0,qau,ea,*straight=NULL,
	 *t1act=NULL,qamold[2],*xbounact=NULL,*bc=NULL,
	 *xforcact=NULL,*xloadact=NULL,*fext=NULL,*clearini=NULL,
         reltime,time,bet=0.,gam=0.,*aux2=NULL,dtime,*fini=NULL,
         *fextini=NULL,*veini=NULL,*accini=NULL,*xstateini=NULL,
	 *ampli=NULL,scal1,*eei=NULL,*t1ini=NULL,pressureratio,
         *xbounini=NULL,dev,*xstiff=NULL,*stx=NULL,*stiini=NULL,
         *enern=NULL,*coefmpc=NULL,*aux=NULL,*xstaten=NULL,
	 *coefmpcref=NULL,*enerini=NULL,*emn=NULL,alpham,betam,
	 *tarea=NULL,*tenv=NULL,*erad=NULL,*fnr=NULL,*fni=NULL,
	 *adview=NULL,*auview=NULL,*qfx=NULL,*cvini=NULL,*cv=NULL,
         *qfn=NULL,*co=NULL,*vold=NULL,*fenv=NULL,sigma=0.,
         *xbodyact=NULL,*cgr=NULL,dthetaref, *vr=NULL,*vi=NULL,
	 *stnr=NULL,*stni=NULL,*vmax=NULL,*stnmax=NULL,*fmpc=NULL,*ener=NULL,
	 *f_cm=NULL, *f_cs=NULL,*adc=NULL,*auc=NULL,*res=NULL,
	 *xstate=NULL,*eenmax=NULL,*adrad=NULL,*aurad=NULL,*bcr=NULL,
	 *xmastnor=NULL,*emeini=NULL,*tinc,*tper,*tmin,*tmax,*tincf,
	 *doubleglob=NULL,*xnoels=NULL,*au=NULL,*resold=NULL,
	 *ad=NULL,*b=NULL,*aub=NULL,*adb=NULL,*pslavsurf=NULL,*pmastsurf=NULL,
	 *x=NULL,*y=NULL,*z=NULL,*xo=NULL,sum1,sum2,flinesearch,
	 *yo=NULL,*zo=NULL,*cdnr=NULL,*cdni=NULL,*fnext=NULL,*fnextini=NULL,
	 allwk=0.,allwkini,energy[4]={0.,0.,0.,0.},energyini[4],
	 energyref,denergymax,dtcont,dtvol,wavespeed[*nmat],emax,r_abs,
         enetoll,dampwk=0.,dampwkini=0.,temax,*tmp=NULL,energystartstep[4],
	 sizemaxinc,*adblump=NULL,*adcpy=NULL,*aucpy=NULL,*rwork=NULL,
	 *sol=NULL,*rgwk=NULL,tol,*sb=NULL,*sx=NULL,*vcontu=NULL;
	 
  // MPADD: initialize rmin to the tolerance
  /*enetoll=0.02;
  r_abs=0.0;
  emax=0.0;
  singularflag=0;*/
  // MPADD end


#ifdef SGI
  ITG token;
#endif
  
  icol=*icolp;irow=*irowp;co=*cop;vold=*voldp;
  ipkon=*ipkonp;lakon=*lakonp;kon=*konp;ielorien=*ielorienp;
  ielmat=*ielmatp;ener=*enerp;xstate=*xstatep;
  
  ipompc=*ipompcp;labmpc=*labmpcp;ikmpc=*ikmpcp;ilmpc=*ilmpcp;
  fmpc=*fmpcp;nodempc=*nodempcp;coefmpc=*coefmpcp;nelemload=*nelemloadp;
  iamload=*iamloadp;sideload=*sideloadp;

  islavsurf=*islavsurfp;pslavsurf=*pslavsurfp;clearini=*clearinip;

  tinc=&timepar[0];
  tper=&timepar[1];
  tmin=&timepar[2];
  tmax=&timepar[3];
  tincf=&timepar[4];
  
  /*if(*ithermal==4){
      uncoupled=1;
      *ithermal=3;
  }else{
      uncoupled=0;
  }*/
  
  /*if(*mortar!=1){
      maxprevcontel=*nslavs;
  }else if(*mortar==1){
      maxprevcontel=*nintpoint;
      if(*nstate_!=0){
	  if(maxprevcontel!=0){
	      NNEW(islavsurfold,ITG,2**ifacecount+2);
	      NNEW(pslavsurfold,double,3**nintpoint);
	      memcpy(&islavsurfold[0],&islavsurf[0],
		     sizeof(ITG)*(2**ifacecount+2));
	      memcpy(&pslavsurfold[0],&pslavsurf[0],
		     sizeof(double)*(3**nintpoint));
	  }
      }
  }*/
  //nslavs_prev_step=*nslavs;

  /* turbulence model 
     iturbulent==0: laminar
     iturbulent==1: k-epsilon
     iturbulent==2: q-omega
     iturbulent==3: SST */
  
  iturbulent=(ITG)physcon[8];
  
  for(k=0;k<3;k++){
      strcpy1(&jobnamef[k*132],&jobnamec[k*132],132);
  }
  
  //qa0=ctrl[20];qau=ctrl[21];ea=ctrl[23];deltmx=ctrl[26];
  //i0ref=ctrl[0];irref=ctrl[1];icref=ctrl[3];

  //sminls=ctrl[28];smaxls=ctrl[29];
  
  memmpc_=mpcinfo[0];mpcfree=mpcinfo[1];//icascade=mpcinfo[2];
  maxlenmpc=mpcinfo[3];

  //alpham=xmodal[0];
  //betam=xmodal[1];

  /* check whether, for a dynamic calculation, damping is involved */

  /*if(*nmethod==4){
      if(*iexpl<=1){
	  if((fabs(alpham)>1.e-30)||(fabs(betam)>1.e-30)){
	      idamping=1;
	      idampingwithoutcontact=1;
	  }else{
	      for(i=0;i<*ne;i++){
		  if(ipkon[i]<0) continue;
		  if(strcmp1(&lakon[i*8],"ED")==0){
		      idamping=1;idampingwithoutcontact=1;break;
		  }
	      }
	  }
      }
  }*/

  /* check whether a sensitivity step may follow (whether design variables
     were defined */

  /*for(i=0;i<*ntie;i++){
      if(strcmp1(&tieset[i*243+80],"D")==0){
	  isensitivity=1;
	  NNEW(adcpy,double,neq[1]);
          /* no asymmetric matrices allowed for sensitivity 
	  NNEW(aucpy,double,nzs[1]);
	  break;
      }
  }*/

  /*if((icascade==2)&&(*iexpl>=2)){
      printf(" *ERROR in nonlingeo: linear and nonlinear MPC's depend on each other\n");
      printf("        This is not allowed in a explicit dynamic calculation\n");
      FORTRAN(stop,());
  }*/
  
  /* check whether the submodel is meant for a fluid-structure
     interaction */
  
  strcpy(fnffrd,jobnamec);
  strcat(fnffrd,"f.frd");
  if((jobnamec[396]!='\0')&&(strcmp1(&jobnamec[396],fnffrd)==0)){
      
      /* fluid-structure interaction: wait till after the compfluid
         call */
      
      NNEW(integerglob,ITG,1);
      NNEW(doubleglob,double,1);
  }else{
      
      /* determining the global values to be used as boundary conditions
	 for a submodel */
      
      getglobalresults(jobnamec,&integerglob,&doubleglob,nboun,iamboun,xboun,
		       nload,sideload,iamload,&iglob,nforc,iamforc,xforc,
                       ithermal,nk,t1,iamt1);
  }

  /* allocating a field for the stiffness matrix */
  
  //NNEW(xstiff,double,(long long)27*mi[0]**ne);
  
  /* allocating force fields */
  
  //NNEW(f,double,neq[1]);
  //NNEW(fext,double,neq[1]);
  
  //NNEW(b,double,neq[1]);
  //NNEW(vini,double,mt**nk);
  
  //NNEW(aux,double,7*maxlenmpc);
  //NNEW(iaux,ITG,2*maxlenmpc);
  
  /* allocating fields for the actual external loading */
  
  NNEW(xbounact,double,*nboun);
  NNEW(xbounini,double,*nboun);
  for(k=0;k<*nboun;++k){xbounact[k]=xbounold[k];}
  NNEW(xforcact,double,*nforc);
  NNEW(xloadact,double,2**nload);
  NNEW(xbodyact,double,7**nbody);
  /* copying the rotation axis and/or acceleration vector */
  for(k=0;k<7**nbody;k++){xbodyact[k]=xbody[k];}
  
  /* assigning the body forces to the elements */ 
  
  if(*nbody>0){
      ifreebody=*ne+1;
      NNEW(ipobody,ITG,2*ifreebody**nbody);
      for(k=1;k<=*nbody;k++){
        FORTRAN(bodyforce,(cbody,ibody,ipobody,nbody,set,istartset,
                    iendset,ialset,&inewton,nset,&ifreebody,&k));
        RENEW(ipobody,ITG,2*(*ne+ifreebody));
      }
      RENEW(ipobody,ITG,2*(ifreebody-1));
      if(inewton==1){NNEW(cgr,double,4**ne);}
  }
  
  /* for mechanical calculations: updating boundary conditions
     calculated in a previous thermal step */
  
  /*if(*ithermal<2) FORTRAN(gasmechbc,(vold,nload,sideload,
				     nelemload,xload,mi));*/
  
  /* for thermal calculations: forced convection and cavity
     radiation*/
  
  if(*ithermal>1){
      NNEW(itg,ITG,*nload+3**nflow);
      NNEW(ieg,ITG,*nflow);
      /* max 6 triangles per face, 4 entries per triangle */
      NNEW(kontri,ITG,24**nload);
      NNEW(nloadtr,ITG,*nload);
      NNEW(nacteq,ITG,4**nk);
      NNEW(nactdog,ITG,4**nk);
      NNEW(v,double,mt**nk);
      FORTRAN(envtemp,(itg,ieg,&ntg,&ntr,sideload,nelemload,
		       ipkon,kon,lakon,ielmat,ne,nload,
                       kontri,&ntri,nloadtr,nflow,ndirboun,nactdog,
                       nodeboun,nacteq,nboun,ielprop,prop,&nteq,
                       v,network,physcon,shcon,ntmat_,co,
                       vold,set,nshcon,rhcon,nrhcon,mi,nmpc,nodempc,
                       ipompc,labmpc,ikboun,&nasym,ttime,&time,iaxial));
      SFREE(v);
      
      if((*mcs>0)&&(ntr>0)){
        NNEW(inocs,ITG,*nk);
        radcyc(nk,kon,ipkon,lakon,ne,cs,mcs,nkon,ialset,istartset,
            iendset,&kontri,&ntri,&co,&vold,&ntrit,inocs,mi);
      }
      else{ntrit=ntri;}
      
      nzsrad=100*ntr;
      NNEW(mast1rad,ITG,nzsrad);
      NNEW(irowrad,ITG,nzsrad);
      NNEW(icolrad,ITG,ntr);
      NNEW(jqrad,ITG,ntr+1);
      NNEW(ipointerrad,ITG,ntr);
      
      if(ntr>0){
        mastructrad(&ntr,nloadtr,sideload,ipointerrad,
                &mast1rad,&irowrad,&nzsrad,
                jqrad,icolrad);
      }
      
      /* determine the network elements belonging to a given node (for usage
         in user subroutine film */

      if((*network>0)||(ntg>0)){
        NNEW(iponoel,ITG,*nk);
        NNEW(inoel,ITG,2**nkon);
        if(*network>0){
            FORTRAN(networkelementpernode,(iponoel,inoel,lakon,ipkon,kon,
            &inoelsize,nflow,ieg,ne,network));
        }
        RENEW(inoel,ITG,2*inoelsize);
      }

      SFREE(ipointerrad);SFREE(mast1rad);
      RENEW(irowrad,ITG,nzsrad);
      
      RENEW(itg,ITG,ntg);
      NNEW(ineighe,ITG,ntg);
      RENEW(kontri,ITG,4*ntrit);
      RENEW(nloadtr,ITG,ntr);
      
      NNEW(adview,double,ntr);
      NNEW(auview,double,2*nzsrad);
      NNEW(tarea,double,ntr);
      NNEW(tenv,double,ntr);
      NNEW(fenv,double,ntr);
      NNEW(erad,double,ntr);
      
      NNEW(ac,double,nteq*nteq);
      NNEW(bc,double,nteq);
      NNEW(ipiv,ITG,nteq);
      NNEW(adrad,double,ntr);
      NNEW(aurad,double,2*nzsrad);
      NNEW(bcr,double,ntr);
      NNEW(ipivr,ITG,ntr);
  }

    /* check for fluid elements */
  
  NNEW(nactdoh,ITG,*ne);
  NNEW(nactdohinv,ITG,*ne);
  nef=0;
  for(i=0;i<*ne;++i){
      if(ipkon[i]<0) continue;
      if(strcmp1(&lakon[8*i],"F")==0){
	  icfd=1;nactdohinv[nef]=i+1;++nef;nactdoh[i]=nef;}
      if(istrainfree==0){
	  if(ielmat[i]<0){istrainfree=1;}
      }
  }
  if((icfd==1)&&(iturbulent>=10)){
    if(iturbulent<20){
      iturbulent=iturbulent-10;
      icfd=2;
    }else{
      iturbulent=iturbulent-20;
      icfd=3;
     }
  }
  
  if(icfd==1){
    NNEW(vel,double,8*nef);
    NNEW(velo,double,8*nef);
    NNEW(veloo,double,8*nef);
    /* checking block structures (CFD calculations) */

    NNEW(ielfa,ITG,24*nef);
    NNEW(nodface,ITG,5*6*nef);//SFREE(nodface);
    NNEW(neiel,ITG,6*nef);
    NNEW(neij,ITG,6*nef);
    NNEW(neifa,ITG,6*nef);
    NNEW(ipoface,ITG,*nk);
    NNEW(ipnei,ITG,*ne+1);
    NNEW(konf,ITG,*nkon);
    memcpy(&konf[0],&kon[0],sizeof(ITG)**nkon);
    DMEMSET(ipoface,0,*nk,0);
    DMEMSET(neiel,0,6*nef,0);
    DMEMSET(ielfa,0,24*nef,0);

    /* gathering topological information (CFD calculations) */

    RENEW(nactdohinv,ITG,nef);
    NNEW(ipkonf,ITG,nef);
    NNEW(lakonf,char,8*nef);
    NNEW(ielmatf,ITG,mi[2]*nef);
    if(*norien>0){NNEW(ielorienf,ITG,mi[2]*nef);}
    NNEW(ifatie,ITG,6*nef);
    NNEW(ifaext,ITG,6*nef);
    NNEW(isolidsurf,ITG,6*nef);
    NNEW(vfa,double,8*6*nef);

    n=0;
    for(i=0;i<*mcs;i++){
      if(floor(cs[17*i+3])>n){n=floor(cs[17*i+3]);}
    }
    NNEW(xo,double,n);NNEW(yo,double,n);NNEW(zo,double,n);	    
    NNEW(x,double,n);NNEW(y,double,n);NNEW(z,double,n);	   
    NNEW(nx,ITG,n);NNEW(ny,ITG,n);NNEW(nz,ITG,n);

    FORTRAN(topocfd,(ne,ipkon,konf,lakon,ipnei,neifa,neiel,ipoface,
		     nodface,ielfa,&nflnei,&nface,ifaext,&nfaext,
		     isolidsurf,&nsolidsurf,set,nset,istartset,iendset,ialset,
		     vel,vold,mi,neij,&nef,nactdoh,ipkonf,lakonf,ielmatf,ielmat,
		     ielorienf,ielorien,norien,cs,mcs,tieset,x,y,z,xo,yo,zo,
		     nx,ny,nz,co,ifatie,velo,veloo,&initial));

    SFREE(xo);SFREE(yo);SFREE(zo);SFREE(x);SFREE(y);SFREE(z);
    SFREE(nx);SFREE(ny);SFREE(nz);

    SFREE(ipoface);
    SFREE(nodface);
    RENEW(neifa,ITG,nflnei);
    RENEW(neiel,ITG,nflnei);
    RENEW(neij,ITG,nflnei);
    RENEW(ielfa,ITG,4*nface);
    RENEW(ifatie,ITG,nface);
    RENEW(ifaext,ITG,nfaext);
    RENEW(isolidsurf,ITG,nsolidsurf);
    RENEW(vfa,double,8*nface);
    RENEW(ipnei,ITG,nef+1);
  }

  if(*ithermal>1){NNEW(qfx,double,3*mi[0]**ne);}

  if((*ithermal==1)||(*ithermal>=3)){
      NNEW(t1ini,double,*nk);
      NNEW(t1act,double,*nk);
      for(k=0;k<*nk;++k){t1act[k]=t1old[k];}
  }

  FORTRAN(checktime,(itpamp,namta,tinc,ttime,amta,tmin,inext,&itp,istep,tper));
  dtheta=(*tinc)/(*tper);

  /* taking care of a small increment at the end of the step
     for face-to-face penalty contact */

  //dthetaref=dtheta;
  /*if((dtheta<=1.e-6)&&(*iexpl<=1)){
      printf("\n *ERROR in nonlingeo\n");
      printf(" increment size smaller than one millionth of step size\n");
      printf(" increase increment size\n\n");
  }*/
  //*tmin=*tmin/(*tper);
  //*tmax=*tmax/(*tper);
  theta=0.;
  
  /* calculating an initial flux norm */
  
  /*if(*ithermal!=2){
      if(qau>1.e-10){qam[0]=qau;}
      else if(qa0>1.e-10){qam[0]=qa0;}
      else if(qa[0]>1.e-10){qam[0]=qa[0];}
      else {qam[0]=1.e-2;}
//      else {qam[0]=1.e0;}
  }
  if(*ithermal>1){
      if(qau>1.e-10){qam[1]=qau;}
      else if(qa0>1.e-10){qam[1]=qa0;}
      else if(qa[1]>1.e-10){qam[1]=qa[1];}
      else {qam[1]=1.e-2;}
//      else {qam[1]=1.e0;}
  }*/

      reltime=theta+dtheta;
      time=reltime**tper;
      dtime=dtheta**tper;
      
      FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,xloadold,xload,
	      xloadact,iamload,nload,ibody,xbody,nbody,xbodyold,xbodyact,
	      t1old,t1,t1act,iamt1,nk,amta,
	      namta,nam,ampli,&time,&reltime,ttime,&dtime,ithermal,nmethod,
              xbounold,xboun,xbounact,iamboun,nboun,
              nodeboun,ndirboun,nodeforc,ndirforc,istep,&iinc,
	      co,vold,itg,&ntg,amname,ikboun,ilboun,nelemload,sideload,mi,
              ntrans,trab,inotr,veold,integerglob,doubleglob,tieset,istartset,
              iendset,ialset,ntie,nmpc,ipompc,ikmpc,ilmpc,nodempc,coefmpc,
              ipobody,iponoel,inoel,ipkon,kon,ielprop,prop,ielmat,
              shcon,nshcon,rhcon,nrhcon,cocon,ncocon,ntmat_,lakon));

    if(icfd==1){
      compfluid(&co,nk,&ipkonf,konf,&lakonf,&sideface,
		ifreestream,&nfreestream,isolidsurf,neighsolidsurf,&nsolidsurf,
		nshcon,shcon,nrhcon,rhcon,&vold,ntmat_,nodeboun,
		ndirboun,nboun,ipompc,nodempc,nmpc,ikmpc,ilmpc,ithermal,
		ikboun,ilboun,&iturbulent,isolver,iexpl,ttime,
		&time,&dtime,nodeforc,ndirforc,xforc,nforc,nelemload,sideload,
		xload,nload,xbody,ipobody,nbody,ielmatf,matname,mi,ncmat_,
		physcon,istep,&iinc,ibody,xloadold,xboun,coefmpc,
		nmethod,xforcold,xforcact,iamforc,iamload,xbodyold,xbodyact,
		t1old,t1,t1act,iamt1,amta,namta,nam,ampli,xbounold,xbounact,
		iamboun,itg,&ntg,amname,t0,&nelemface,&nface,cocon,ncocon,
		xloadact,
		tper,jmax,jout,set,nset,istartset,iendset,ialset,prset,prlab,
		nprint,trab,inotr,ntrans,filab,labmpc,sti,norien,orab,jobnamef,
		tieset,ntie,mcs,ics,cs,nkon,&mpcfree,&memmpc_,fmpc,&nef,&inomat,
		qfx,neifa,neiel,ielfa,ifaext,vfa,vel,ipnei,&nflnei,&nfaext,
		typeboun,neij,tincf,nactdoh,nactdohinv,ielorienf,jobnamec,
		ifatie,nstate_,xstate,orname,kon,ctrl,kode,velo,veloo,
		&initial);

	  /* determining the global values to be used as boundary conditions
	     for a submodel */
	  
	  /*
	  getglobalresults(jobnamec,&integerglob,&doubleglob,nboun,iamboun,
			   xboun,nload,sideload,iamload,&iglob,nforc,iamforc,
                           xforc,ithermal,nk,t1,iamt1);
      }*/
    }

  //SFREE(f);SFREE(b);SFREE(fext);
  SFREE(xbounact);SFREE(xforcact);SFREE(xloadact);SFREE(xbodyact);
  if(*nbody>0) SFREE(ipobody);if(inewton==1){SFREE(cgr);}
  SFREE(ampli);SFREE(xbounini);//SFREE(xstiff);
  SFREE(integerglob);SFREE(doubleglob);
  if((*ithermal==1)||(*ithermal>=3)){SFREE(t1act);SFREE(t1ini);}

  if(icfd==1){
      SFREE(neifa);SFREE(neiel);SFREE(neij);SFREE(ielfa);SFREE(ifaext);
      SFREE(vfa);SFREE(nactdoh);SFREE(nactdohinv);SFREE(konf);
      SFREE(ipkonf);SFREE(lakonf);SFREE(ielmatf);SFREE(ifatie);
      SFREE(ipnei);SFREE(isolidsurf);
      if(*norien>0) SFREE(ielorienf);
  }
  return;
}      