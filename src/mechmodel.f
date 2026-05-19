!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine mechmodel(elconloc,elas,emec,kode,emec0,ithermal,
     &     icmd,beta,stre,xkl,ckl,vj,xikl,vij,plconloc,xstate,xstateini,
     &     ielas,amat,t0l,t1l,dtime,time,ttime,iel,iint,nstate_,mi,
     &     iorien,pgauss,orab,eloc,mattyp,pnewdt,istep,iinc,ipkon,
     &     nmethod,iperturb,depvisc,nlgeom_undo,mpconloc,rdploc,gsloc)
!
!     kode=-1: Arruda-Boyce
!          -2: Mooney-Rivlin
!          -3: Neo-Hooke
!          -4: Ogden (N=1)
!          -5: Ogden (N=2)
!          -6: Ogden (N=3)
!          -7: Polynomial (N=1)
!          -8: Polynomial (N=2)
!          -9: Polynomial (N=3)
!          -10: Reduced Polynomial (N=1)
!          -11: Reduced Polynomial (N=2)
!          -12: Reduced Polynomial (N=3)
!          -13: Van der Waals (not implemented yet)
!          -14: Yeoh
!          -15: Hyperfoam (N=1)
!          -16: Hyperfoam (N=2)
!          -17: Hyperfoam (N=3)
!          -50: deformation plasticity
!          -51: incremental plasticity (no viscosity)
!          -52: viscoplasticity
!          -70: metal powder
!       < -100: user material routine with -kode-100 user
!               defined constants with keyword *USER MATERIAL
!
      implicit none
!
      character*80 amat
!
      integer kode,ithermal,icmd,ielas,iel,iint,nstate_,mi(*),iorien,
     &  mattyp,istep,iinc,ipkon(*),nmethod,iperturb(*),nlgeom_undo
!
      real*8 elconloc(*),elas(21),emec(*),emec0(*),beta(*),stre(*),
     &  ckl(*),vj,plconloc(*),t0l,t1l,xkl(*),xikl(*),vij,depvisc,
     &  dtime,didc(27),d2idc2(243),dibdc(27),d2ibdc2(243),
     &  dudc(9),d2udc2(81),dldc(27),d2ldc2(243),dlbdc(27),d2lbdc2(243),
     &  pgauss(3),orab(7,*),time,ttime,eloc(6),pnewdt,mpconloc(*),
     &  rdploc(10,3),gsloc(30,5)
!
      real*8 xstate(nstate_,mi(1),*),xstateini(nstate_,mi(1),*)
!
      if(kode.gt.0) then
         call linel(kode,mattyp,beta,emec,stre,elas,elconloc,
     &  iorien,orab,pgauss)
      elseif(kode.gt.-50) then
         mattyp=3
         call rubber(elconloc,elas,emec,kode,didc,d2idc2,
     &     dibdc,d2ibdc2,dudc,d2udc2,dldc,d2ldc2,dlbdc,d2lbdc2,
     &     ithermal,icmd,beta,stre)
      elseif(kode.eq.-50) then
         mattyp=3
         call defplas(elconloc,elas,emec,ithermal,icmd,beta,stre,
     &     ckl,vj)
      elseif(kode.gt.-60) then
         mattyp=3
         if(iperturb(2).eq.1) then
            call incplas(elconloc,plconloc,xstate,xstateini,elas,emec,
     &           ithermal,icmd,beta,stre,vj,kode,ielas,amat,t1l,dtime,
     &           time,ttime,iel,iint,nstate_,mi(1),eloc,pgauss,nmethod,
     &           pnewdt,depvisc)
         else
            call incplas_lin(elconloc,plconloc,xstate,xstateini,elas,
     &           emec,
     &           ithermal,icmd,beta,stre,vj,kode,ielas,amat,t1l,dtime,
     &           time,ttime,iel,iint,nstate_,mi(1),eloc,pgauss,nmethod,
     &           pnewdt,depvisc)
         endif
      elseif(kode.eq.-70) then
         mattyp=3
         call transformation(amat,iel,iint,kode,elconloc,plconloc,
     &        mpconloc,emec,emec0,ckl,
     &        beta,xikl,vij,xkl,vj,ithermal,t1l,dtime,time,ttime,icmd,
     &        ielas,mi(1),nstate_,xstateini,xstate,stre,elas,
     &        pgauss,istep,iinc,pnewdt,nmethod,iperturb)
      elseif(kode.gt.-100) then
         mattyp=3
         call transformation(amat,iel,iint,kode,elconloc,plconloc,
     &        mpconloc,emec,emec0,ckl,
     &        beta,xikl,vij,xkl,vj,ithermal,t1l,dtime,time,ttime,icmd,
     &        ielas,mi(1),nstate_,xstateini,xstate,stre,elas,
     &        pgauss,istep,iinc,pnewdt,nmethod,iperturb,rdploc,gsloc)
      else
         mattyp=3
         call umat_main(amat,iel,iint,kode,elconloc,emec,emec0,beta,
     &        xikl,vij,xkl,vj,ithermal,t0l,t1l,dtime,time,ttime,icmd,
     &        ielas,mi(1),nstate_,xstateini,xstate,stre,elas,iorien,
     &        pgauss,orab,pnewdt,istep,iinc,ipkon,nmethod,iperturb,
     &        depvisc,eloc,nlgeom_undo)
      endif
      
!      if (iint.eq.1) then
!        open(1111, file='mechmodel.txt')
!        write(1111,*)"iinc", iinc
!        write(1111,*)"stress"
!        write(1111,*)stre(1:6)
!        write(1111,*)elas
!        write(1111,*)""
!      endif

      return
      end
