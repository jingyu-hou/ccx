!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!     
      subroutine applybounp(nodeboun,ndirboun,nboun,xbounact,
     &     ithermal,nk,iponoel,inoel,vold,vcontu,t1act,isolidsurf,
     &     nsolidsurf,xsolidsurf,nfreestream,ifreestream,turbulent,
     &     vcon,shcon,nshcon,rhcon,nrhcon,ielmat,ntmat_,physcon,v,
     &     ipompc,nodempc,coefmpc,nmpc,inomat,mi)
!     
!     applies velocity boundary conditions
!     
      implicit none
!     
      integer turbulent,
     &     nrhcon(*),mi(*),ielmat(mi(3),*),ntmat_,nodeboun(*),
     &     isolidsurf(*),
     &     ndirboun(*),nshcon(*),nk,i,nboun,node,imat,ithermal(*),
     &     iponoel(*),
     &     inoel(3,*),nsolidsurf,ifreenode,ifreestream(*),nfreestream,k,
     &     index,ipompc(*),nodempc(3,*),nmpc,ist,ndir,inomat(*),ndiri,
     &     nodei
!     
      real*8 rhcon(0:1,ntmat_,*),vold(0:mi(2),*),xbounact(*),shcon,
     &     vcontu(2,*),t1act(*),temp,r,dvi,xsolidsurf(*),reflength,
     &     refkin,reftuf,refvel,cp,vcon(0:4,*),physcon(*),v(0:mi(2),*),
     &     coefmpc(*),fixed_pres,size,correction,residu
!     
!     inserting the pressure boundary conditions
!     
      do i=1,nboun
        if(ndirboun(i).ne.4) cycle
!     
        node=nodeboun(i)
        v(4,node)=xbounact(i)-vold(4,node)
      enddo
!     
!     inserting the pressure MPC conditions
!     
      do i=1,nmpc
        ist=ipompc(i)
        ndir=nodempc(2,ist)
        if(ndir.ne.4) cycle
        node=nodempc(1,ist)
!     
!     check whether fluid MPC
!     
        imat=inomat(node)
        if(imat.eq.0) cycle
!     
        index=nodempc(3,ist)
        residu=0.d0
        if(index.ne.0) then
          do
            nodei=nodempc(1,index)
            ndiri=nodempc(2,index)
!     
            residu=residu+coefmpc(index)*
     &           (vold(ndiri,nodei)+v(ndiri,nodei))
            index=nodempc(3,index)
            if(index.eq.0) exit
          enddo
        endif
!     
        v(ndir,node)=-residu/coefmpc(ist)-vold(ndir,node)
      enddo
!     
      return
      end
