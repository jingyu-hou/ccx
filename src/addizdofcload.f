!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!  
      subroutine addizdofcload(nodeforc,ndirforc,nactdof,mi,izdof,
     &  nzdof,iforc,iznode,nznode,nk,imdnode,nmdnode,xforc,ntrans,
     &  inotr)
!
!     adds the dof in which a point force was applied to iznode, izdof
!     and to ** imdnode if user-defined load **
!     (needed in dyna.c and steadystate.c)
!
      implicit none
!
      integer nodeforc(2,*),ndirforc(*),iforc,node,j,jdof,mi(*),nk,
     &  nactdof(0:mi(2),*),izdof(*),nzdof,iznode(*),nznode,nodebasis,
     &  imdnode(*),nmdnode,ntrans,itr,inotr(2,*)
!
      real*8 xforc(*)
!
      node=nodeforc(1,iforc)
!
!     adding the nodes in the basis sector to iznode
!
      nodebasis=mod(node,nk)
      call addimd(iznode,nznode,nodebasis)
c!
c!     user-defined load
c!
c      if((xforc(iforc).lt.1.2357111318d0).and.
c     &     (xforc(iforc).gt.1.2357111316d0)) then
c         call addimd(imdnode,nmdnode,node)
c      endif
!
      if(ntrans.eq.0) then
         itr=0
      else
         itr=inotr(1,node)
      endif
!
      if(itr.eq.0) then
!        
!        no local transformation
!
         j=ndirforc(iforc)
!
!        C-convention!
!
         jdof=nactdof(j,node)-1
         if(jdof.gt.0) call addimd(izdof,nzdof,jdof)
      else
!
!        local transformation: loop over all dofs
!
         do j=1,3
            jdof=nactdof(j,node)-1
            if(jdof.gt.0) call addimd(izdof,nzdof,jdof)
         enddo
      endif
!
      return
      end

