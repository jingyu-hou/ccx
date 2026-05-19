!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!  
      subroutine addimdnodecload(nodeforc,iforc,imdnode,nmdnode,xforc,
     &              ikmpc,ilmpc,ipompc,
     &              nodempc,nmpc,imddof,nmddof,
     &              nactdof,mi,imdmpc,nmdmpc,imdboun,nmdboun,
     &              ikboun,nboun,ilboun,ithermal)
!
!     adds the dof in which a user-defined point force was applied to imdnode
!     (needed in dyna.c and steadystate.c)
!
      implicit none
!
      integer nodeforc(2,*),iforc,node,imdnode(*),nmdnode,ikmpc(*),
     &  ilmpc(*),ipompc(*),nodempc(3,*),nmpc,imddof(*),nmddof,
     &  mi(*),nactdof(0:mi(2),*),imdmpc(*),nmdmpc,imdboun(*),nmdboun,
     &  ikboun(*),nboun,ilboun(*),ithermal,k
!
      real*8 xforc(*)
!
      node=nodeforc(1,iforc)
!
!     user-defined load
!
      if((xforc(iforc).lt.1.2357111318d0).and.
     &     (xforc(iforc).gt.1.2357111316d0)) then
!
         call addimd(imdnode,nmdnode,node)
!
!        add the degrees of freedom corresponding to the node
!
         if(ithermal.ne.2) then
            do k=1,3
               call addimdnodedof(node,k,ikmpc,ilmpc,ipompc,
     &              nodempc,nmpc,imdnode,nmdnode,imddof,nmddof,
     &              nactdof,mi,imdmpc,nmdmpc,imdboun,nmdboun,
     &              ikboun,nboun,ilboun)
            enddo
         else
            k=0
            call addimdnodedof(node,k,ikmpc,ilmpc,ipompc,
     &           nodempc,nmpc,imdnode,nmdnode,imddof,nmddof,
     &           nactdof,mi,imdmpc,nmdmpc,imdboun,nmdboun,ikboun,
     &           nboun,ilboun)
         endif
      endif
!
      return
      end

