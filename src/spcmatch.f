!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine spcmatch(xboun,nodeboun,ndirboun,nboun,xbounold,
     &   nodebounold,ndirbounold,nbounold,ikboun,ilboun,vold,reorder,
     &   nreorder,mi)
!
!     matches SPC boundary conditions of one step with those of
!     the previous step
!
      implicit none
!
      integer nodeboun(*),ndirboun(*),nboun,nodebounold(*),ilboun(*),
     &  ndirbounold(*),nbounold,i,kflag,idof,id,nreorder(*),ikboun(*),
     &  mi(*)
!
      real*8 xboun(*),xbounold(*),vold(0:mi(2),*),reorder(*)
!
      kflag=2
!
      do i=1,nboun
         nreorder(i)=0
      enddo
!
      do i=1,nbounold
         idof=8*(nodebounold(i)-1)+ndirbounold(i)
         if(nboun.gt.0) then
            call nident(ikboun,idof,nboun,id)
         else
            id=0
         endif
         if((id.gt.0).and.(ikboun(id).eq.idof)) then
            reorder(ilboun(id))=xbounold(i)
            nreorder(ilboun(id))=1
         endif
      enddo
!
      do i=1,nboun
         if(nreorder(i).eq.0) then
            if(ndirboun(i).gt.4) then
               reorder(i)=0.d0
            else
               reorder(i)=vold(ndirboun(i),nodeboun(i))
            endif
         endif
      enddo
!
      do i=1,nboun
         xbounold(i)=reorder(i)
      enddo
!
      return
      end

