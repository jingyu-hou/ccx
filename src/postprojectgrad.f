!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine postprojectgrad(ndesi,nodedesi,dgdxglob,nactive,
     &   nobject,nnlconst,ipoacti,nk,iconst,objectset,iconstacti,
     &   inameacti)         
!
!     calculates the projected gradient
!
      implicit none
!
      character*81 objectset(4,*)
!
      integer ndesi,nodedesi(*),irow,icol,nactive,nobject,nnlconst,
     &   ipoacti(*),nk,ipos,node,iconst,i,m,iconstacti(*),
     &   inameacti(*)
!
      real*8 dgdxglob(2,nk,nobject),scalar,dd,len
!
!     calculation of final projected gradient
!
      if(nactive.gt.0) then
         do irow=1,ndesi
            node=nodedesi(irow)
            dgdxglob(2,node,nobject)=dgdxglob(2,node,1)
     &          -dgdxglob(2,node,nobject)
         enddo
         objectset(1,nobject)(1:11)='PROJECTGRAD'        
!
         write(*,*)
         write(*,*) '*INFO in postprojectgrad:'
         write(*,*) '      at least 1 constraint active:'
         write(*,*) '      projected gradient has been '
         write(*,*) '      calculated based on the '
         write(*,*) '      constraints:'
         if(nnlconst.eq.nactive) then
            do i=1,nactive
               write(*,'(7x,a12)') objectset(1,ipoacti(i))
            enddo
            write(*,*)
         elseif((nnlconst.lt.nactive).and.(nnlconst.gt.0)) then
            do i=1,nnlconst
               write(*,'(7x,a12)') objectset(1,ipoacti(i))
            enddo
            write(*,'(7x,a12)') objectset(1,inameacti(i))
            write(*,*)
         elseif(nnlconst.eq.0) then
            write(*,'(7x,a12)') objectset(1,inameacti(i))
            write(*,*)
         endif
      else
         write(*,*)
         write(*,*) '*INFO in postprojectgrad:'
         write(*,*) '      no constraint active:'
         write(*,*) '      no projected gradient '
         write(*,*) '      calculated' 
         write(*,*)
      endif    
!
      return        
      end

