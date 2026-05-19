!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine scalesen(dgdxglob,nobject,nk,nodedesi,ndesi,
     &   objectset,iscaleflag)
!
!     Scaling the sensitivities      
!
!     iscaleflag=0: greatest vector value is scaled to 1
!     iscaleflag=1: length of the vector is scaled to 1     
!
      implicit none
!
      character*81 objectset(4,*)
!
      integer nobject,nk,nodedesi(*),i,ndesi,m,iscaleflag,kflag,node
!
      real*8 dgdxglob(2,nk,nobject),dd,len
!
!
      kflag=0
      if(iscaleflag.eq.0) then
         if(objectset(1,nobject)(1:11).eq.'PROJECTGRAD') then
            kflag=1
         endif
!     
         len=0.d0
         dd=0.d0
         if(kflag.eq.1) then
            do i=1,ndesi
               node=nodedesi(i)
               len=len+dgdxglob(2,node,nobject)**2
               dd=max(dd,abs(dgdxglob(1,node,nobject)))
            enddo
            if(dd.ne.0.d0) then
               do i=1,ndesi
                  node=nodedesi(i)
                  dgdxglob(1,node,nobject)=dgdxglob(1,node,nobject)/dd
               enddo
            endif
         else
            do i=1,ndesi
               node=nodedesi(i)
               len=len+dgdxglob(2,node,1)**2
            enddo
         endif
         len=dsqrt(len)
c         write(5,*)
c         write(5,*) 'LENGTH OF DESCENT GRADIENT VECTOR:'
c         write(5,*)
c         write(5,'(7x,e14.7)') len
         do m=1,nobject
            if(objectset(1,m)(1:9).eq.'THICKNESS') cycle
            dd=0.d0
            do i=1,ndesi
               node=nodedesi(i)
               dd=max(dd,abs(dgdxglob(2,node,m)))
            enddo
            do i=1,ndesi
               node=nodedesi(i)
               dgdxglob(2,node,m)=dgdxglob(2,node,m)/dd
               if(objectset(1,m)(1:11).eq.'PROJECTGRAD') then
                  dgdxglob(1,node,m)=dgdxglob(1,node,m)/dd
               endif
            enddo
         enddo
      elseif(iscaleflag.eq.1) then
         do m=1,nobject
            if(objectset(1,m)(1:9).eq.'THICKNESS') cycle
            if(objectset(1,m)(1:9).eq.'FIXGROWTH') cycle
            if(objectset(1,m)(1:12).eq.'FIXSHRINKAGE') cycle
            dd=0.d0
            do i=1,ndesi
               node=nodedesi(i)
               dd=dd+dgdxglob(2,node,m)**2
            enddo
            dd=dsqrt(dd)
            do i=1,ndesi
               node=nodedesi(i)
               dgdxglob(2,node,m)=dgdxglob(2,node,m)/dd
            enddo
         enddo
      endif
!     
      return        
      end
      



