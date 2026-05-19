!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine extrapol_vel5(ielfa,ipnei,vel,xlet,gradvfa,xxj,
     &  nef,nfacea,nfaceb,ncfd)
!
!     correct the facial velocity gradients:
!     Moukalled et al. p 289
!
      implicit none
!
      integer ielfa(4,*),ipnei(*),nef,nfacea,nfaceb,i,k,l,indexf,iel1,
     &  iel2,ncfd
!
      real*8 vel(nef,0:7),xlet(*),gradvfa(3,3,*),xxj(3,*),dd
!
!
!
      do i=nfacea,nfaceb
         iel2=ielfa(2,i)
         if(iel2.gt.0) then
            iel1=ielfa(1,i)
            indexf=ipnei(iel1)+ielfa(4,i)
            do k=1,ncfd
               dd=(vel(iel2,k)-vel(iel1,k))/xlet(indexf)
     &              -gradvfa(k,1,i)*xxj(1,indexf)
     &              -gradvfa(k,2,i)*xxj(2,indexf)
     &              -gradvfa(k,3,i)*xxj(3,indexf)
               do l=1,ncfd
                  gradvfa(k,l,i)=gradvfa(k,l,i)+dd*xxj(l,indexf)
               enddo
            enddo
         endif
      enddo
!
      return
      end
