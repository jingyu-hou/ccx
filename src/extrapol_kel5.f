!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine extrapol_kel5(ielfa,ipnei,vel,xlet,gradkfa,xxj,
     &  nef,nfacea,nfaceb,ncfd)
!
!     correct the facial turbulent kinetic energy gradients:
!     Moukalled et al. p 289
!
      implicit none
!
      integer ielfa(4,*),ipnei(*),nef,nfacea,nfaceb,i,k,iel1,iel2,
     &  indexf,ncfd
!
      real*8 vel(nef,0:7),xlet(*),gradkfa(3,*),xxj(3,*),dd
!
!
!
      do i=nfacea,nfaceb
         iel2=ielfa(2,i)
         if(iel2.gt.0) then
            iel1=ielfa(1,i)
            indexf=ipnei(iel1)+ielfa(4,i)
            dd=(vel(iel2,6)-vel(iel1,6))/xlet(indexf)
     &        -gradkfa(1,i)*xxj(1,indexf)
     &        -gradkfa(2,i)*xxj(2,indexf)
     &        -gradkfa(3,i)*xxj(3,indexf)
            do k=1,ncfd
               gradkfa(k,i)=gradkfa(k,i)+dd*xxj(k,indexf)
            enddo
         endif
      enddo
!            
      return
      end
