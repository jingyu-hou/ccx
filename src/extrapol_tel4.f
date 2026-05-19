!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine extrapol_tel4(ielfa,vfa,vfap,gradtfa,rf,ifabou,
     &  ipnei,vel,xxi,xle,nef,nfacea,nfaceb)
!
!     extrapolation of temperature values to the faces (taking the
!     skewness of the elements into account)
!
!     for temperature calculations the external faces have
!     as boundary condition either a specified temperature or
!     a specified flux. If the user did not apply any of these,
!     a zero specified flux is implicitly assumed
!
      implicit none
!
      integer ielfa(4,*),ifabou(*),nfacea,nfaceb,nef,i,iel1,iel2,
     &  indexf,ipnei(*)
!
      real*8 vfap(0:7,*),vel(nef,0:7),
     &  vfa(0:7,*),rf(3,*),gradtfa(3,*),xle(*),xxi(3,*)
!
!
!
!     Moukalled et al. p 279
!
      do i=nfacea,nfaceb
         iel1=ielfa(1,i)
         iel2=ielfa(2,i)
         if(iel2.gt.0) then
!     
!     interpolation
!     
            vfa(0,i)=vfap(0,i)+gradtfa(1,i)*rf(1,i)
     &           +gradtfa(2,i)*rf(2,i)
     &           +gradtfa(3,i)*rf(3,i)
         elseif(ielfa(3,i).gt.0) then
!     
!     no implicit zero gradient
!     
            if(ifabou(-iel2).gt.0) then
!     
!     temperature given
!     
               vfa(0,i)=vfap(0,i)
            else
!     
!     flux given: gradient=flux
!     
               indexf=ipnei(iel1)+ielfa(4,i)
               vfa(0,i)=vel(iel1,0)
     &              +(gradtfa(1,i)*xxi(1,indexf)+
     &              gradtfa(2,i)*xxi(2,indexf)+
     &              gradtfa(3,i)*xxi(3,indexf))*xle(indexf)
            endif
         else
!     
!     zero gradient
!     
            vfa(0,i)=vel(iel1,0)
         endif
      enddo
!            
      return
      end
