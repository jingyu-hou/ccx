!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine calcflux(area,vfa,xxna,ipnei,nef,neifa,flux,xxj,
     &  gradpfa,xlet,xle,vel,advfa,ielfa,neiel,ifabou,hfa,nefa,
     &  nefb)
!
!     calculate the mass flow using the newly calculated velocity
!     and based on the Rhie-Chow interpolation
!
!     vfa(1..3,*) is only used for external faces
!
      implicit none
!
      integer i,indexf,ipnei(*),ifa,nef,neifa(*),ielfa(4,*),
     &  neiel(*),iwall,iel,ipointer,ifabou(*),nefa,nefb
!
      real*8 area(*),vfa(0:7,*),xxna(3,*),flux(*),xxj(3,*),gradpfa(3,*),
     &  xlet(*),xle(*),vel(nef,0:7),advfa(*),coef,hfa(3,*)
!
!
!
      do i=nefa,nefb
         do indexf=ipnei(i)+1,ipnei(i+1)
            ifa=neifa(indexf)
            iel=neiel(indexf)
!
            if(iel.eq.0) then
!
!              external face
!
               iwall=0
               ipointer=abs(ielfa(2,ifa))
               if(ipointer.gt.0) then
                  iwall=ifabou(ipointer+5)
               endif
!     
               if(iwall.eq.0) then
!
!                 no wall nor sliding conditions
!
                  if(ielfa(3,ifa).gt.0) then
!
!                    inlet: velocity fixed
!                     
                     flux(indexf)=vfa(5,ifa)*
     &                            (vfa(1,ifa)*xxna(1,indexf)+
     &                             vfa(2,ifa)*xxna(2,indexf)+
     &                             vfa(3,ifa)*xxna(3,indexf))
                  else
!
!                    outlet or 1-layer
!
                     coef=advfa(ifa)*
     &                    ((vfa(4,ifa)-vel(i,4))/xle(indexf)-
     &                     (gradpfa(1,ifa)*xxj(1,indexf)+
     &                      gradpfa(2,ifa)*xxj(2,indexf)+
     &                      gradpfa(3,ifa)*xxj(3,indexf)))
                     flux(indexf)=vfa(5,ifa)*
     &                  ((hfa(1,ifa)-coef*xxj(1,indexf))*xxna(1,indexf)+
     &                   (hfa(2,ifa)-coef*xxj(2,indexf))*xxna(2,indexf)+
     &                   (hfa(3,ifa)-coef*xxj(3,indexf))*xxna(3,indexf))
                  endif
               else
!
!                 wall or sliding: velocity of the wall
!
                  flux(indexf)=0.d0
c                  flux(indexf)=vfa(5,ifa)*
c     &                 (vfa(1,ifa)*xxna(1,indexf)+
c     &                  vfa(2,ifa)*xxna(2,indexf)+
c     &                  vfa(3,ifa)*xxna(3,indexf))
               endif
            else
!
!              internal face
!
               coef=advfa(ifa)*
     &              ((vel(iel,4)-vel(i,4))/xlet(indexf)-
     &              (gradpfa(1,ifa)*xxj(1,indexf)+
     &              gradpfa(2,ifa)*xxj(2,indexf)+
     &              gradpfa(3,ifa)*xxj(3,indexf)))
               flux(indexf)=vfa(5,ifa)*
     &              ((hfa(1,ifa)-coef*xxj(1,indexf))*xxna(1,indexf)+
     &              (hfa(2,ifa)-coef*xxj(2,indexf))*xxna(2,indexf)+
     &              (hfa(3,ifa)-coef*xxj(3,indexf))*xxna(3,indexf))
             endif
         enddo
      enddo
!  
      return
      end
