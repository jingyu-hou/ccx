!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine extrapol_kel3(ielfa,xrlfa,icyclic,ifatie,gradkfa,
     &  gradkel,c,ipnei,xxi,nfacea,nfaceb,ncfd)
!
!     interpolate/extrapolate the turbulent kinetic energy gradient
!     from the center of the elements to the center of the faces
!           
      implicit none
!
      integer ielfa(4,*),icyclic,ifatie(*),ipnei(*),nfacea,nfaceb,
     &  iel1,iel2,i,l,indexf,ncfd
!
      real*8 xrlfa(3,*),gradkfa(3,*),gradkel(3,*),c(3,3),xxi(3,*),
     &  gradnor,xl1,xl2
!
!
! 
      do i=nfacea,nfaceb
         iel1=ielfa(1,i)
         xl1=xrlfa(1,i)
         iel2=ielfa(2,i)
         if(iel2.gt.0) then
!
!           face between two elements
!
            xl2=xrlfa(2,i)
            if((icyclic.eq.0).or.(ifatie(i).eq.0)) then
               do l=1,ncfd
                  gradkfa(l,i)=xl1*gradkel(l,iel1)+
     &                 xl2*gradkel(l,iel2)
               enddo
            elseif(ifatie(i).gt.0) then
               do l=1,ncfd
                  gradkfa(l,i)=xl1*gradkel(l,iel1)+xl2*
     &                  (gradkel(1,iel2)*c(l,1)+
     &                   gradkel(2,iel2)*c(l,2)+
     &                   gradkel(3,iel2)*c(l,3))
               enddo
            else
               do l=1,ncfd
                  gradkfa(l,i)=xl1*gradkel(l,iel1)+xl2*
     &                  (gradkel(1,iel2)*c(1,l)+
     &                   gradkel(2,iel2)*c(2,l)+
     &                   gradkel(3,iel2)*c(3,l))
               enddo
            endif
         elseif(ielfa(3,i).gt.0) then
!
!           the above condition implies iel2!=0
!           if iel2 were zero, no b.c. would apply and
!           ielfa(3,i) would be zero or negative
!     
!           boundary face; no zero gradient
!
            do l=1,ncfd
               gradkfa(l,i)=xl1*gradkel(l,iel1)+
     &              xrlfa(3,i)*gradkel(l,abs(ielfa(3,i)))
            enddo
         else
!     
!           boundary face; zero gradient
!   
            indexf=ipnei(iel1)+ielfa(4,i)
            gradnor=gradkel(1,iel1)*xxi(1,indexf)+
     &              gradkel(2,iel1)*xxi(2,indexf)+
     &              gradkel(3,iel1)*xxi(3,indexf)
            do l=1,ncfd
                  gradkfa(l,i)=gradkel(l,iel1)
     &                        -gradnor*xxi(l,indexf)
            enddo
         endif
      enddo
!            
      return
      end
