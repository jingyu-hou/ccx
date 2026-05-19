!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine extrapol_pel3(ielfa,xrlfa,icyclic,ifatie,gradpfa,
     &  gradpel,c,ipnei,xxi,nfacea,nfaceb,ncfd)
!
!     interpolate/extrapolate the pressure gradient from the
!     center of the elements to the center of the faces
!
      implicit none
!
      integer ielfa(4,*),icyclic,ifatie(*),ipnei(*),nfacea,nfaceb,
     &  iel1,iel2,i,l,indexf,ncfd
!
      real*8 xrlfa(3,*),gradpfa(3,*),gradpel(3,*),c(3,3),xxi(3,*),
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
!           face in between two elements
!
            xl2=xrlfa(2,i)
            if((icyclic.eq.0).or.(ifatie(i).eq.0)) then
               do l=1,ncfd
                  gradpfa(l,i)=xl1*gradpel(l,iel1)+
     &                 xl2*gradpel(l,iel2)
               enddo
            elseif(ifatie(i).gt.0) then
               do l=1,ncfd
                  gradpfa(l,i)=xl1*gradpel(l,iel1)+xl2*
     &                  (gradpel(1,iel2)*c(l,1)+
     &                   gradpel(2,iel2)*c(l,2)+
     &                   gradpel(3,iel2)*c(l,3))
               enddo
            else
               do l=1,ncfd
                  gradpfa(l,i)=xl1*gradpel(l,iel1)+xl2*
     &                  (gradpel(1,iel2)*c(1,l)+
     &                   gradpel(2,iel2)*c(2,l)+
     &                   gradpel(3,iel2)*c(3,l))
               enddo
            endif
         elseif(ielfa(3,i).ne.0) then
!     
!           boundary face; more than one layer; extrapolation
!     
            do l=1,ncfd
               gradpfa(l,i)=xl1*gradpel(l,iel1)+
     &              xrlfa(3,i)*gradpel(l,abs(ielfa(3,i)))
            enddo
         else
!     
!           boundary face; one layer
!   
            indexf=ipnei(iel1)+ielfa(4,i)
            gradnor=gradpel(1,iel1)*xxi(1,indexf)+
     &              gradpel(2,iel1)*xxi(2,indexf)+
     &              gradpel(3,iel1)*xxi(3,indexf)
            do l=1,ncfd
                  gradpfa(l,i)=gradpel(l,iel1)
     &                        -gradnor*xxi(l,indexf)
            enddo
         endif
      enddo
!            
      return
      end
