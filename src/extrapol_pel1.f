!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine extrapol_pel1(ielfa,xrlfa,vfap,vel,ifabou,xbounact,
     &  nef,nfacea,nfaceb)
!
!     extrapolation of pressure element values to the faces
!
      implicit none
!
      integer ielfa(4,*),ifabou(*),nfacea,nfaceb,nef,i,iel1,iel2,
     &  ibou
!
      real*8 xrlfa(3,*),vfap(0:7,*),vel(nef,0:7),xbounact(*),xl1
!
!
!     
      do i=nfacea,nfaceb
         iel1=ielfa(1,i)
         xl1=xrlfa(1,i)
         iel2=ielfa(2,i)
         if(iel2.gt.0) then
!
!           face between two elements: interpolation
!
            vfap(4,i)=xl1*vel(iel1,4)+xrlfa(2,i)*vel(iel2,4)
!
        elseif(ielfa(3,i).ne.0) then
!
!           boundary face; more than one layer
!            
            ibou=0
            if(iel2.lt.0) then
               if(ifabou(-iel2+4).gt.0) then
                  ibou=ifabou(-iel2+4)
               endif
            endif
!
            if(ibou.gt.0) then
!
!              pressure boundary condition
!
               vfap(4,i)=xbounact(ibou)
            else
!
!              extrapolation
!
               vfap(4,i)=xl1*vel(iel1,4)
     &                 +xrlfa(3,i)*vel(abs(ielfa(3,i)),4)
           endif
         else
!
!           boundary face; one layer
!
            vfap(4,i)=vel(iel1,4)
         endif
      enddo
!            
      return
      end
