!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine extrapol_pel4(ielfa,vfa,vfap,gradpfa,rf,ifabou,
     &  ipnei,vel,xxi,xle,nef,nfacea,nfaceb)
!
!     inter/extrapolation of p at the center of the elements
!     to the center of the faces: taking the skewness of the 
!     elements into account (using the gradient at the center
!     of the faces)
!
      implicit none
!
      integer ielfa(4,*),ifabou(*),ipnei(*),nef,nfacea,nfaceb,i,indexf,
     &  iel1,iel2,ibou
!
      real*8 vfa(0:7,*),vfap(0:7,*),gradpfa(3,*),rf(3,*),vel(nef,0:7),
     &  xxi(3,*),xle(*)
!
!
!
!        Moukalled et al. p 279
!
      do i=nfacea,nfaceb
         iel1=ielfa(1,i)
         iel2=ielfa(2,i)
         if(iel2.gt.0) then
!
!           face between two elements
!
            vfa(4,i)=vfap(4,i)+gradpfa(1,i)*rf(1,i)
     &           +gradpfa(2,i)*rf(2,i)
     &           +gradpfa(3,i)*rf(3,i)
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
!              pressure given
!     
               vfa(4,i)=vfap(4,i)
            else
!     
!              extrapolation
!     
               vfa(4,i)=vfap(4,i)+gradpfa(1,i)*rf(1,i)
     &              +gradpfa(2,i)*rf(2,i)
     &              +gradpfa(3,i)*rf(3,i)
            endif
         else
!     
!           boundary face; one layer
!     
c            indexf=ipnei(iel1)+ielfa(4,i)
            vfa(4,i)=vel(iel1,4)
c     &           +(gradpfa(1,i)*xxi(1,indexf)+
c     &           gradpfa(2,i)*xxi(2,indexf)+
c     &           gradpfa(3,i)*xxi(3,indexf))*xle(indexf)
         endif
      enddo
!            
      return
      end
