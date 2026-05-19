!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine extrapol_tel2(ipnei,neifa,vfa,area,xxna,volume,gradtel,
     &  nefa,nefb,ncfd)
!
!     calculate the gradient of the temperature at the center of
!     the elements
!
      implicit none
!
      integer ipnei(*),neifa(*),nefa,nefb,ifa,i,l,indexf,ncfd
!
      real*8 vfa(0:7,*),area(*),xxna(3,*),volume(*),gradtel(3,*)
!
!
!
      do i=nefa,nefb
!
!        initialization
!     
         do l=1,ncfd
            gradtel(l,i)=0.d0
         enddo
!
         do indexf=ipnei(i)+1,ipnei(i+1)
            ifa=neifa(indexf)
            do l=1,ncfd
               gradtel(l,i)=gradtel(l,i)+
     &              vfa(0,ifa)*xxna(l,indexf)
            enddo
         enddo
!     
!        dividing by the volume of the element
!     
         do l=1,ncfd
            gradtel(l,i)=gradtel(l,i)/volume(i)
         enddo
      enddo
!            
      return
      end
