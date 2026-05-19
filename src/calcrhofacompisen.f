!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine calcrhofacompisen(nface,vfa,shcon,ielmat,ntmat_,
     &  mi,ielfa,cvfa,velo,nef)
!
!     calculation of the density at the face centers
!     (compressible fluids)
!
      implicit none
!
      integer nface,i,imat,ntmat_,mi(*),nef,velo(nef,0:7),
     &  ielmat(mi(3),*),ielfa(4,*),j
!
      real*8 t1l,vfa(0:7,*),shcon(0:3,ntmat_,*),cvfa(*) 
!     
      do i=1,nface
         t1l=vfa(0,i)
         j=ielfa(1,i)
!
!        take the material of the first adjacent element
!
         imat=ielmat(1,ielfa(1,i))
         vfa(5,i)=vfa(4,i)/(shcon(3,1,imat)*
c     &       (10.5d0-(vfa(1,i)**2+vfa(2,i)**2+vfa(3,i)**2)/2.d0))
     &       (5.98696d0-(vfa(1,i)**2+vfa(2,i)**2+vfa(3,i)**2)/2.d0))
c     &       (1.41827d0-(vfa(1,i)**2+vfa(2,i)**2+vfa(3,i)**2)/2.d0))
c         vfa(5,i)=vfa(4,i)/(shcon(3,1,imat)*
c     &       (velo(j,0)+(velo(j,1)**2+velo(j,2)**2+velo(j,3)**2)/2.d0-
c     &(vfa(1,i)**2+vfa(2,i)**2+vfa(3,i)**2)/2.d0))
      enddo
!            
      return
      end
