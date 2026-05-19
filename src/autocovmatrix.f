!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine autocovmatrix(co,ad,au,jqs,irows,ndesi,nodedesi,
     &  physcon)         
!
!     calculates the values of the autocovariance matrix
!
      implicit none
!
      integer jqs(*),irows(*),ndesi,nodedesi(*),idof,j,jdof,node1,
     &  node2
!
      real*8 co(3,*),ad(*),au(*),physcon(*),dist,corrlength,sigma
!
      corrlength=physcon(13)
      sigma=physcon(12)
!
      do idof=1,ndesi
         ad(idof)=sigma*sigma
         do j=jqs(idof),jqs(idof+1)-1
            jdof=irows(j)
            node1=nodedesi(idof)
            node2=nodedesi(jdof)
            dist=dsqrt((co(1,node1)-co(1,node2))**2+
     &                 (co(2,node1)-co(2,node2))**2+
     &                 (co(3,node1)-co(3,node2))**2)
!
!           assign the value to the autocovariance matrix
!
            au(j)=ad(idof)*dexp(-(dist/corrlength)**2)
         enddo
      enddo
!
      return        
      end




