!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine objective_freq(dgdx,
     &  df,vold,ndesi,iobject,mi,nactdofinv,
     &  jqs,irows)
!
      implicit none
!
      integer ndesi,iobject,mi(*),idesvar,j,idir,
     &  jqs(*),irows(*),nactdofinv(*),node,idof,inode,mt
!      
      real*8 dgdx(ndesi,*),df(*),vold(0:mi(2),*)
!
      intent(in) df,vold,ndesi,iobject,mi,nactdofinv,jqs,irows
!
      intent(inout) dgdx
!
!     ----------------------------------------------------------------
!     Calculation of the total differential:
!     dgdx = dgdx + vold^(T) * ( df )
!     ----------------------------------------------------------------
!     
      mt=mi(2)+1
!
      do idesvar=1,ndesi
         do j=jqs(idesvar),jqs(idesvar+1)-1
            idof=irows(j)
            inode=nactdofinv(idof)
            node=inode/mt+1
            idir=inode-mt*(inode/mt)
            dgdx(idesvar,iobject)=dgdx(idesvar,iobject)
     &                           +vold(idir,node)*df(j)
         enddo
      enddo
!      
      return
      end
