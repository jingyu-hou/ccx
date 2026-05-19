!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine objective_freq_cs(dgdx,
     &  df,vold,ndesi,iobject,mi,nactdofinv,
     &  jqs,irows,nk,nzss)
!
      implicit none
!
      integer ndesi,iobject,mi(*),idesvar,j,idir,nk,nzss,
     &  jqs(*),irows(*),nactdofinv(*),node,idof,inode,mt
!      
      real*8 dgdx(ndesi,*),df(*),vold(0:mi(2),*)
!
      intent(in) df,vold,ndesi,iobject,mi,nactdofinv,jqs,irows,nk,nzss
!
      intent(inout) dgdx
!
!     ----------------------------------------------------------------
!     Calculation of the frequency sensitivity in the cyclic symmetric case
!
!     dgdx = -vold_R^T*df_R+vold_I^T*df_I
!
!     notice that df was calculated in mafillsmcsse as 
!     - (K-lambda*M)*eigenvector
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
     &          -vold(idir,node)*df(j)-vold(idir,nk+node)*df(nzss+j)
         enddo
      enddo
!      
      return
      end
