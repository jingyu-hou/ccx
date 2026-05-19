!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine objective_stress_tot(dgdx,df,ndesi,iobject,jqs,
     &   irows,dgdu)
!
      implicit none
!
      integer ndesi,iobject,idesvar,j,jqs(*),irows(*),idof
!      
      real*8 dgdx(ndesi,*),df(*),dgdu(*)
!
      intent(inout) dgdu,dgdx
!
!     ----------------------------------------------------------------
!     Calculation of the total differential:
!     non-linear:  dgdx = dgdx + dgdu * ( df )
!     ----------------------------------------------------------------
!
!     Calculation of the total differential:    
!
      do idesvar=1,ndesi
         do j=jqs(idesvar),jqs(idesvar+1)-1
            idof=irows(j)
            dgdx(idesvar,iobject)=dgdx(idesvar,iobject) 
     &           +dgdu(idof)*df(j)
         enddo
      enddo
!      
      return
      end
