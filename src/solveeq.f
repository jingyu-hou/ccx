!
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine solveeq(adb,aub,adl,addiv,b,sol,aux,icol,irow,jq,
     &  neq,nzs,nzl)
!
!     solving a system of equations by iteratively solving the
!     lumped version
!     The diagonal terms f the original system are stored in adb,
!     the off-diagonal terms in aub
!     Ref: The Finite Element Method for Fluid Dynamics,
!          O.C. Zienkiewicz, R.L. Taylor & P. Nithiarasu
!          6th edition (2006) ISBN 0 7506 6322 7
!          p. 61
!
      implicit none
!
      integer icol(*),irow(*),jq(*),neq,nzs,nzl,i,j,k,maxit
!
      real*8 adb(*),aub(*),adl(*),addiv(*),b(*),sol(*),aux(*),p
!
      data maxit /1/
!
!     first iteration
!
      do i=1,neq
         sol(i)=b(i)*adl(i)
c         write(*,*) 'solveeq ',i,b(i),adl(i)
      enddo
      if(maxit.eq.1) return
!
!     iterating maxit times
!
      do k=2,maxit
!
!        multiplying the difference of the original matrix
!        with the lumped matrix with the actual solution 
!
c         call opfem(neq,p,sol,aux,adb,aub,icol,irow,nzl)
         call op(neq,sol,aux,adb,aub,jq,irow)
!
         do i=1,neq
            sol(i)=(b(i)-aux(i))*adl(i)
         enddo
!
      enddo
!
      return
      end
