!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine calceigenvalues(c,al)
!
!     calculates the eigenvalues al of the symmetric 3x3 matrix c
!     the eigenvalues are sorted in increasing order
!
      implicit none
!
      integer idummy(3),three,kflag,i
!
      real*8 c(3,3),al(3),v1,v2,v3,bb,cc,cm,cn,tt,pi
!     
      data kflag /1/
      data three /3/
!
!     calculation of the eigenvalues of c
!     Simo & Hughes Computational Inelasticity p 244
!
      pi=4.d0*datan(1.d0)
!
      v1=c(1,1)+c(2,2)+c(3,3)
      v2=c(2,2)*c(3,3)+c(1,1)*c(3,3)+c(1,1)*c(2,2)-
     &     (c(2,3)*c(2,3)+c(1,3)*c(1,3)+c(1,2)*c(1,2))
      v3=c(1,1)*(c(2,2)*c(3,3)-c(2,3)*c(2,3))
     &     -c(1,2)*(c(1,2)*c(3,3)-c(1,3)*c(2,3))
     &     +c(1,3)*(c(1,2)*c(2,3)-c(1,3)*c(2,2))
!
      bb=v2-v1*v1/3.d0
      cc=-2.d0*v1**3/27.d0+v1*v2/3.d0-v3
      if(dabs(bb).le.1.d-10) then
         if(dabs(cc).gt.1.d-10) then
            al(1)=-cc**(1.d0/3.d0)
         else
            al(1)=0.d0
         endif
         al(2)=al(1)
         al(3)=al(1)
      else
         cm=2.d0*dsqrt(-bb/3.d0)
         cn=3.d0*cc/(cm*bb)
         if(dabs(cn).gt.1.d0) then
            if(cn.gt.1.d0) then
               cn=1.d0
            else
               cn=-1.d0
            endif
         endif
         tt=datan2(dsqrt(1.d0-cn*cn),cn)/3.d0
         al(1)=cm*dcos(tt)
         al(2)=cm*dcos(tt+2.d0*pi/3.d0)
         al(3)=cm*dcos(tt+4.d0*pi/3.d0)
      endif
      do i=1,3
         al(i)=al(i)+v1/3.d0
      enddo
!     
!     sorting
!     
      call dsort(al,idummy,three,kflag)
!
      return
      end
