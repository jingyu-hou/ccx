      subroutine bmatrix(shp,xkl,nope,b)
!      subroutine bmatrix(shp,xkl,nope,iint,b,e,c)
      implicit none

      real*8 shp(4,26),xkl(3,3),ba(6,3),b(6,24),
     &       id(6,6),m(6,1),n(1,8),e(8,8),c(24,8),weight
      integer nope,iint,i,j,k,inode
      include "gauss.f"
      
      m=reshape((/1.d0,1.d0,1.d0,0.d0,0.d0,0.d0/),(/6,1/))
      
!      weight=weight3d2(iint)
      
      do inode=1,nope
         do i=1,3
            do j=1,3
               ba(i,j)=xkl(j,i)*shp(i,inode)
               k=i+1
               if (k.gt.3) then
                  k=k-3
               endif
               ba(i+3,j)=xkl(j,i)*shp(k,inode)+xkl(j,k)*shp(i,inode)
            enddo
         enddo
      
         do i=1,6
            do j=1,3
               b(i,(inode-1)*3+j)=ba(i,j)
            enddo
         enddo
      enddo

!      id=0.d0
!      do i=1,6
!         id(i,i)=1.d0
!      enddo
!
!      id=id-1.d0/3.d0*matmul(m,transpose(m))
!
!      do i=1,nope
!         n(1,i)=shp(4,i)
!         n(1,i)=1.d0/8.d0
!      enddo
!      
!      e=e+matmul(transpose(n),n)*weight/8.d0
!      c=c+matmul(transpose(bb),matmul(m,n))*weight/8.d0
      

      return
      end subroutine
