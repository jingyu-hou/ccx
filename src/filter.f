!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine filter(dgdxglob,nobject,nk,nodedesi,ndesi,
     &           objectset,xo,yo,zo,x,y,z,nx,ny,nz,neighbor,r,
     &           ndesia,ndesib,xdesi,distmin)               
!
!     Filtering of sensitivities      
!
      implicit none
!
      character*81 objectset(4,*)

      integer nobject,nk,nodedesi(*),nnodesinside,i,actdir,
     &        ndesi,j,m,neighbor(ndesi+6),nx(ndesi),
     &        ny(ndesi),nz(ndesi),istat,ndesia,ndesib
!
      real*8 dgdxglob(2,nk,nobject),xo(ndesi),yo(ndesi),zo(ndesi),
     &       x(ndesi),y(ndesi),z(ndesi),filterrad,r(ndesi+6),
     &       filterval(ndesi),nominator,denominator,distmin,
     &       xdesi(3,ndesi),scalar
!
!     Calculate filtered sensitivities
!
!     Check if direction weighting is turned on
!
      if(objectset(2,1)(14:16).eq.'DIR') then
         actdir=1
      else
         actdir=0
      endif
!   
!     Assign filter radius (taken from first defined object function)
!
      read(objectset(2,1)(21:40),'(f20.0)',iostat=istat) filterrad     
!        
      do j=ndesia,ndesib
!     
         call near3d_se(xo,yo,zo,x,y,z,nx,ny,nz,xo(j),yo(j),zo(j),
     &        ndesi,neighbor,r,nnodesinside,filterrad)
!  
!        Calculate function value of the filterfunction CONST 
!
         if(objectset(2,1)(1:5).eq.'CONST') then
            do i=1,nnodesinside
               filterval(i)=1.0d0
            enddo
!  
!        Calculate function value of the filterfunction LIN 
!
         elseif(objectset(2,1)(1:3).eq.'LIN') then
            do i=1,nnodesinside
               filterval(i)=(1-dsqrt(r(i))/filterrad)*filterrad
            enddo
!
!        Calculate function value of the filterfunction QUAD
!
         elseif(objectset(2,1)(1:4).eq.'QUAD') then
            do i=1,nnodesinside
               filterval(i)=(-1*(1+dsqrt(r(i))/filterrad)
     &                      *(-1+dsqrt(r(i))/filterrad))*filterrad
            enddo
!
!        Calculate function value of the filterfunction CUB
!
         elseif(objectset(2,1)(1:5).eq.'CUBIC') then
            do i=1,nnodesinside
               filterval(i)=(2*(dsqrt(r(i))/filterrad)**3
     &                      -3*(dsqrt(r(i))/filterrad)**2+1)*filterrad
            enddo 
         endif
!  
!        Calculate filtered sensitivity
! 
         do m=1,nobject
            if(objectset(1,m)(1:9).eq.'THICKNESS') cycle             
            if(objectset(1,m)(1:9).eq.'FIXGROWTH') cycle             
            if(objectset(1,m)(1:12).eq.'FIXSHRINKAGE') cycle
            nominator=0.d0
            denominator=0.d0
            do i=1,nnodesinside
               if(actdir.eq.1) then          
                  scalar=(xdesi(1,j)*xdesi(1,neighbor(i))
     &                +xdesi(2,j)*xdesi(2,neighbor(i))
     &                +xdesi(3,j)*xdesi(3,neighbor(i)))/(distmin**2)
                  if(scalar.lt.0.d0) then
                     scalar=0.d0
                  endif
                  nominator=nominator+filterval(i)*scalar*
     &                dgdxglob(1,nodedesi(neighbor(i)),m)
                  denominator=denominator+filterval(i)      
               else
                  nominator=nominator+filterval(i)*
     &                dgdxglob(1,nodedesi(neighbor(i)),m)
                  denominator=denominator+filterval(i)
               endif
            enddo 
            dgdxglob(2,nodedesi(j),m)=nominator/denominator
         enddo
      enddo
!
      return        
      end




