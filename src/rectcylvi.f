!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine rectcylvi(co,v,fn,stn,qfn,een,cs,n,icntrl,t,filab,
     &  imag,mi,emn)
!
!     cf. subroutine rectcyl
!     In the present routine, the imaginary part of the 
!     displacements and stresses for all sectors are transformed 
!     from the cylindrical into the rectangular system
!
      implicit none
!
      character*87 filab(*)
      integer i,j,n,icntrl,imag,mi(*)
      real*8 co(3,*),v(0:mi(2),*),fn(0:mi(2),*),stn(6,*),een(6,*),
     &  a(3,3),emn(6,*),
     &  xr,xt,xz,b(3,3),cs(17,*),t(3),u(3),qfn(3,*),csab(7),
     &  xn(3),r(3),z,theta,rr,c(3,3),ctm,ct,st,ddx,ddy,dd
!
      do i=1,7
         csab(i)=cs(5+i,1)
      enddo
!
      do i=1,n
         j=i
         call transformatrix(csab,co(1,i),a)
!     
         if((filab(1)(1:3).eq.'U  ').or.
     &        (filab(11)(1:4).eq.'PU'))  then 
            xr=v(1,j)*a(1,1)+v(2,j)*a(1,2)+v(3,j)*a(1,3)
            xt=v(1,j)*a(2,1)+v(2,j)*a(2,2)+v(3,j)*a(2,3)
            xz=v(1,j)*a(3,1)+v(2,j)*a(3,2)+v(3,j)*a(3,3)
            v(1,j)=xr
            v(2,j)=xt
            v(3,j)=xz
         endif
!     
         if((filab(3)(1:4).eq.'S   ').or.
     &        (filab(18)(1:4).eq.'PHS ')) then
            b(1,1)=stn(1,j)*a(1,1)+stn(4,j)*a(1,2)+stn(5,j)*a(1,3)
            b(1,2)=stn(1,j)*a(2,1)+stn(4,j)*a(2,2)+stn(5,j)*a(2,3)
            b(1,3)=stn(1,j)*a(3,1)+stn(4,j)*a(3,2)+stn(5,j)*a(3,3)
            b(2,1)=stn(4,j)*a(1,1)+stn(2,j)*a(1,2)+stn(6,j)*a(1,3)
            b(2,2)=stn(4,j)*a(2,1)+stn(2,j)*a(2,2)+stn(6,j)*a(2,3)
            b(2,3)=stn(4,j)*a(3,1)+stn(2,j)*a(3,2)+stn(6,j)*a(3,3)
            b(3,1)=stn(5,j)*a(1,1)+stn(6,j)*a(1,2)+stn(3,j)*a(1,3)
            b(3,2)=stn(5,j)*a(2,1)+stn(6,j)*a(2,2)+stn(3,j)*a(2,3)
            b(3,3)=stn(5,j)*a(3,1)+stn(6,j)*a(3,2)+stn(3,j)*a(3,3)
!     
            stn(1,j)=a(1,1)*b(1,1)+a(1,2)*b(2,1)+a(1,3)*b(3,1)
            stn(2,j)=a(2,1)*b(1,2)+a(2,2)*b(2,2)+a(2,3)*b(3,2)
            stn(3,j)=a(3,1)*b(1,3)+a(3,2)*b(2,3)+a(3,3)*b(3,3)
            stn(4,j)=a(1,1)*b(1,2)+a(1,2)*b(2,2)+a(1,3)*b(3,2)
            stn(5,j)=a(1,1)*b(1,3)+a(1,2)*b(2,3)+a(1,3)*b(3,3)
            stn(6,j)=a(2,1)*b(1,3)+a(2,2)*b(2,3)+a(2,3)*b(3,3)
         endif
!
         if(filab(4)(1:4).eq.'E   ') then
            b(1,1)=een(1,j)*a(1,1)+een(4,j)*a(1,2)+een(5,j)*a(1,3)
            b(1,2)=een(1,j)*a(2,1)+een(4,j)*a(2,2)+een(5,j)*a(2,3)
            b(1,3)=een(1,j)*a(3,1)+een(4,j)*a(3,2)+een(5,j)*a(3,3)
            b(2,1)=een(4,j)*a(1,1)+een(2,j)*a(1,2)+een(6,j)*a(1,3)
            b(2,2)=een(4,j)*a(2,1)+een(2,j)*a(2,2)+een(6,j)*a(2,3)
            b(2,3)=een(4,j)*a(3,1)+een(2,j)*a(3,2)+een(6,j)*a(3,3)
            b(3,1)=een(5,j)*a(1,1)+een(6,j)*a(1,2)+een(3,j)*a(1,3)
            b(3,2)=een(5,j)*a(2,1)+een(6,j)*a(2,2)+een(3,j)*a(2,3)
            b(3,3)=een(5,j)*a(3,1)+een(6,j)*a(3,2)+een(3,j)*a(3,3)
!     
            een(1,j)=a(1,1)*b(1,1)+a(1,2)*b(2,1)+a(1,3)*b(3,1)
            een(2,j)=a(2,1)*b(1,2)+a(2,2)*b(2,2)+a(2,3)*b(3,2)
            een(3,j)=a(3,1)*b(1,3)+a(3,2)*b(2,3)+a(3,3)*b(3,3)
            een(4,j)=a(1,1)*b(1,2)+a(1,2)*b(2,2)+a(1,3)*b(3,2)
            een(5,j)=a(1,1)*b(1,3)+a(1,2)*b(2,3)+a(1,3)*b(3,3)
            een(6,j)=a(2,1)*b(1,3)+a(2,2)*b(2,3)+a(2,3)*b(3,3)
         endif
!     
         if(filab(5)(1:4).eq.'RF  ') then
            xr=fn(1,j)*a(1,1)+fn(2,j)*a(1,2)+fn(3,j)*a(1,3)
            xt=fn(1,j)*a(2,1)+fn(2,j)*a(2,2)+fn(3,j)*a(2,3)
            xz=fn(1,j)*a(3,1)+fn(2,j)*a(3,2)+fn(3,j)*a(3,3)
            fn(1,j)=xr
            fn(2,j)=xt
            fn(3,j)=xz
         endif
!     
         if(filab(9)(1:4).eq.'HFL ') then
            xr=qfn(1,j)*a(1,1)+qfn(2,j)*a(1,2)+qfn(3,j)*a(1,3)
            xt=qfn(1,j)*a(2,1)+qfn(2,j)*a(2,2)+qfn(3,j)*a(2,3)
            xz=qfn(1,j)*a(3,1)+qfn(2,j)*a(3,2)+qfn(3,j)*a(3,3)
            qfn(1,j)=xr
            qfn(2,j)=xt
            qfn(3,j)=xz
         endif
!     
         if(filab(32)(1:4).eq.'ME  ') then
            b(1,1)=emn(1,j)*a(1,1)+emn(4,j)*a(1,2)+emn(5,j)*a(1,3)
            b(1,2)=emn(1,j)*a(2,1)+emn(4,j)*a(2,2)+emn(5,j)*a(2,3)
            b(1,3)=emn(1,j)*a(3,1)+emn(4,j)*a(3,2)+emn(5,j)*a(3,3)
            b(2,1)=emn(4,j)*a(1,1)+emn(2,j)*a(1,2)+emn(6,j)*a(1,3)
            b(2,2)=emn(4,j)*a(2,1)+emn(2,j)*a(2,2)+emn(6,j)*a(2,3)
            b(2,3)=emn(4,j)*a(3,1)+emn(2,j)*a(3,2)+emn(6,j)*a(3,3)
            b(3,1)=emn(5,j)*a(1,1)+emn(6,j)*a(1,2)+emn(3,j)*a(1,3)
            b(3,2)=emn(5,j)*a(2,1)+emn(6,j)*a(2,2)+emn(3,j)*a(2,3)
            b(3,3)=emn(5,j)*a(3,1)+emn(6,j)*a(3,2)+emn(3,j)*a(3,3)
!     
            emn(1,j)=a(1,1)*b(1,1)+a(1,2)*b(2,1)+a(1,3)*b(3,1)
            emn(2,j)=a(2,1)*b(1,2)+a(2,2)*b(2,2)+a(2,3)*b(3,2)
            emn(3,j)=a(3,1)*b(1,3)+a(3,2)*b(2,3)+a(3,3)*b(3,3)
            emn(4,j)=a(1,1)*b(1,2)+a(1,2)*b(2,2)+a(1,3)*b(3,2)
            emn(5,j)=a(1,1)*b(1,3)+a(1,2)*b(2,3)+a(1,3)*b(3,3)
            emn(6,j)=a(2,1)*b(1,3)+a(2,2)*b(2,3)+a(2,3)*b(3,3)
         endif
!     
      enddo
!     
      return
      end
