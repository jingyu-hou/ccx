!
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine e_c3d_vlhs(co,nk,konl,lakonl,sm,nelem,ipvar,var,
     &     smscalel)
!
!     computation of the velocity element matrix for the element with
!     the topology in konl
!
      implicit none
!
      character*8 lakonl
!
      integer konl(20),nk,nelem,i,j,nmethod,ii,jj,ii1,jj1,kk,
     &  nope,mint3d,iflag,ipvar(*),index
!
      real*8 co(3,*),xl(3,20),shp(4,20),xi,et,ze,xsj,smscalel,
     &   sm(78,78),weight,var(*)
!
!
!      
      include "gauss.f"
!
      data iflag /2/
!
      if(lakonl(4:4).eq.'8') then
         nope=8
      elseif(lakonl(4:4).eq.'4') then
         nope=4
      elseif(lakonl(4:4).eq.'6') then
         nope=6
      endif
!
      if(lakonl(4:5).eq.'8R') then
         mint3d=1
      elseif(lakonl(4:4).eq.'8') then
        mint3d=8
      elseif(lakonl(4:4).eq.'4') then
         mint3d=1
      elseif(lakonl(4:5).eq.'6 ') then
         mint3d=2
      endif
!
!     computation of the coordinates of the local nodes
!
      do i=1,nope
        do j=1,3
          xl(j,i)=co(j,konl(i))
        enddo
      enddo
!
!     initialisation of sm
!
      do i=1,3*nope
         do j=1,3*nope
            sm(i,j)=0.d0
         enddo
      enddo
!
!     computation of the matrix: loop over the Gauss points
!
      index=ipvar(nelem)
      do kk=1,mint3d
         if(lakonl(4:5).eq.'8R') then
            weight=weight3d1(kk)
         elseif(lakonl(4:4).eq.'8') then
            weight=weight3d2(kk)
         elseif(lakonl(4:4).eq.'4') then
            weight=weight3d4(kk)
         elseif(lakonl(4:5).eq.'6 ') then
            weight=weight3d7(kk)
         endif
!
!        copying the shape functions, their derivatives and the
!        Jacobian determinant from field var
!
         do jj=1,nope
            do ii=1,4
               index=index+1
               shp(ii,jj)=var(index)
            enddo
         enddo
         index=index+1
         xsj=var(index)
         
         index=index+nope+14
!
         weight=weight*xsj
c         weight=weight*xsj*smscalel
!     
         jj1=1
         do jj=1,nope
!     
            ii1=1
            do ii=1,jj
!     
!              lhs velocity matrix
!     
               sm(ii1,jj1)=sm(ii1,jj1)
     &              +shp(4,ii)*shp(4,jj)*weight
               sm(ii1+1,jj1+1)=sm(ii1,jj1)
               sm(ii1+2,jj1+2)=sm(ii1,jj1)
!     
               ii1=ii1+3
            enddo
            jj1=jj1+3
         enddo
      enddo
!
      return
      end

