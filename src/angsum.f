!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!  
      subroutine angsum(lakon,kon,ipkon,neigh,ipneigh,co,node,itypflag,
     &     angle)
!
!     computes the sum of all spaceangles of the element edges
!     adjacent to a node
!
!     author: Sascha Merz
!
      implicit none
!
      integer kon(*),ipkon(*),ielem,j,k,lneigh4(3,4),indexe,nnr,
     & neigh(2,*),ipneigh(*),index,m,lneigh8(3,8),nvertex,node,itypflag
!
      real*8 co(3,*),angle,cotet(3,4),spaceangle
!     
      data lneigh4 /2,3,4,1,3,4,1,2,4,1,2,3/
!
      data lneigh8 /2,4,5,1,3,6,2,4,7,1,3,8,
     &             1,6,8,2,5,7,3,6,8,4,5,7/
!
      character*8 lakon(*)
!
      index=ipneigh(node)
!
      angle=0.d0
      do j=1,3
         cotet(j,1)=co(j,node)
      enddo
      do
         if(index.eq.0) exit
         ielem=neigh(1,index)
!
         if(lakon(ielem)(1:5).eq.'C3D20'.and.itypflag.eq.1) then
            nvertex=8
         elseif(lakon(ielem)(1:5).eq.'C3D10'.and.itypflag.eq.2) then
            nvertex=4
         elseif(lakon(ielem)(1:4).eq.'C3D8'.and.itypflag.eq.3) then
            nvertex=8
         else
            index=neigh(2,index)
            cycle
         endif
!
         indexe=ipkon(ielem)
         do m=1,nvertex
            if(kon(indexe+m).eq.node) exit
         enddo
         do j=1,3
            if(nvertex.eq.4) then
               nnr=kon(indexe+lneigh4(j,m))
            elseif(nvertex.eq.8) then
               nnr=kon(indexe+lneigh8(j,m))            
            endif
            do k=1,3
               cotet(k,j+1)=co(k,nnr)
            enddo
         enddo
         angle=angle+spaceangle(cotet)
         index=neigh(2,index)
      enddo
!
      return
      end
!      
      real*8 function spaceangle(cotet)
!
      implicit none
!
      integer i,j
      real*8 vector(3,3),ca,cb,cc,ca2,cb2,cc2,sa,sb,sc,cosa,sina,
     &     cotanb,cotanc,a,b,c,cotet(3,4),absval
!     calculate normal vectors
      do i=1,3
!     i is vector 1, 2 and 3
         do j=1,3
!     j is x,y,z
            vector(j,i)=cotet(j,i+1)-cotet(j,1)
         enddo
        absval=dsqrt(vector(1,i)*vector(1,i)
     &       +vector(2,i)*vector(2,i)
     &       +vector(3,i)*vector(3,i))
         do j=1,3
            vector(j,i)=vector(j,i)/absval
!            write(*,*) 'vektor ij',i,j,' ist',vector(i,j)
         enddo
      enddo
!
      ca=vector(1,1)*vector(1,2)+vector(2,1)*vector(2,2)+
     &     vector(3,1)*vector(3,2)
      cb=vector(1,2)*vector(1,3)+vector(2,2)*vector(2,3)+
     &     vector(3,2)*vector(3,3)
      cc=vector(1,1)*vector(1,3)+vector(2,1)*vector(2,3)+
     &     vector(3,1)*vector(3,3)
!
      ca2=min(ca*ca,1.d0)
      cb2=min(cb*cb,1.d0)
      cc2=min(cc*cc,1.d0)
!
      sa=dsqrt(1.d0-ca2)
      sb=dsqrt(1.d0-cb2)
      sc=dsqrt(1.d0-cc2)
!
      if((dabs(sa).lt.1.d-8).or.(dabs(sb).lt.1.d-8).or.
     &   (dabs(sc).lt.1.d-8)) then
        spaceangle=0.d0
        return
      endif
!
!      sa=dsqrt(1.d0-ca*ca)
!      sb=dsqrt(1.d0-cb*cb)
!      sc=dsqrt(1.d0-cc*cc)
!
      cosa=(ca-cb*cc)/(sb*sc)
      sina=dsqrt(1.d0-cosa*cosa)
      cotanb=(sc*cb/sb-cosa*cc)/sina
      cotanc=(sb*cc/sc-cosa*cb)/sina
!
      a=dacos(cosa)
      b=datan(1.d0/cotanb)
      c=datan(1.d0/cotanc)
!
      if(b.lt.0) b=b+3.141592653589793d0
      if(c.lt.0) c=c+3.141592653589793d0
!
      spaceangle=a+b+c-3.141592653589793d0
!
      return
      end
