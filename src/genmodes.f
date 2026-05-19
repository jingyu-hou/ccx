!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
c     Bernhardi start
      subroutine genmodes(i,kon,ipkon,lakon,ne,nk,nk_,co)
!
!     generate nodes for incompatible modes
!
      implicit none
!
      character*8 lakon(*)
!
      real*8 co(3,*),coords(3)
!
      integer i,kon(*),ipkon(*),ne,nope,nopeexp,
     &  nk,nk_,j,indexe,k,nodeb(8,3)
!
      indexe=ipkon(i)
!
      if(lakon(i)(1:5).eq.'C3D8I')then
         nope=8
         nopeexp=3
      else
         write(*,*) "*ERROR in genmodes: wrong element type, element=",
     &               lakon(i)
         call exit(201)
      endif
!
!     generating additional nodes for the incompatible element. 
!      
!     determining the mean value of the coordinates of the element
!      
      do k=1,3
         coords(k)=0.d0
         do j=1,nope
            coords(k)=coords(k)+co(k,kon(indexe+j))
         enddo
         coords(k)=coords(k)/8.d0
      enddo
!
      do j=1,nopeexp
         nk=nk+1
           if(nk.gt.nk_) then
              write(*,*) '*ERROR in genmodes: increase nk_'
              call exit(201)
           endif
         kon(indexe+nope+j)=nk
         do k=1,3
            co(k,nk)=coords(k)
         enddo
      enddo
!
      return
      end
c     Bernhardi end

