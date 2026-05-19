!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine extrapol_vel2(ipnei,neifa,vfa,area,xxna,volume,gradvel,
     &  nefa,nefb,ncfd)
!
!     calculate the gradient of the velocities at the center of
!     the elements
!
      implicit none
!
      integer ipnei(*),neifa(*),nefa,nefb,i,k,l,ifa,indexf,ncfd
!
      real*8 vfa(0:7,*),area(*),xxna(3,*),volume(*),gradvel(3,3,*)
!
!
!
      do i=nefa,nefb
!
!           initialization
!
         do k=1,ncfd
            do l=1,ncfd
               gradvel(k,l,i)=0.d0
            enddo
         enddo
!
         do indexf=ipnei(i)+1,ipnei(i+1)
            ifa=neifa(indexf)
            do k=1,ncfd
               do l=1,ncfd
                  gradvel(k,l,i)=gradvel(k,l,i)+
     &                 vfa(k,ifa)*xxna(l,indexf)
               enddo
            enddo
         enddo
!     
!     dividing by the volume of the element
!     
         do k=1,ncfd
            do l=1,ncfd
               gradvel(k,l,i)=gradvel(k,l,i)/volume(i)
            enddo
         enddo
      enddo
c      do i=1,120
c         write(*,*) 'extrapol_vel2',i,gradvel(1,1,i),gradvel(1,2,i)
c      enddo
!
      return
      end
