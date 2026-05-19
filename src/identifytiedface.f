!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine identifytiedface(tieset,ntie,set,nset,ifaceslave,kind)
!
!     identifies slave nodes in tied slave faces
!
      implicit none
!
      character*1 kind
      character*81 tieset(3,*),slavset,set(*)
!
      integer ifaceslave(*),i,j,nset,ipos,ntie
!
!     nodes per face for tet elements
!
      do i=1,ntie
         if(tieset(1,i)(81:81).ne.kind) cycle
         slavset=tieset(2,i)
         ipos=index(slavset,' ')
         slavset(ipos:ipos)='T'
         do j=1,nset
            if(set(j).eq.slavset) exit
         enddo
         if(j.gt.nset) then
            slavset(ipos:ipos)='S'
            do j=1,nset
               if(set(j).eq.slavset) then
                  exit
               endif
            enddo
            if(j.gt.nset) then
               write(*,*) 
     &           '*ERROR in identifytiedface: ',
     &           'tied contact nodal slave surface ',
     &              slavset
               write(*,*) '       does not exist'
               call exit(201)
            else
               tieset(2,i)(ipos:ipos)='S'
               ifaceslave(i)=0
            endif
         else
            tieset(2,i)(ipos:ipos)='T'
            ifaceslave(i)=1
         endif
      enddo
      return
      end

