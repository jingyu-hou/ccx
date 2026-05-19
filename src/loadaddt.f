!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine loadaddt(nelement,label,valfilm,valtemp,nelemload,
     &  sideload,xload,nload,nload_,iamload,iamptemp,
     &  iampfilm,nam,node,iload)
!
!     adds a thermal dload condition to the data base
!
      implicit none
!
      character*20 label,sideload(*)
!
      integer nelemload(2,*),iamload(2,*),id,iload,
     &  nelement,nload,nload_,j,iamptemp,nam,iampfilm,node
!
      real*8 xload(2,*),valfilm,valtemp
!
      call nident2(nelemload,nelement,nload,id)
      if(id.gt.0) then
!
!        it is possible that several *DLOAD, *FILM or
!        *RADIATE boundary conditions are applied to one
!        and the same element
!
         if(nelemload(1,id).eq.nelement) then
            do
               if(sideload(id).eq.label) then
!     
!     loading on same element face detected: values
!     are replaced
!     
                  xload(1,id)=valfilm
                  xload(2,id)=valtemp
                  nelemload(2,id)=node
                  if(nam.gt.0) then
                     iamload(1,id)=iampfilm
                     iamload(2,id)=iamptemp
                  endif
                  iload=id
                  return
cc               elseif(sideload(id).lt.label) then
cc                  exit
               endif
               id=id-1
               if((id.eq.0).or.(nelemload(1,id).ne.nelement)) then
                  exit
               endif
            enddo
         endif
      endif
!
!     loading a element face on which no previous loading
!     was applied
!
!     loading conditions on one and the same element are
!     alphabetized based on field sideload
!
      nload=nload+1
      if(nload.gt.nload_) then
         write(*,*) '*ERROR in loadadd: increase nload_'
         call exit(201)
      endif
!
!     shifting existing loading
!
      do j=nload,id+2,-1
         nelemload(1,j)=nelemload(1,j-1)
         nelemload(2,j)=nelemload(2,j-1)
         sideload(j)=sideload(j-1)
         xload(1,j)=xload(1,j-1)
         xload(2,j)=xload(2,j-1)
         if(nam.gt.0) then
            iamload(1,j)=iamload(1,j-1)
            iamload(2,j)=iamload(2,j-1)
         endif
      enddo
!
!     inserting new loading
!
      nelemload(1,id+1)=nelement
      nelemload(2,id+1)=node
      sideload(id+1)=label
      xload(1,id+1)=valfilm
      xload(2,id+1)=valtemp
      if(nam.gt.0) then
         iamload(1,id+1)=iampfilm
         iamload(2,id+1)=iamptemp
      endif
      iload=id+1
!
      return
      end

