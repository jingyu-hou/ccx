!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine identifytransform(nelement,label,nelemload,sideload,
     &  nload,loadid)
!
!     checks whether a transformation was applied to a given face,
!     and if so, which one (only for CFD)
!
      implicit none
!
      character*20 label,sideload(*)
!
      integer nelemload(2,*),nelement,nload,id,loadid
!
      loadid=0
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
               if (sideload(id).eq.label) then
                  loadid=id
               elseif(sideload(id).lt.label) then
                  exit
               endif
               id=id-1
               if((id.eq.0).or.(nelemload(1,id).ne.nelement)) then
                  exit
               endif
            enddo
         endif
      endif
!
      return
      end

