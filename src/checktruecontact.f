!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine checktruecontact(ntie,tieset,tietol,elcon,itruecontact,
     &  ncmat_,ntmat_)
!
!     check whether for face-to-face penalty contact the
!     surface behavior definition is such that true contact is 
!     defined and not just tied contact
!
      implicit none
!
      character*81 tieset(3,*)
!
      integer itruecontact,i,ntie,imat,ncmat_,ntmat_
!
      real*8 tietol(3,*),elcon(0:ncmat_,ntmat_,*)
!
!     if at least one tied contact is present, itruecontact
!     is set to zero and no check is performed whether tension
!     occurs in the contact areas
!
      itruecontact=1
      do i=1,ntie
         if(tieset(1,i)(81:81).eq.'C') then
            imat=int(tietol(2,i))
            if(int(elcon(3,1,imat)).eq.4) then
               itruecontact=0
               exit
            endif
         endif
      enddo
!
      return
      end
