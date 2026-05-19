!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine writeelem(i,lakon)
!
!     this routine is called if an inconsistency is noticed between
!     the element count and the number of elements stored in the frd-
!     file. Such an inconsistency will lead to a crash while reading
!     a binary frd-file
!
      implicit none
!
      character*8 lakon(*)
!
      integer i
!
      write(*,*) '*ERROR in writeelem:'
      write(*,*) '       element ',i+1,' with label ',lakon(i+1)
      write(*,*) '       is not stored in the frd-file. Yet, '
      write(*,*) '       it is taken into account in the element'
      write(*,*) '       count: inconsistency. Please contact the'
      write(*,*) '       author of WeICME'
      call exit(201)
!
      return
      end

