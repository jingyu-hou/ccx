!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine readsen(g0,dgdx,ndesi,nobject,nodedesi,jobnamef)
!
!     reads "raw" sensitivities to file jobname.sen
!
      implicit none
!
      character*132 jobnamef(*),cfile
!
      integer ndesi,nobject,nodedesi(*),i,j,idummy
!
      real*8 g0(*),dgdx(ndesi,*)
!
      intent(in) ndesi,nobject,nodedesi,jobnamef
!
      intent(inout) g0,dgdx
!
!     storing the objectives
!
      do i=1,132
         cfile(i:i)=' '
      enddo
      do i=1,132
         if(jobnamef(1)(i:i).eq.' ') exit
         cfile(i:i)=jobnamef(1)(i:i)
      enddo
      cfile(i:i+4)='.sen0'
      open(27,file=cfile,status='unknown')
!
c      read(27,*) g0(1)
      read(27,*) (g0(j),j=1,nobject)
!
      close(27)
!
!     storing the sensitivity of the objectives
!
      cfile(i+4:i+4)='1'
      open(27,file=cfile,status='unknown')
!
      do i=1,ndesi
         read(27,*) idummy,(dgdx(i,j),j=1,nobject)
         if(idummy.ne.nodedesi(i)) then
            write(*,*) '*ERROR in readsen: design nodes not'
            write(*,*) '       in correct ascending order in'
            write(*,*) '       file',cfile
         endif
      enddo
!
      close(27)
!
      return
      end

