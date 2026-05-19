!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine writesen(g0,dgdx,ndesi,nobject,nodedesi,jobnamef)
!
!     writes "raw" sensitivities to file jobname.sen
!
      implicit none
!
      character*132 jobnamef(*),cfile
!
      integer ndesi,nobject,nodedesi(*),i,j
!
      real*8 g0(*),dgdx(ndesi,*)
!
      intent(in) g0,dgdx,ndesi,nobject,nodedesi,jobnamef
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
      write(27,*) (g0(j),j=1,nobject)
!
      close(27)
!
!     storing the sensitivity of the objectives
!
      cfile(i+4:i+4)='1'
      open(27,file=cfile,status='unknown')
!
      do i=1,ndesi
         write(27,*) nodedesi(i),(dgdx(i,j),j=1,nobject)
      enddo
!
      close(27)
!
      return
      end

