!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine writeview(ntr,adview,auview,fenv,nzsrad,
     &  jobnamef)
!     
!     writing the viewfactors to file
!
      implicit none
!
      character*80 version
      character*132 jobnamef(*),fnvw
!     
      integer ntr,nzsrad,i,k
!
      real*8 adview(*),auview(*),fenv(*)
!     
      write(*,*) 'Writing the viewfactors to file'
      write(*,*)
!     
      if(jobnamef(3)(1:1).eq.' ') then
         do i=1,132
            if(jobnamef(1)(i:i).eq.' ') exit
         enddo
         i=i-1
         fnvw=jobnamef(1)(1:i)//'.vwf'
      else
         fnvw=jobnamef(3)
      endif
      open(10,file=fnvw,status='unknown',form='unformatted',
     &     access='sequential',err=10)
!
      call getversion(version)
!     
      write(10) version
      write(10) (adview(k),k=1,ntr)
      write(10) (auview(k),k=1,2*nzsrad)
      write(10)(fenv(k),k=1,ntr)
      close(10)
!     
      return
!
 10   write(*,*) '*ERROR in writeview: could not open file ',fnvw
      call exit(201)
      end
      
