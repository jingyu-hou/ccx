!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine noanalysiss(inpc,textpart,nmethod,iperturb,istep,
     &  istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc,tper,ier)
!
!     reading the input deck: *NO ANALYSIS
!
      implicit none
!
      character*1 inpc(*)
      character*132 textpart(16)
!
      integer nmethod,iperturb,istep,istat,n,key,iline,ipol,inl,
     &  ipoinp(2,*),inp(3,*),ipoinpc(0:*),ier
!
      real*8 tper
!
      if(istep.lt.1) then
         write(*,*)
     &      '*ERROR reading *NO ANALYSIS: *NO ANALYSIS can only be used'
         write(*,*) '  within a STEP'
         ier=1
         return
      endif
!
      write(*,*) '*WARNING: no analysis option was chosen'
!
      nmethod=0
      iperturb=0
      tper=1.d0
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end

