!
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine test_cast(kode,time,nk,mi,inomat,vold,coolingrate,
     &  gtemp,ntmat_,ithermal,istep,iinc,icounter)
!
      implicit none
!
      character*3 m1,m3
      character*8 fmat
      character*132 text
!
      integer kode,nk,mi(*),inomat(*),ntmat_,ithermal(*),
     &        nout,i,istep,iinc,icounter
!
      real*8 time,vold(0:mi(2),*),coolingrate(*),gtemp(*)
!
      save nout
!
      m1=' -1'
      m3=' -3'
!
      nout=0
      do i=1,nk
         if(inomat(i).le.0) cycle
         nout=nout+1
      enddo
!
      if(time.le.0.d0) then
         fmat(1:8)='(e12.5) '
      elseif((dlog10(time).ge.0.d0).and.(dlog10(time).lt.11.d0)) then
         fmat(1:5)='(f12.'
         write(fmat(6:7),'(i2)') 11-int(dlog10(time)+1.d0)
         fmat(8:8)=')'
      else
         fmat(1:8)='(e12.5) '
      endif
!
      text='    1PSTEP'
c      write(text(25:36),'(i12)') kode
c      write(13,'(a132)') text
      write(text(25:36),'(3i12)') icounter
      write(text(37:48),'(i12)') iinc
      write(text(49:60),'(i12)') istep
      icounter=icounter+1
      write(13,'(a80)') text
!
!     write Tdot data 
!
      text=
     & '  100CL       .00000E+00                                 3    1'
      text(75:75)='1'
      write(text(25:36),'(i12)') nout
      write(text(8:12),'(i5)') 100+kode
      write(text(13:24),fmat) time
      write(text(59:63),'(i5)') kode
      write(13,'(a132)') text
      text=' -4  TDOT        1    1'
      write(13,'(a132)') text
      text=' -5  TS          1    1    0    0'
      write(13,'(a132)') text
!
      do i=1,nk
         if(inomat(i).le.0) cycle       
         write(13,100) m1,i,coolingrate(i)
      enddo
!
      write(13,'(a3)') m3
!
!     write G data 
!
      text=
     & '  100CL       .00000E+00                                 3    1'
      text(75:75)='1'
      write(text(25:36),'(i12)') nout
      write(text(8:12),'(i5)') 100+kode
      write(text(13:24),fmat) time
      write(text(59:63),'(i5)') kode
      write(13,'(a132)') text
      text=' -4  GTEMP       1    1'
      write(13,'(a132)') text
      text=' -5  TS          1    1    0    0'
      write(13,'(a132)') text
!
      do i=1,nk
         if(inomat(i).le.0) cycle     
         write(13,100) m1,i,gtemp(i)
      enddo
!
      write(13,'(a3)') m3
!
 100  format(a3,i10,1p,6e12.5)
!
      return
      end
!