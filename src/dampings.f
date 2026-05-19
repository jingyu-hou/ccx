!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine dampings(inpc,textpart,xmodal,istep,
     &  istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc,irstrt,ier,
     &  dacon,nmat)
!
!     reading the input deck: *DAMPING
!
      implicit none
!
      character*1 inpc(*)
      character*132 textpart(16)
!
      integer istep,istat,n,key,iline,ipol,inl,ipoinp(2,*),
     &  inp(3,*),ipoinpc(0:*),i,irstrt(*),ier,nmat
!
      real*8 xmodal(*),dacon(*)
!
      if((istep.gt.0).and.(irstrt(1).ge.0)) then
         write(*,*) '*ERROR reading *DAMPING: *DAMPING should be placed'
         write(*,*) '       before all step definitions'
         ier=1
         return
      endif
!
      if(nmat.eq.0) then
         write(*,*) 
     &       '*ERROR reading *DAMPING: *DAMPING should be preceded'
         write(*,*) '  by a *MATERIAL card'
         ier=1
         return
      endif
!
      do i=2,n
         if(textpart(i)(1:6).eq.'ALPHA=') then
            read(textpart(i)(7:26),'(f20.0)',iostat=istat) xmodal(1)
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*DAMPING%",ier)
               return
            endif
         elseif(textpart(i)(1:5).eq.'BETA=') then
            read(textpart(i)(6:25),'(f20.0)',iostat=istat) xmodal(2)
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*DAMPING%",ier)
               return
            endif
         elseif(textpart(i)(1:11).eq.'STRUCTURAL=') then
            read(textpart(i)(12:31),'(f20.0)',iostat=istat) dacon(nmat)
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*DAMPING%",ier)
               return
            endif
         else
            write(*,*) 
     &        '*WARNING reading *DAMPING: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*DAMPING%")
         endif
      enddo
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end

