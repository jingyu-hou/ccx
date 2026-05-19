!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine modaldampings(inpc,textpart,nmethod,xmodal,istep,
     &  istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc,ier)
!
!     reading the input deck: *MODAL DAMPING
!
      implicit none
!
      logical rayleigh
!
      character*1 inpc(*)
      character*132 textpart(16)
!
      integer nmethod,istep,istat,n,key,iline,ipol,inl,ipoinp(2,*),
     &  inp(3,*),ipoinpc(0:*),i,lowfrequ,highfrequ,k,ier
!
      real*8 xmodal(*),zeta
!
      if(istep.lt.1) then
         write(*,*) 
     &      '*ERROR reading *MODAL DAMPING: *MODAL DAMPING can only'
         write(*,*) '  be used within a STEP'
         ier=1
         return
      endif
!
      rayleigh=.false.
      do i=2,n
         if(textpart(i)(1:8).eq.'RAYLEIGH') then
            rayleigh=.true.
         else
            write(*,*) 
     &      '*WARNING reading *MODAL DAMPING: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*MODAL DAMPING%")
         endif
      enddo
      if(rayleigh) then
         xmodal(11)=-0.5d0
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) then
            write(*,*) '*ERROR reading *MODAL DAMPING: definition 
     &                  not complete'
            write(*,*) '       '
            call inputerror(inpc,ipoinpc,iline,
     &           "*MODAL DAMPING%",ier)
            return
         endif
         read(textpart(3)(1:20),'(f20.0)',iostat=istat) xmodal(1)
         if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*MODAL DAMPING%",ier)
            return
         endif
         read(textpart(4)(1:20),'(f20.0)',iostat=istat) xmodal(2)
         if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*MODAL DAMPING%",ier)
            return
         endif
!
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
!
      else
         do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) exit
!
            read(textpart(1)(1:10),'(i10)',iostat=istat) lowfrequ
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*MODAL DAMPING%",ier)
               return
            endif
            read(textpart(2)(1:10),'(i10)',iostat=istat) highfrequ
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*MODAL DAMPING%",ier)
               return
            endif
            read(textpart(3)(1:20),'(f20.0)',iostat=istat) zeta
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*MODAL DAMPING%",ier)
               return
            endif
!
            if(highfrequ<lowfrequ) highfrequ=lowfrequ  
            do k=lowfrequ,highfrequ
               xmodal(11+k)=zeta
            enddo
         enddo  
      endif
      return
      end

