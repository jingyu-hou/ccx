!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine liquidfractions(inpc,textpart,flcon,nflcon,lh,
     &  nmat,ntmat_,irstrt,istep,istat,n,iline,ipol,inl,ipoinp,inp,
     &  ipoinpc,ier)
!
!     reading the input deck: *LIQUIDFRACTION
!
      implicit none
!
      character*1 inpc(*)
      character*132 textpart(16)
!
      integer nflcon(*),nmat,ntmat,ntmat_,istep,istat,n,ipoinpc(0:*),
     &  key,irstrt(*),iline,ipol,inl,ipoinp(2,*),inp(3,*),i,ier
!
      real*8 flcon(0:1,ntmat_,*),lh(*)
!
      ntmat=0
      lh(nmat)=0.d0
!
      if((istep.gt.0).and.(irstrt(1).ge.0)) then
         write(*,*) 
     &      '*ERROR reading *fl: *fl should be placed'
         write(*,*) '  before all step definitions'
         ier=1
         return
      endif
!
      if(nmat.eq.0) then
         write(*,*) 
     &    '*ERROR reading *fl: *fl should be preceded'
         write(*,*) '  by a *MATERIAL card'
         ier=1
         return
      endif
!
      do i=2,n
         if(textpart(i)(1:11).eq.'LATANTHEAT=') then
            read(textpart(i)(12:25),'(f20.0)',iostat=istat) lh(nmat)
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*LIQUIDFRACTION%",ier)
               return
            endif
         else
            write(*,*) 
     &        '*WARNING reading *fl: parameter not recognized:'
            write(*,*) '         ',
     &        textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*fl%")
         endif
      enddo
!
      do
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) return
         ntmat=ntmat+1
         nflcon(nmat)=ntmat
         if(ntmat.gt.ntmat_) then
            write(*,*) '*ERROR reading *fl: increase ntmat_'
            ier=1
            return
         endif
         read(textpart(1)(1:20),'(f20.0)',iostat=istat) 
     &            flcon(1,ntmat,nmat)
         if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*fl%",ier)
            return
         endif
         read(textpart(2)(1:20),'(f20.0)',iostat=istat) 
     &            flcon(0,ntmat,nmat)
         if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*fl%",ier)
            return
         endif
      enddo
!
      return
      end

