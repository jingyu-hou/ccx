!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine valuesatinfinitys(inpc,textpart,physcon,
     &  istep,istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc,ier)
!
!     reading the input deck: *VALUES AT INFINITY
!
!     physcon(4): static temperature at infinity
!     physcon(5): norm of the velocity at infinity
!     physcon(6): static pressure at infinity
!     physcon(7): density at infinity
!     physcon(8): length of the computational domain
!
      implicit none
!
      character*1 inpc(*)
      character*132 textpart(16)
!
      integer i,istep,istat,n,key,iline,ipol,inl,ipoinp(2,*),inp(3,*),
     &  ipoinpc(0:*),ier
!
      real*8 physcon(*)
!
      if(istep.gt.0) then
         write(*,*) 
     &   '*ERROR reading *VALUES AT INFINITY: *VALUES AT INFINITY'
         write(*,*) '        should only be used before the first STEP'
         ier=1
         return
      endif
!
      do i=2,n
         write(*,*) 
     &  'WARNING reading *VALUES AT INFINITY: parameter not recognized:'
         write(*,*) '         ',
     &        textpart(i)(1:index(textpart(i),' ')-1)
         call inputwarning(inpc,ipoinpc,iline,
     &"*VALUES AT INFINITY%")
      enddo
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      do i=1,5
         read(textpart(i),'(f20.0)',iostat=istat) physcon(3+i)
         if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*VALUES AT INFINITY%",ier)
            return
         endif
      enddo
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end







