!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine specificgasconstants(inpc,textpart,shcon,nshcon,
     &  nmat,ntmat_,irstrt,istep,istat,n,iline,ipol,inl,ipoinp,
     &  inp,ipoinpc,ier)
!
!     reading the input deck: *SPECIFIC GAS CONSTANT
!
      implicit none
!
      character*1 inpc(*)
      character*132 textpart(16)
!
      integer nshcon(*),nmat,ntmat_,istep,istat,n,ipoinpc(0:*),
     &  key,irstrt(*),iline,ipol,inl,ipoinp(2,*),inp(3,*),i,ier
!
      real*8 shcon(0:3,ntmat_,*)
!
      if((istep.gt.0).and.(irstrt(1).ge.0)) then
         write(*,*) 
     &   '*ERROR reading *SPECIFIC GAS CONSTANT: *SPECIFIC GAS CONSTANT'
         write(*,*) '  should be placed before all step definitions'
         ier=1
         return
      endif
!
      if(nmat.eq.0) then
         write(*,*) 
     &   '*ERROR reading *SPECIFIC GAS CONSTANT: *SPECIFIC GAS CONSTANT'
         write(*,*) '  should be preceded by a *MATERIAL card'
         ier=1
         return
      endif
!
      do i=2,n
         write(*,*) 
     &     '*WARNING in specificgasconstants: parameter not recognized:'
         write(*,*) '         ',
     &        textpart(i)(1:index(textpart(i),' ')-1)
         call inputwarning(inpc,ipoinpc,iline,
     &"*SPECIFIC GAS CONSTANT%")
      enddo
!
      do
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) return
         read(textpart(1)(1:20),'(f20.0)',iostat=istat) 
     &        shcon(3,1,nmat)
         if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*SPECIFIC GAS CONSTANT%",ier)
            return
         endif
      enddo
!
      return
      end

