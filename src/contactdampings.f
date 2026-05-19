!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine contactdampings(inpc,textpart,elcon,nelcon,
     &  nmat,ntmat_,ncmat_,irstrt,istep,istat,n,iline,ipol,inl,ipoinp,
     &  inp,ipoinpc,imat,ier)
!
!     reading the input deck: *CONTACT DAMPING
!
      implicit none
!
      character*1 inpc(*)
      character*132 textpart(16)
!
      integer nelcon(2,*),nmat,ntmat_,istep,istat,ipoinpc(0:*),
     &  n,key,i,ncmat_,irstrt(*),iline,ipol,inl,ipoinp(2,*),inp(3,*),
     &  imat,ier
!
      real*8 elcon(0:ncmat_,ntmat_,*)
!
      if((istep.gt.0).and.(irstrt(1).ge.0)) then
         write(*,*) '*ERROR reading *CONTACT DAMPING:'
         write(*,*) '       *CONTACT DAMPING should be placed'
         write(*,*) '       before all step definitions'
         ier=1
         return
      endif
!
      if(imat.eq.0) then
         write(*,*) '*ERROR reading *CONTACT DAMPING:'
         write(*,*) '       *CONTACT DAMPING should be preceded'
         write(*,*) '       by a *SURFACE INTERACTION card'
         ier=1
         return
      endif
!
!     default: no tangential damping
!
      elcon(8,1,imat)=0.d0
!
      do i=2,n
         if(textpart(i)(1:16).eq.'TANGENTFRACTION=') then
            read(textpart(i)(17:36),'(f20.0)',iostat=istat) 
     &               elcon(8,1,imat)
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*CONTACT DAMPING%",ier)
               return
            endif
         else
            write(*,*) 
     &   '*WARNING reading *CONTACT DAMPING: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*CONTACT DAMPING%")
         endif
      enddo
!
      nelcon(1,imat)=max(nelcon(1,imat),8)
      nelcon(2,imat)=1
!
!     no temperature dependence allowed; last line is decisive
!
      do
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) return
         read(textpart(1)(1:20),'(f20.0)',iostat=istat)
     &        elcon(5,1,imat)
         if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &        "*CONTACT DAMPING%",ier)
            return
         endif
         elcon(0,1,imat)=0.d0
      enddo
!     
      return
      end

