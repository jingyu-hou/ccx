!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine gapheatgenerations(inpc,textpart,elcon,nelcon,
     &  imat,ntmat_,ncmat_,irstrt,istep,istat,n,iline,ipol,inl,ipoinp,
     &  inp,ipoinpc,nstate_,ier)
!
!     reading the input deck: *GAP HEAT GENERATION
!
      implicit none
!
      logical user
!
      character*1 inpc(*)
      character*132 textpart(16)
!
      integer nelcon(2,*),imat,ntmat_,istep,istat,ipoinpc(0:*),
     &  n,key,i,ncmat_,irstrt(*),iline,ipol,inl,ipoinp(2,*),inp(3,*),
     &  nstate_,nline,ier
!
      real*8 elcon(0:ncmat_,ntmat_,*)
!
      user=.false.
!
      if((istep.gt.0).and.(irstrt(1).ge.0)) then
         write(*,*) '*ERROR reading *GAP HEAT GENERATION:'
         write(*,*) '       *GAP HEAT GENERATION should be placed'
         write(*,*) '       before all step definitions'
         ier=1
         return
      endif
!
      if(imat.eq.0) then
         write(*,*) '*ERROR reading *GAP HEAT GENERATION:'
         write(*,*) '       *GAP HEAT GENERATION should be preceded'
         write(*,*) '       by a *SURFACE INTERACTION card'
         ier=1
         return
      endif
!
      nstate_=max(nstate_,9)
!
      if(nelcon(1,imat).ne.-51) nelcon(1,imat)=max(nelcon(1,imat),11)
      nelcon(2,imat)=1
!
      do i=2,n
         if(textpart(i)(1:4).eq.'USER') then
            user=.true.
         else
            write(*,*) '*WARNING reading *GAP HEAT GENERATION:'
            write(*,*) '         parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*MATERIAL%")
         endif
      enddo
!
!     no temperature dependence allowed; last line is decisive
!
      nline=0
      do
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) exit
         nline=nline+1
         do i=1,3
            read(textpart(i)(1:20),'(f20.0)',iostat=istat)
     &           elcon(8+i,1,imat)
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*GAP HEAT GENERATION%",ier)
               return
            endif
         enddo
         if(elcon(9,1,imat).le.0.d0) then
            write(*,*) 
     &       '*ERROR reading *GAP HEAT GENERATION: fraction of'
            write(*,*) '       dissipated energy converted into heat'
            write(*,*) '       must be strictly positive'
            call inputerror(inpc,ipoinpc,iline,
     &           "*GAP HEAT GENERATION%",ier)
            return
         endif
         if((elcon(10,1,imat).lt.0.d0).or.
     &      (elcon(10,1,imat).gt.1.d0)) then
            write(*,*) 
     &         '*ERROR reading *GAP HEAT GENERATION: weighting factor'
            write(*,*) '       for the distribution of heat between'
            write(*,*) '       the slave and master surface must '
            write(*,*) '       be contained in [0,1]'
            call inputerror(inpc,ipoinpc,iline,
     &           "*GAP HEAT GENERATION%",ier)
            return
         endif
         elcon(0,1,imat)=0.d0
      enddo
!
      if((.not.user).and.(nline.eq.0)) then
         write(*,*) 
     &        '*ERROR reading *GAP HEAT GENERATION: no data given'
         ier=1
         return
      endif
!
!     user subroutine: vnorm=-0.01
!     no user subroutine: vnorm from input deck: vnorm >= 0.
!                         vnorm calculated from differential velocity: vnorm=-1.
!
      if(user) then
         elcon(11,1,imat)=-0.01d0
      elseif(elcon(11,1,imat).lt.0.d0) then
         elcon(11,1,imat)=-1.d0
      endif
!     
      return
      end

