!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!  
      subroutine amplitudes(inpc,textpart,amname,amta,namta,nam,
     &  nam_,namtot_,irstrt,istep,istat,n,iline,ipol,inl,ipoinp,inp,
     &  ipoinpc,namtot,ier)
!
!     reading the input deck: *AMPLITUDE
!
      implicit none
!
      logical user
!
      character*1 inpc(*)
      character*80 amname(*)
      character*132 textpart(16)
!
      integer namta(3,*),nam,nam_,istep,istat,n,key,i,namtot,
     &  namtot_,irstrt(*),iline,ipol,inl,ipoinp(2,*),inp(3,*),ipos,
     &  ipoinpc(0:*),ier
!
      real*8 amta(2,*),x,y,shiftx,shifty
!
      user=.false.
!
      shiftx=0.d0
      shifty=0.d0
!
      if((istep.gt.0).and.(irstrt(1).ge.0)) then
         write(*,*) '*ERROR reading *AMPLITUDE: *AMPLITUDE should be'
         write(*,*) '  placed before all step definitions'
         ier=1
         return
      endif
!
      nam=nam+1
      if(nam.gt.nam_) then
         write(*,*) '*ERROR reading *AMPLITUDE: increase nam_'
         ier=1
         return
      endif
      namta(3,nam)=nam
      amname(nam)='
     &                           '
!
      do i=2,n
         if(textpart(i)(1:5).eq.'NAME=') then
            amname(nam)=textpart(i)(6:85)
            if(textpart(i)(86:86).ne.' ') then
               write(*,*) '*ERROR reading *AMPLITUDE: amplitude'
               write(*,*) '      name is too long'
               write(*,*) '      (more than 80 characters)'
               write(*,*) '      amplitude name:',textpart(i)(1:132)
               ier=1
               return
            endif
         elseif(textpart(i)(1:14).eq.'TIME=TOTALTIME') then
            namta(3,nam)=-nam
         elseif(textpart(i)(1:4).eq.'USER') then
            namta(1,nam)=0
            namta(2,nam)=0
            user=.true.
         elseif(textpart(i)(1:18).eq.'DEFINITION=TABULAR') then
            cycle
         elseif(textpart(i)(1:14).eq.'VALUE=RELATIVE') then
            cycle
         elseif(textpart(i)(1:6).eq.'SHIFTX') then
               read(textpart(i)(8:27),'(f20.0)',iostat=istat) shiftx
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*AMPLITUDE%",ier)
                  return
               endif
         elseif(textpart(i)(1:6).eq.'SHIFTY') then
               read(textpart(i)(8:27),'(f20.0)',iostat=istat) shifty
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*AMPLITUDE%",ier)
                  return
               endif
         else
            write(*,*) 
     &        '*WARNING reading *AMPLITUDE: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*AMPLITUDE%")
         endif
      enddo
!
      if(amname(nam).eq.'                                               
     &                                 ') then
         write(*,*) '*ERROR reading *AMPLITUDE: Amplitude has no name'
         call inputerror(inpc,ipoinpc,iline,
     &        "*AMPLITUDE%",ier)
         return
      endif
!
      if(.not.user) then
         namta(1,nam)=namtot+1
      endif
!
      do
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) exit
         do i=1,4
            if(textpart(2*i-1)(1:1).ne.' ') then  
               namtot=namtot+1
               if(namtot.gt.namtot_) then
                  write(*,*) 
     &               '*ERROR reading *AMPLITUDE: increase namtot_'
                  ier=1
                  return
               endif
               read(textpart(2*i-1),'(f20.0)',iostat=istat) x
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*AMPLITUDE%",ier)
                  return
               endif
               read(textpart(2*i),'(f20.0)',iostat=istat) y
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*AMPLITUDE%",ier)
                  return
               endif
               amta(1,namtot)=x+shiftx
               amta(2,namtot)=y+shifty
               namta(2,nam)=namtot
            else
               exit
            endif
         enddo
      enddo
!
      if(namta(1,nam).gt.namta(2,nam)) then
         ipos=index(amname(nam),' ')
         write(*,*) 
     &      '*WARNING reading *AMPLITUDE: *AMPLITUDE definition ',
     &       amname(nam)(1:ipos-1) 
         write(*,*) '         has no data points'
         nam=nam-1
c      else
c         call reorderampl(amname,namta,nam)
      endif
!
      return
      end

