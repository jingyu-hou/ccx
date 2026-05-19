!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine outputs(inpc,textpart,jout,itpamp,istep,istat,n,iline,
     &  ipol,inl,ipoinp,inp,ipoinpc,ier)
!
!     reading the *OUTPUT card in the input deck
!
      implicit none
!
      character*1 inpc(*)
      character*132 textpart(16)
!
      integer istep,istat,n,key,ii,jout(2),joutl,iline,ipol,inl,
     &  ipoinp(2,*),inp(3,*),ipoinpc(0:*),itpamp,ier
!
      if(istep.lt.1) then
         write(*,*) '*ERROR reading *OUTPUT'
         write(*,*) '       *OUTPUT'
         write(*,*) '       should only be used within a *STEP' 
         write(*,*) '       definition'
         ier=1
         return
      endif
!
      do ii=2,n
        if(textpart(ii)(1:10).eq.'FREQUENCY=') then
           read(textpart(ii)(11:20),'(i10)',iostat=istat) joutl
           if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*OUTPUT %",ier)
              return
           endif
           if(joutl.eq.0) then
              do
                 call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &                inl,ipoinp,inp,ipoinpc)
                 if((key.eq.1).or.(istat.lt.0)) return
              enddo
           endif
           if(joutl.gt.0) then
              jout(1)=joutl
              itpamp=0
           endif
        elseif(textpart(ii)(1:11).eq.'FREQUENCYF=') then
           read(textpart(ii)(12:21),'(i10)',iostat=istat) joutl
           if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*OUTPUT %",ier)
              return
           endif
           if(joutl.eq.0) then
              do
                 call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &                inl,ipoinp,inp,ipoinpc)
                 if((key.eq.1).or.(istat.lt.0)) return
              enddo
           endif
           if(joutl.gt.0) then
              jout(2)=joutl
              itpamp=0
           endif
        else
            write(*,*) 
     &             '*WARNING reading *NODE FILE or *EL FILE:' 
            write(*,*) '         parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(ii)(1:index(textpart(ii),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*OUTPUT %")
        endif
      enddo
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end
