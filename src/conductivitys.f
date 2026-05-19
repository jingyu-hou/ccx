!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine conductivitys(inpc,textpart,cocon,ncocon,
     &  nmat,ntmat_,irstrt,istep,istat,n,iline,ipol,inl,ipoinp,inp,
     &  ipoinpc,ier)
!
!     reading the input deck: *CONDUCTIVITY
!
      implicit none
!
      character*1 inpc(*)
      character*132 textpart(16)
!
      integer ncocon(2,*),nmat,ntmat,ntmat_,istep,istat,n,
     &  ipoinpc(0:*),ier,
     &  i,ityp,key,irstrt(*),iline,ipol,inl,ipoinp(2,*),inp(3,*)
!
      real*8 cocon(0:6,ntmat_,*)
!
      ntmat=0
!
      if((istep.gt.0).and.(irstrt(1).ge.0)) then
         write(*,*) 
     &     '*ERROR reading *CONDUCTIVITY: *CONDUCTIVITY should be'
         write(*,*) '  placed before all step definitions'
         ier=1
         return
      endif
!
      if(nmat.eq.0) then
         write(*,*)
     &      '*ERROR reading *CONDUCTIVITY: *CONDUCTIVITY should be'
         write(*,*) '  preceded by a *MATERIAL card'
         ier=1
         return
      endif
!
      ityp=1
!
      do i=2,n
         if(textpart(i)(1:5).eq.'TYPE=') then
            if(textpart(i)(6:8).eq.'ISO') then
               ityp=1
            elseif(textpart(i)(6:10).eq.'ORTHO') then
               ityp=3
            elseif(textpart(i)(6:10).eq.'ANISO') then
               ityp=6
            endif
         else
            write(*,*) 
     &       '*WARNING reading *CONDUCTIVITY: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*CONDUCTIVITY%")
         endif
      enddo
!
      ncocon(1,nmat)=ityp
!
      if(ityp.eq.1) then
         do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) return
            ntmat=ntmat+1
            ncocon(2,nmat)=ntmat
            if(ntmat.gt.ntmat_) then
               write(*,*) 
     &          '*ERROR reading *CONDUCTIVITY: increase ntmat_'
               ier=1
               return
            endif
            do i=1,1
               read(textpart(i)(1:20),'(f20.0)',iostat=istat)
     &                 cocon(i,ntmat,nmat)
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*CONDUCTIVITY%",ier)
                  return
               endif
            enddo
            read(textpart(2)(1:20),'(f20.0)',iostat=istat) 
     &            cocon(0,ntmat,nmat)
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*CONDUCTIVITY%",ier)
               return
            endif
         enddo
      elseif(ityp.eq.3) then
         do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) return
            ntmat=ntmat+1
            ncocon(2,nmat)=ntmat
            if(ntmat.gt.ntmat_) then
               write(*,*) 
     &           '*ERROR reading *CONDUCTIVITY: increase ntmat_'
               ier=1
               return
            endif
            do i=1,3
               read(textpart(i)(1:20),'(f20.0)',iostat=istat)
     &               cocon(i,ntmat,nmat)
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*CONDUCTIVITY%",ier)
                  return
               endif
            enddo
            read(textpart(4)(1:20),'(f20.0)',iostat=istat) 
     &           cocon(0,ntmat,nmat)
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*CONDUCTIVITY%",ier)
               return
            endif
         enddo
      elseif(ityp.eq.6) then
         do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) return
            ntmat=ntmat+1
            ncocon(2,nmat)=ntmat
            if(ntmat.gt.ntmat_) then
               write(*,*) 
     &           '*ERROR reading *CONDUCTIVITY: increase ntmat_'
               ier=1
               return
            endif
            do i=1,6
               read(textpart(i)(1:20),'(f20.0)',iostat=istat)
     &                   cocon(i,ntmat,nmat)
               if(istat.gt.0) then
                  call inputerror(inpc,ipoinpc,iline,
     &                 "*CONDUCTIVITY%",ier)
                  return
               endif
            enddo
            read(textpart(7)(1:20),'(f20.0)',iostat=istat) 
     &            cocon(0,ntmat,nmat)
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*CONDUCTIVITY%",ier)
               return
            endif
         enddo
      endif
!
      return
      end

