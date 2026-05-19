!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine drxs(inpc,textpart,mgrain,nmgrain,
     &  nmat,ncmat_,irstrt,istep,istat,n,iperturb,iline,ipol,
     &  inl,ipoinp,inp,ipoinpc,ier,nstate_) 
!
!     reading the input deck: *METAL POWDER
!
      implicit none
!
      character*1 inpc(*)
      character*132 textpart(16)
!
      integer nmat,istep,istat,ier,
     &  n,key,i,iperturb(2),ncmat_,irstrt(*),iline,ipol,inl,
     &  ipoinp(2,*),inp(3,*),ipoinpc(0:*),id,k,nrhcon(*),nstate_
!
      real*8 mgrain(20,10,*)
!
      iperturb(1)=3
      iperturb(2)=1
      write(*,*) '*INFO reading *METAL POWDER: nonlinear'
      write(*,*) '      geometric effects are turned on'
      write(*,*)
!
      if((istep.gt.0).and.(irstrt(1).ge.0)) then
         write(*,*) '*ERROR reading *drx:'
         write(*,*) '       *drx'
         write(*,*) '  should be placed before all step definitions'
         ier=1
         return
      endif
!
      if(nmat.eq.0) then
         write(*,*) '*ERROR reading *drx:'
         write(*,*) '       *drx'
         write(*,*) '  should bepreceded by a *MATERIAL card'
         ier=1
         return
      endif
!
      do i=2,n
         write(*,*) 
     &        '*WARNING reading *drx:'
         write(*,*) '         parameter not recognized:'
         write(*,*) '         ',
     &        textpart(i)(1:index(textpart(i),' ')-1)
         call inputwarning(inpc,ipoinpc,iline,
     &"DRX%")
      enddo
!
!

      do
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) exit
         ntmat=ntmat+1
         if(ntmat.gt.ntmat_) then
            write(*,*) '*ERROR reading *METAL POWDER:
     &                 increase ntmat_'
            ier=1
            return
         endif
         do i=1,5
            read(textpart(i)(1:20),'(f20.0)',iostat=istat) 
     &            elcon(i+7,ntmat,nmat)
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*METAL POWDER%",ier)
               return
            endif
         enddo
         if(textpart(6)(1:1).ne.' ') then
            read(textpart(6)(1:20),'(f20.0)',iostat=istat)
     &            temperature
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*METAL POWDER%",ier)
               return
            endif
         else
            temperature=0.d0
         endif
         elcon(13,ntmat,nmat)=temperature
      enddo
!
!     interpolating the powder data at the elastic temperature
!     data points
!
      write(*,*) '*INFO: grain size prediction is activated'
      write(*,*)
!
      nelcon(1,nmat)=-70
      nstate_=max(nstate_,16)

      do i=1,nelcon(2,nmat)
         t1l=elcon(0,i,nmat)
         call ident2(elcon(13,1,nmat),t1l,ntmat,ncmat_+1,id)
         if(ntmat.eq.0) then
            continue
         elseif((ntmat.eq.1).or.(id.eq.0)) then
            elcon(3,i,nmat)=elcon(8,1,nmat)
            elcon(4,i,nmat)=elcon(9,1,nmat)
            elcon(5,i,nmat)=elcon(10,1,nmat)
            elcon(6,i,nmat)=elcon(11,1,nmat)
            elcon(7,i,nmat)=elcon(12,1,nmat)
         elseif(id.eq.ntmat) then
            elcon(3,i,nmat)=elcon(8,id,nmat)
            elcon(4,i,nmat)=elcon(9,id,nmat)
            elcon(5,i,nmat)=elcon(10,id,nmat)
            elcon(6,i,nmat)=elcon(11,id,nmat)
            elcon(7,i,nmat)=elcon(12,id,nmat)
         else
            do k=3,7
               elcon(k,i,nmat)=elcon(k+5,id,nmat)+
     &            (elcon(k+5,id+1,nmat)-elcon(k+5,id,nmat))*
     &            (t1l-elcon(13,id,nmat))/
     &            (elcon(13,id+1,nmat)-elcon(13,id,nmat))
            enddo
         endif
         write(*,*) t1l,(elcon(k,i,nmat),k=3,7)
      enddo
    
      return
      end

