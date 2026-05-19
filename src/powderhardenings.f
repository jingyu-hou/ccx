!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine powderhardenings(inpc,textpart,mpcon,nmpcon,
     &  nmat,nmpmat_,irstrt,istep,istat,n,iperturb,iline,ipol,
     &  inl,ipoinp,inp,ipoinpc,ier) 
!
!     reading the input deck: *POWDER HARDENDING
!
      implicit none
!
      character*1 inpc(*)
      character*132 textpart(16)
!
      integer nmpcon(2,*),nmat,nmpmat,istep,istat,ier,
     &  n,key,i,iperturb(2),iend,nmpmat_,irstrt(*),iline,ipol,inl,
     &  ipoinp(2,*),inp(3,*),ipoinpc(0:*)
!
      real*8 mpcon(2,nmpmat_,*)
!
      nmpmat=0
      iperturb(1)=3
      iperturb(2)=1
!
      if((istep.gt.0).and.(irstrt(1).ge.0)) then
         write(*,*) '*ERROR reading *POWDER HARDENDING:'
         write(*,*) '       *POWDER HARDENDING'
         write(*,*) '  should be placed in *MATERIAL card'
         ier=1
         return
      endif
!
      if(nmat.eq.0) then
         write(*,*) '*ERROR reading *POWDER HARDENDING:'
         write(*,*) '       *POWDER HARDENDING'
         write(*,*) '  should bepreceded by a *MATERIAL card'
         ier=1
         return
      endif
!
      do i=2,n
         write(*,*) 
     &        '*WARNING reading *POWDER HARDENDING:'
         write(*,*) '         parameter not recognized:'
         write(*,*) '         ',
     &        textpart(i)(1:index(textpart(i),' ')-1)
         call inputwarning(inpc,ipoinpc,iline,
     &"METAL POWDER%")
      enddo
!
      nmpcon(1,nmat)=2
!
      iend=2
      do
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) return
         nmpmat=nmpmat+1
         nmpcon(2,nmat)=nmpmat
         if(nmpmat.gt.400) then
            write(*,*) '*ERROR reading *POWDER HARDENDING:'
            write(*,*) '       increase nmpmat_'
            ier=1
            return
         endif
         do i=1,iend
            read(textpart(i)(1:20),'(f20.0)',iostat=istat) 
     &              mpcon(i,nmpmat,nmat)
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "POWDER HARDENING%",ier)
               return
            endif
         enddo
      enddo
!
      return
      end

