!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine deformationplasticitys(inpc,textpart,elcon,nelcon,
     &  nmat,ntmat_,ncmat_,irstrt,istep,istat,n,iperturb,iline,ipol,
     &  inl,ipoinp,inp,ipoinpc,ier)
!
!     reading the input deck: *DEFORMATION PLASTICITY
!
      implicit none
!
      character*1 inpc(*)
      character*132 textpart(16)
!
      integer nelcon(2,*),nmat,ntmat,ntmat_,istep,istat,ier,
     &  n,key,i,iperturb(2),iend,ncmat_,irstrt(*),iline,ipol,inl,
     &  ipoinp(2,*),inp(3,*),ipoinpc(0:*)
!
      real*8 elcon(0:ncmat_,ntmat_,*)
!
      ntmat=0
      iperturb(1)=3
      iperturb(2)=1
      write(*,*) '*INFO reading *DEFORMATION PLASTICITY: nonlinear'
      write(*,*) '      geometric effects are turned on'
      write(*,*)
!
      if((istep.gt.0).and.(irstrt(1).ge.0)) then
         write(*,*) '*ERROR reading *DEFORMATION PLASTICITY:'
         write(*,*) '       *DEFORMATION PLASTICITY'
         write(*,*) '  should be placed before all step definitions'
         ier=1
         return
      endif
!
      if(nmat.eq.0) then
         write(*,*) '*ERROR reading *DEFORMATION PLASTICITY:'
         write(*,*) '       *DEFORMATION PLASTICITY'
         write(*,*) '  should bepreceded by a *MATERIAL card'
         ier=1
         return
      endif
!
      do i=2,n
         write(*,*) 
     &        '*WARNING reading *DEFORMATION PLASTICITY:'
         write(*,*) '         parameter not recognized:'
         write(*,*) '         ',
     &        textpart(i)(1:index(textpart(i),' ')-1)
         call inputwarning(inpc,ipoinpc,iline,
     &"DEFORMATION PLASTICITY%")
      enddo
!
      nelcon(1,nmat)=-50
!
      iend=5
      do
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) return
         ntmat=ntmat+1
         nelcon(2,nmat)=ntmat
         if(ntmat.gt.ntmat_) then
            write(*,*) '*ERROR reading *DEFORMATION PLASTICITY:'
            write(*,*) '       increase ntmat_'
            ier=1
            return
         endif
         do i=1,iend
            read(textpart(i)(1:20),'(f20.0)',iostat=istat) 
     &              elcon(i,ntmat,nmat)
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "DEFORMATION PLASTICITY%",ier)
               return
            endif
         enddo
         read(textpart(6)(1:20),'(f20.0)',iostat=istat) 
     &              elcon(0,ntmat,nmat)
         if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "DEFORMATION PLASTICITY%",ier)
            return
         endif
      enddo
!
      return
      end

