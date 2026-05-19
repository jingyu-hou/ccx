!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine gaps(inpc,textpart,nelcon,nmat,ntmat_,npmat_,
     &        plicon,nplicon,
     &        ncmat_,elcon,matname,irstrt,istep,istat,n,iline,ipol,
     &        inl,ipoinp,inp,nmat_,set,istartset,iendset,ialset,
     &        nset,ielmat,ielorien,ipoinpc,mi,ier)
!
!     reading the input deck: *GAP
!
      implicit none
!
      character*1 inpc(*)
      character*80 matname(*)
      character*81 set(*),elset
      character*132 textpart(16)
!
      integer mi(*),nelcon(2,*),nmat,ntmat_,ntmat,npmat_,npmat,istep,
     &  n,key,i,nplicon(0:ntmat_,*),ncmat_,istat,istartset(*),
     &  iendset(*),irstrt(*),iline,ipol,inl,ipoinp(2,*),inp(3,*),nmat_,
     &  ialset(*),ipos,nset,j,k,ielmat(mi(3),*),ielorien(mi(3),*),
     &  ipoinpc(0:*),ier 
!
      real*8 plicon(0:2*npmat_,ntmat_,*),temperature,
     &  elcon(0:ncmat_,ntmat_,*)
!
      ntmat=0
!
      if((istep.gt.0).and.(irstrt(1).ge.0)) then
         write(*,*) '*ERROR reading *GAP: *GAP should be placed'
         write(*,*) '  before all step definitions'
         ier=1
         return
      endif
!
      nmat=nmat+1
      if(nmat.gt.nmat_) then
         write(*,*) '*ERROR reading *GAP: increase nmat_'
         ier=1
         return
      endif
      matname(nmat)(1:3)='GAP'
      do i=4,80
         matname(nmat)(i:i)=' '
      enddo
!
      do i=2,n
         if(textpart(i)(1:6).eq.'ELSET=') then
            elset=textpart(i)(7:86)
            elset(81:81)=' '
            ipos=index(elset,' ')
            elset(ipos:ipos)='E'
         else
            write(*,*) 
     &        '*WARNING reading *GAP: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*GAP%")
         endif
      enddo
!
!     6 parameters
!
      nelcon(1,nmat)=6
!     
      do
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) exit
         ntmat=ntmat+1
         nelcon(2,nmat)=ntmat
         if(ntmat.gt.ntmat_) then
            write(*,*) '*ERROR reading *GAP: increase ntmat_'
            ier=1
            return
         endif
!     
!     defaults for spring constant (force vs. displacement)
!     and force at infinite clearance
!     
         elcon(5,ntmat,nmat)=1.d12
         elcon(6,ntmat,nmat)=1.d-3
!     
!        reading the initial clearance and the normal direction
!
         do i=1,min(4,n)
            read(textpart(i)(1:20),'(f20.0)',iostat=istat)
     &           elcon(i,ntmat,nmat)
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &                         "*GAP%",ier)
               return
            endif
         enddo
!
!        reading entry 6 and 7 (spring constant and force at
!        infinite clearance)
!
         do i=6,min(7,n)
            read(textpart(i)(1:20),'(f20.0)',iostat=istat)
     &           elcon(i-1,ntmat,nmat)
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &                         "*GAP%",ier)
               return
            endif
         enddo
!
         elcon(0,ntmat,nmat)=0.d0
      enddo
!
      if(ntmat.eq.0) then
         write(*,*) '*ERROR reading *GAP: *GAP card without data'
         ier=1
         return
      endif
      do i=1,nset
         if(set(i).eq.elset) exit
      enddo
      if(i.gt.nset) then
         elset(ipos:ipos)=' '
         write(*,*) '*ERROR reading *GAP: element set ',elset
         write(*,*) '       has not yet been defined. '
         call inputerror(inpc,ipoinpc,iline,
     &        "*GAP%",ier)
         return
      endif
!
!     assigning the elements of the set the appropriate material
!
      do j=istartset(i),iendset(i)
         if(ialset(j).gt.0) then
            ielmat(1,ialset(j))=nmat
            ielorien(1,ialset(j))=0
         else
            k=ialset(j-2)
            do
               k=k-ialset(j)
               if(k.ge.ialset(j-1)) exit
               ielmat(1,k)=nmat
               ielorien(1,k)=0
            enddo
         endif
      enddo
!
      return
      end

