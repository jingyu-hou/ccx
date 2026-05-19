!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine surfaceinteractions(inpc,textpart,matname,nmat,nmat_,
     &  irstrt,istep,istat,n,iline,ipol,inl,ipoinp,inp,nrhcon,ipoinpc,
     &  imat,ier,npair,pair2mat,pairmasterslave)
!
!     reading the input deck: *SURFACE INTERACTION
!
      implicit none
!
      character*1 inpc(*)
      character*80 matname(*)
      character*132 textpart(16)
      character*80 pairmasterslave(2,nmat_)
!
      integer nmat,nmat_,istep,istat,n,key,i,irstrt(*),iline,ipol,inl,
     &  ipoinp(2,*),inp(3,*),nrhcon(*),ipoinpc(0:*),imat,ier,
     &  npair,pair2mat(nmat_)
!
      if((istep.gt.0).and.(irstrt(1).ge.0)) then
         write(*,*) '*ERROR reading *SURFACE INTERACTION:'
         write(*,*) '       *SURFACE INTERACTION should be placed'
         write(*,*) '       before all step definitions'
         ier=1
         return
      endif
!
      nmat=nmat+1
      npair=npair+1
      pair2mat(npair)=nmat
      pairmasterslave(1,npair)="NONEMASTER"
      pairmasterslave(2,npair)="NONSLAVE"
      if(nmat.gt.nmat_) then
         write(*,*) 
     &       '*ERROR reading *SURFACE INTERACTION: increase nmat_'
         ier=1
         return
      endif
      imat=nmat
!
      do i=2,n
         if(textpart(i)(1:5).eq.'NAME=') then
            matname(nmat)=textpart(i)(6:85)
            if(textpart(i)(86:86).ne.' ') then
               write(*,*) '*ERROR reading *SURFACE INTERACTION:'
               write(*,*) '       name too long'
               write(*,*) '       (more than 80 characters)'
               write(*,*) '       interaction name:',textpart(i)(1:132)
               ier=1
               return
            endif
c            exit
         elseif(textpart(i)(1:7).eq.'MASTER=') then
            pairmasterslave(1,npair)=textpart(i)(8:85)
         elseif(textpart(i)(1:6).eq.'SLAVE=') then
            pairmasterslave(2,npair)=textpart(i)(7:85)
         else
            write(*,*) '*WARNING reading *SURFACE INTERACTION:'
            write(*,*) '         parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*SURFACE INTERACTION%")
         endif
      enddo
!
!     a fictitious nonzero number of density values is stored in nrhcon
!     for contact calculations in which all materials are required to 
!     have a density assigned (e.g. dynamic calculations). This is needed
!     since a surface interaction is internally treated as material
!
      nrhcon(nmat)=-1
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end

