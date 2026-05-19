!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine materials(inpc,textpart,matname,nmat,nmat_,
     &  irstrt,istep,istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc,
     &  imat,ier)
!
!     reading the input deck: *MATERIAL
!
      implicit none
!
      character*1 inpc(*)
      character*80 matname(*)
      character*132 textpart(16)
!
      integer nmat,nmat_,istep,istat,n,key,i,irstrt(*),iline,ipol,inl,
     &  ipoinp(2,*),inp(3,*),ipoinpc(0:*),imat,ier
!
      if((istep.gt.0).and.(irstrt(1).ge.0)) then
         write(*,*) 
     &       '*ERROR reading *MATERIAL: *MATERIAL should be placed'
         write(*,*) '  before all step definitions'
         ier=1
         return
      endif
!
      nmat=nmat+1
      if(nmat.gt.nmat_) then
         write(*,*) '*ERROR reading *MATERIAL: increase nmat_'
         ier=1
         return
      endif
!
      imat=nmat
!
      do i=2,n
         if(textpart(i)(1:5).eq.'NAME=') then
            matname(nmat)=textpart(i)(6:85)
            if(textpart(i)(86:86).ne.' ') then
               write(*,*) 
     &             '*ERROR reading *MATERIAL: material name too long'
               write(*,*) '       (more than 80 characters)'
               write(*,*) '       material name:',textpart(i)(1:132)
               ier=1
               return
            endif
            exit
         else
            write(*,*) 
     &        '*WARNING reading *MATERIAL: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*MATERIAL%")
         endif
      enddo
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end

