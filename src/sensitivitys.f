!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine sensitivitys(inpc,textpart,nmethod,
     &  istep,istat,n,iline,ipol,inl,ipoinp,
     &  inp,tieset,ipoinpc,ntie,tinc,tper,tmin,tmax,tincf,isens,
     &  objectset,ier)
!
!     reading the input deck: *SENSITIVITY
!
      implicit none
!
      logical iread,iwrite
!
      character*1 inpc(*)
      character*81 tieset(3,*),objectset(4,*)
      character*132 textpart(16)
!
      integer nmethod,istep,istat,n,key,i,isens,
     &  iline,ipol,inl,ipoinp(2,*),inp(3,*),
     &  ipoinpc(0:*),ntie,ier
!
      real*8 tinc,tper,tmin,tmax,tincf
!
      if(isens.eq.1) then
         write(*,*) '*ERROR reading *SENSITIVITY:'
         write(*,*) '       no more than one *SENSITIVITY'
         write(*,*) '       is allowed per input deck'
         ier=1
         return
      endif
!
      if(istep.lt.1) then
         write(*,*) '*ERROR reading *SENSITIVITY: *SENSITIVITY can
     &only be used within a STEP'     
         ier=1
         return
      endif
!
      if(istep.lt.2) then
         write(*,*) '*ERROR reading *SENSITIVITY: *SENSITIVITY'
         write(*,*) '       requires a previous *STATIC, *GREEN or'
         write(*,*) '       *FREQUENCY step'
         ier=1
         return
      endif
!
      tinc=0.d0
      tper=0.d0
      tmin=0.d0
      tmax=0.d0
      tincf=0.d0
!
      iwrite=.false.
      iread=.false.
!
      do i=2,n
         if(textpart(i)(1:4).eq.'READ') then
            if(iwrite) then
               write(*,*) '*ERROR reading *SENSITIVITY:'
               write(*,*) '       WRITE and READ are mutually'
               write(*,*) '       exclusive'
               call inputerror(inpc,ipoinpc,iline,
     &              "*SENSITIVITY%",ier)
               return
            endif
            objectset(1,1)(81:81)='R'
            iread=.true.
         elseif(textpart(i)(1:5).eq.'WRITE') then
            if(iread) then
               write(*,*) '*ERROR reading *SENSITIVITY:'
               write(*,*) '       WRITE and READ are mutually'
               write(*,*) '       exclusive'
               call inputerror(inpc,ipoinpc,iline,
     &              "*SENSITIVITY%",ier)
               return
            endif
            objectset(1,1)(81:81)='W'
            iread=.false.
         else
            write(*,*) 
     &        '*WARNING reading *SENSITIVITY: parameter not 
     &recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*SENSITIVITY%")
         endif
      enddo
!
      nmethod=12
!
!     check whether design variables were defined
!
      do i=1,ntie
         if(tieset(1,i)(81:81).eq.'D') exit
      enddo
      if(i.gt.ntie) then
         write(*,*) '*ERROR reading *SENSITIVITY'
         write(*,*) '       no design variables were defined'
         ier=1
         return
      endif
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
!
!
      return
      end

