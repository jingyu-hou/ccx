!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine basemotions(inpc,textpart,amname,nam,ibasemotion,
     &  xboun,ndirboun,iamboun,typeboun,nboun,istep,istat,n,iline,ipol,
     &  inl,ipoinp,inp,ipoinpc,iamplitudedefault,ier)
!
!     reading the input deck: *BASE MOTION
!
      implicit none
!
      character*1 typeboun(*),type,inpc(*)
      character*80 amname(*),amplitude
      character*132 textpart(16)
!
      integer iamplitude,idof,i,j,n,nam,ibasemotion,iamboun(*),nboun,
     &  ndirboun(*),istep,istat,iline,ipol,inl,ipoinp(2,*),inp(3,*),
     &  key,ipoinpc(0:*),iamplitudedefault,ier
!
      real*8 xboun(*)
!
      type='A'
      iamplitude=iamplitudedefault
      idof=0
!
      if(istep.lt.1) then
         write(*,*) '*ERROR reading *BASE MOTION:'
         write(*,*) '       *BASE MOTION should only be used'
         write(*,*) '       within a STEP'
         ier=1
         return
      endif
!
      do i=2,n
         if(textpart(i)(1:4).eq.'DOF=') then
            read(textpart(i)(5:14),'(i10)',iostat=istat) idof
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*BASE MOTION%",ier)
               return
            endif
         elseif(textpart(i)(1:10).eq.'AMPLITUDE=') then
            read(textpart(i)(11:90),'(a80)') amplitude
            do j=1,nam
               if(amname(j).eq.amplitude) then
                  iamplitude=j
                  exit
               endif
            enddo
            if(j.gt.nam) then
               write(*,*) '*ERROR reading *BASE MOTION:'
               write(*,*) '       nonexistent amplitude'
               write(*,*) '  '
               call inputerror(inpc,ipoinpc,iline,
     &              "*BASE MOTION%",ier)
               return
            endif
            iamplitude=j
         elseif(textpart(i)(1:5).eq.'TYPE=') then
            if(textpart(i)(6:17).eq.'DISPLACEMENT') then
               type='B'
            elseif(textpart(i)(6:17).eq.'ACCELERATION') then
               type='A'
            else
               write(*,*) '*ERROR reading *BASE MOTION:'
               write(*,*) '       invalid TYPE'
               call inputerror(inpc,ipoinpc,iline,
     &              "*BASE MOTION%",ier)
               return
            endif
         else
            write(*,*) 
     &        '*WARNING reading *BASE MOTION: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*BASE MOTION%")
         endif
      enddo
!
      if(idof.eq.0) then
         write(*,*) '*ERROR reading *BASE MOTION'
         write(*,*) '       no degree of freedom specified'
         ier=1
         return
      elseif(iamplitude.eq.0) then
         write(*,*) '*ERROR reading *BASE MOTION'
         write(*,*) '       no amplitude specified'
         ier=1
         return
      endif
!
      if(ibasemotion.eq.0) then
!
!        no previous *BASE MOTION within the actual step
!
         ibasemotion=1
         do i=1,nboun
            if(ndirboun(i).eq.idof) then
               xboun(i)=1.d0
               iamboun(i)=iamplitude
               typeboun(i)=type
            else
               xboun(i)=0.d0
            endif
         enddo
      else
!
!        previous *BASE MOTION within the actual step
!
         do i=1,nboun
            if(ndirboun(i).eq.idof) then
               xboun(i)=1.d0
               iamboun(i)=iamplitude
               typeboun(i)=type
            endif
         enddo
      endif
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end

