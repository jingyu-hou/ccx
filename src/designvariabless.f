!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine designvariabless(inpc,textpart,tieset,tietol,istep,
     &       istat,n,iline,ipol,inl,ipoinp,inp,ntie,ntie_,ipoinpc,
     &       set,nset,ier)
!
!     reading the input deck: *DESIGNVARIABLES
!
      implicit none
!
      character*1 inpc(*)
      character*81 tieset(3,*),set(*)
      character*132 textpart(16)
!
      integer istep,istat,n,i,key,ipos,iline,ipol,inl,ipoinp(2,*),
     &  inp(3,*),ntie,ntie_,ipoinpc(0:*),nset,itype,ier
!
      real*8 tietol(3,*)
!
!     Check of correct position in Inputdeck
!
      if(istep.gt.0) then
         write(*,*) '*ERROR reading *DESIGN VARIABLES: *DESIGNVARIABLES'
         write(*,*) ' should be placed before all step definitions'
         ier=1
         return
      endif
!
!     Check of correct number of ties
!
      ntie=ntie+1
      if(ntie.gt.ntie_) then
         write(*,*) '*ERROR reading *DESIGN VARIABLES: increase ntie_'
         ier=1
         return
      endif
!
!     Read in *DESIGNVARIABLES
!
      itype=0
      do i=2,n
         if(textpart(i)(1:5).eq.'TYPE=') then
            read(textpart(i)(6:85),'(a80)',iostat=istat) 
     &           tieset(1,ntie)(1:80)
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*DESIGNVARIABLE%",ier)
               return
            endif
            itype=1
         endif
       enddo  
!
      if(itype.eq.0) then
         write(*,*) 
     &'*ERROR reading *DESIGN VARIABLES: type is lacking'
         call inputerror(inpc,ipoinpc,iline,
     &        "*DESIGNVARIABLE%",ier)
         return
      endif
!
!     Add "D" at the end of the name of the designvariable keyword
!      
      tieset(1,ntie)(81:81)='D' 
!
      if(tieset(1,ntie)(1:10).eq.'COORDINATE') then
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) then
            write(*,*)
     &'*ERROR reading *DESIGN VARIABLES: definition'
            write(*,*) '      is not complete.'
            ier=1
            return
         endif
!
!        Read the name of the design variable node set
!
         tieset(2,ntie)(1:81)=textpart(1)(1:81)
         ipos=index(tieset(2,ntie),' ')
         tieset(2,ntie)(ipos:ipos)='N'
!
!        Check existence of the node set
!
         do i=1,nset
            if(set(i).eq.tieset(2,ntie)) exit
         enddo
         if(i.gt.nset) then
            write(*,*) '*ERROR reading *DESIGN VARIABLES'
            write(*,*) 'node set ',tieset(2,ntie)(1:ipos-1),
     &           'does not exist. Card image:'
            call inputerror(inpc,ipoinpc,iline,
     &           "*DESIGN VARIABLES%",ier)
            return
         endif
      endif
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end



