!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine nodalthicknesss(inpc,textpart,set,istartset,iendset,
     &  ialset,nset,thickn,nk,istep,istat,n,iline,ipol,inl,ipoinp,
     &  inp,iaxial,ipoinpc,ier)
!
!     reading the input deck: *NODAL THICKNESS
!
      implicit none
!
      character*1 inpc(*)
      character*81 set(*),noset
      character*132 textpart(16)
!
      integer istartset(*),iendset(*),ialset(*),nset,nk,istep,istat,n,
     &  key,i,j,k,l,ipos,iline,ipol,inl,ipoinp(2,*),inp(3,*),iaxial,
     &  ipoinpc(0:*),ier
!
      real*8 thickn(2,*),thickness1,thickness2
!
      if(istep.gt.0) then
         write(*,*) '*ERROR reading *NODAL THICKNESS: *NODAL THICKNESS'
         write(*,*) '      should be placed before all step definitions'
         ier=1
         return
      endif
!
      do i=2,n
         write(*,*) 
     &    '*WARNING reading *NODAL THICKNESS: parameter not recognized:'
         write(*,*) '         ',
     &        textpart(i)(1:index(textpart(i),' ')-1)
         call inputwarning(inpc,ipoinpc,iline,
     &"*NODAL THICKNESS%")
      enddo
!
      do
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) return
         read(textpart(2)(1:20),'(f20.0)',iostat=istat) thickness1
         if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*NODAL THICKNESS%",ier)
            return
         endif
         if(iaxial.eq.180) thickness1=thickness1/iaxial
         if(n.eq.2) then
            thickness2=0.d0
         else
            read(textpart(2)(1:20),'(f20.0)',iostat=istat) thickness2
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &              "*NODAL THICKNESS%",ier)
               return
            endif
         endif
         read(textpart(1)(1:10),'(i10)',iostat=istat) l
         if(istat.eq.0) then
            thickn(1,l)=thickness1
            thickn(2,l)=thickness2
         else
            read(textpart(1)(1:80),'(a80)',iostat=istat) noset
            noset(81:81)=' '
            ipos=index(noset,' ')
            noset(ipos:ipos)='N'
            do i=1,nset
               if(set(i).eq.noset) exit
            enddo
            if(i.gt.nset) then
               noset(ipos:ipos)=' '
               write(*,*) '*ERROR reading *NODAL THICKNESS: node set ',
     &                noset
               write(*,*) '  has not yet been defined. '
               call inputerror(inpc,ipoinpc,iline,
     &              "*NODAL THICKNESS%",ier)
               return
            endif
            do j=istartset(i),iendset(i)
               if(ialset(j).gt.0) then
                  thickn(1,ialset(j))=thickness1
                  thickn(2,ialset(j))=thickness2
               else
                  k=ialset(j-2)
                  do
                     k=k-ialset(j)
                     if(k.ge.ialset(j-1)) exit
                     thickn(1,k)=thickness1
                     thickn(2,k)=thickness2
                  enddo
               endif
            enddo
         endif
      enddo
!
      return
      end

