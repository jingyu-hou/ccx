!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine userelements(textpart,n,iuel,nuel,inpc,ipoinpc,
     &  iline,ier,ipoinp,inp,inl,ipol)
!
!     reading the input deck: *USER ELEMENT
!
      implicit none
!
      character*1 inpc(*)
      character*132 textpart(16)
!
      integer n,iuel(4,*),nuel,i,j,istat,number,ipoinpc(0:*),iline,
     &  four,nodes,intpoints,maxdof,id,ier,key,ipoinp(2,*),inp(3,*),
     &  inl,ipol
!
      four=4
!
      do i=2,n
         if(textpart(i)(1:6).eq.'TYPE=U') then
            number=ichar(textpart(i)(7:7))*256**3+
     &             ichar(textpart(i)(8:8))*256**2+
     &             ichar(textpart(i)(9:9))*256+
     &             ichar(textpart(i)(10:10))
c            read(textpart(i)(7:16),'(i10)',iostat=istat) number
c            if(istat.gt.0) then
c               call inputerror(inpc,ipoinpc,iline,
c     &                          "*USER ELEMENT%",ier)
c                return
c             endif
         elseif(textpart(i)(1:6).eq.'NODES=') then
            read(textpart(i)(7:16),'(i10)',iostat=istat) nodes
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &                         "*USER ELEMENT%",ier)
               return
            endif
         elseif(textpart(i)(1:18).eq.'INTEGRATIONPOINTS=') then
            read(textpart(i)(19:28),'(i10)',iostat=istat) intpoints
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &                         "*USER ELEMENT%",ier)
               return
            endif
         elseif(textpart(i)(1:7).eq.'MAXDOF=') then
            read(textpart(i)(8:17),'(i10)',iostat=istat) maxdof
            if(istat.gt.0) then
               call inputerror(inpc,ipoinpc,iline,
     &                         "*USER ELEMENT%",ier)
               return
            endif
         endif
      enddo
!
!     check range
!
c      if(number.gt.9999) then
c         write(*,*) '*ERROR reading *USER ELEMENT'
c         write(*,*) '       element number ',number,' exceeds 9999'
c         ier=1
c         return
c      endif
!
      if(intpoints.gt.255) then
         write(*,*) '*ERROR reading *USER ELEMENT'
         write(*,*) '       number of integration points ',intpoints,
     &       ' exceeds 255'
         ier=1
         return
      endif
!
      if(maxdof.gt.255) then
         write(*,*) '*ERROR reading *USER ELEMENT'
         write(*,*) '       highest degree of freedom ',maxdof,
     &       ' exceeds 255'
         ier=1
         return
      endif
!
      if(nodes.gt.255) then
         write(*,*) '*ERROR reading *USER ELEMENT'
         write(*,*) '       number of nodes ',nodes,' exceeds 255'
         ier=1
         return
      endif
!
!     storing the element information in iuel
!
      call nidentk(iuel,number,nuel,id,four)
!
      if(id.gt.0) then
         if(iuel(1,id).eq.number) then
            write(*,*) '*ERROR reading *USER ELEMENT'
            write(*,*) '       element number was already defined'
            ier=1
            return
         endif
      endif
!
      nuel=nuel+1
      do i=nuel,id+2,-1
         do j=1,4
            iuel(j,i)=iuel(j,i-1)
         enddo
      enddo
      iuel(1,id+1)=number
      iuel(2,id+1)=intpoints
      iuel(3,id+1)=maxdof
      iuel(4,id+1)=nodes
!
      return
      end







