!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine getnewline(inpc,textpart,istat,n,key,iline,
     &  ipol,inl,ipoinp,inp,ipoinpc)
!
      implicit none
!
!     parser for the input file (original order)
!
      character*1 inpc(*)
      character*132 textpart(16)
      character*1320 text
!
      integer istat,n,key,iline,ipol,inl,ipoinp(2,*),inp(3,*),
     &  ipoinpc(0:*),i,j,nentries
!
      parameter(nentries=17)
!
!     reading a new line
!
      if(iline.eq.inp(2,inl)) then
         if(inp(3,inl).eq.0) then
            do
               ipol=ipol+1
               if(ipol.gt.nentries) then
                  istat=-1
                  return
               elseif(ipoinp(1,ipol).ne.0) then
                  exit
               endif
            enddo
            inl=ipoinp(1,ipol)
            iline=inp(1,inl)
         else
            inl=inp(3,inl)
            iline=inp(1,inl)
         endif
      else
         iline=iline+1
      endif
      j=0
      do i=ipoinpc(iline-1)+1,ipoinpc(iline)
         j=j+1
         text(j:j)=inpc(i)
      enddo
      text(j+1:j+1)=' '
!
      istat=0
      key=0
!
!     only free format is supported
!
      if((text(1:1).eq.'*').and.(text(2:2).ne.'*')) then
         key=1
      endif
!
      call splitline(text,textpart,n)
!
      return
      end



