!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine bodyforce(cbody,ibody,ipobody,nbody,set,istartset,
     &  iendset,ialset,inewton,nset,ifreebody,k)
!
!     assigns the body forces to the elements by use of field ipobody
!
      implicit none
!
      character*81 cbody(*),elset,set(*)
!
      integer ibody(3,*),ipobody(2,*),i,j,l,istartset(*),nbody,
     &  iendset(*),ialset(*),kindofbodyforce,inewton,nset,istat,
     &  ifreebody,k,index
!
      elset=cbody(k)
      kindofbodyforce=ibody(1,k)
      if(kindofbodyforce.eq.3) inewton=1
!
!     check whether element number or set name
!
      read(elset,'(i21)',iostat=istat) l
      if(istat.eq.0) then
         if(ipobody(1,l).eq.0) then
            ipobody(1,l)=k
         else
!
            index=l
            do
               if(ipobody(1,index).eq.k) exit
               if(ipobody(2,index).eq.0) then
                  ipobody(2,index)=ifreebody
                  ipobody(1,ifreebody)=k
                  ipobody(2,ifreebody)=0
                  ifreebody=ifreebody+1
                  exit
               endif
               index=ipobody(2,index)
            enddo
         endif
         return
      endif
!
!     set name
!
      do i=1,nset
         if(set(i).eq.elset) exit
      enddo
!     
      do j=istartset(i),iendset(i)
         if(ialset(j).gt.0) then
            l=ialset(j)
            if(ipobody(1,l).eq.0) then
               ipobody(1,l)=k
            else
!
               index=l
               do
                  if(ipobody(1,index).eq.k) exit
                  if(ipobody(2,index).eq.0) then
                     ipobody(2,index)=ifreebody
                     ipobody(1,ifreebody)=k
                     ipobody(2,ifreebody)=0
                     ifreebody=ifreebody+1
                     exit
                  endif
                  index=ipobody(2,index)
               enddo
            endif
         else
            l=ialset(j-2)
            do
               l=l-ialset(j)
               if(l.ge.ialset(j-1)) exit
               if(ipobody(1,l).eq.0) then
                  ipobody(1,l)=k
               else
!
                  index=l
                  do
                     if(ipobody(1,index).eq.k) exit
                     if(ipobody(2,index).eq.0) then
                        ipobody(2,index)=ifreebody
                        ipobody(1,ifreebody)=k
                        ipobody(2,ifreebody)=0
                        ifreebody=ifreebody+1
                        exit
                     endif
                     index=ipobody(2,index)
                  enddo
               endif
            enddo
         endif
      enddo
!     
      return
      end

