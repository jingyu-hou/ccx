!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine mpcrem(i,mpcfree,nodempc,nmpc,ikmpc,ilmpc,labmpc,
     &  coefmpc,ipompc)
!
!     removes multiple point constraint i
!
      implicit none
!
      character*20 labmpc(*)
!
      integer nodempc(3,*),node,nmpc,i,j,index,mpcfree,mpcfreeold,
     &  ikmpc(*),ilmpc(*),idof,id,ipompc(*),idir
!
      real*8 coefmpc(*)
!
      mpcfreeold=mpcfree
      index=ipompc(i)
      ipompc(i)=0
      node=nodempc(1,index)
      idir=nodempc(2,index)
      idof=8*(node-1)+idir
      call nident(ikmpc,idof,nmpc,id)
c      mpcfree=nodempc(3,index)
      mpcfree=index
!
!     removing the MPC from fields nodempc and coefmpc
!
      do
         nodempc(1,index)=0
         nodempc(2,index)=0
         coefmpc(index)=0
         if(nodempc(3,index).ne.0) then
            index=nodempc(3,index)
         else
            nodempc(3,index)=mpcfreeold
            exit
         endif
      enddo
!
!     decrementing nmpc
!
      nmpc=nmpc-1
!
!     shifting fields ikmpc,ilmpc
!
      do j=id,nmpc
         ikmpc(j)=ikmpc(j+1)
         ilmpc(j)=ilmpc(j+1)
      enddo
      ikmpc(nmpc+1)=0
      ilmpc(nmpc+1)=0
!
!     shifting fields ipompc,labmpc
!
      do j=i,nmpc
         ipompc(j)=ipompc(j+1)
         labmpc(j)=labmpc(j+1)
      enddo
      ipompc(nmpc+1)=0
!
!     updating ilmpc
!
      do j=1,nmpc
         if(ilmpc(j).gt.i) ilmpc(j)=ilmpc(j)-1
      enddo
!
      return
      end

