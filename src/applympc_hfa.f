!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine applympc_hfa(nface,ielfa,is,ie,ifabou,ipompc,hfa,
     &  coefmpc,nodempc,ipnei,neifa,labmpc,xbounact,nactdoh)
!
!     applies MPC's to the faces
!
      implicit none
!
      character*20 labmpc(*)
!
      integer i,j,nface,ielfa(4,*),ipointer,is,ie,ifabou(*),mpc,
     &  ipompc(*),index,iel,iface,nodempc(3,*),ipnei(*),neifa(*),
     &  nactdoh(*),ielorig
!
      real*8 coefmpc(*),denominator,hfa(3,*),sum,xbounact(*)
!
      do i=1,nface
         if(ielfa(2,i).ge.0) cycle
         ipointer=-ielfa(2,i)
         do j=is,ie
            if(ifabou(ipointer+j).ge.0) cycle
            mpc=-ifabou(ipointer+j)
            index=ipompc(mpc)
            denominator=coefmpc(index)
            sum=0.d0
            index=nodempc(3,index)
            do
               if(index.eq.0) exit
               if(nodempc(1,index).lt.0) then
!
!                 a negative number refers to a boundary
!                 condition (fields nodeboun, ndirboun..)
!                 resulting from a SPC in local coordinates
!                  
                  sum=sum+coefmpc(index)*xbounact(-nodempc(1,index))
               else
!
!                 face term
!
                  ielorig=int(nodempc(1,index)/10.d0)
                  iel=nactdoh(ielorig)
                  iface=nodempc(1,index)-10*ielorig
                  sum=sum+coefmpc(index)
     &                 *hfa(nodempc(2,index),neifa(ipnei(iel)+iface))
               endif
               index=nodempc(3,index)
            enddo
            hfa(j,i)=-sum/denominator
         enddo
      enddo
!     
      return
      end
