!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine applympc(nface,ielfa,is,ie,ifabou,ipompc,vfa,coefmpc,
     &  nodempc,ipnei,neifa,labmpc,xbounact,nactdoh,ifaext,nfaext)
!
!     applies MPC's to the faces
!
      implicit none
!
      character*20 labmpc(*)
!
      integer i,j,nface,ielfa(4,*),ipointer,is,ie,ifabou(*),mpc,
     &  ipompc(*),index,iel,iface,nodempc(3,*),ipnei(*),neifa(*),
     &  nactdoh(*),ielorig,ifaext(*),nfaext,k
!
      real*8 coefmpc(*),denominator,vfa(0:7,*),sum,xbounact(*),
     &  coefnorm
!
!
!
      do k=1,nfaext
         i=ifaext(k)
         if(ielfa(2,i).ge.0) cycle
         ipointer=-ielfa(2,i)
         do j=is,ie
            if(ifabou(ipointer+j).ge.0) cycle
            mpc=-ifabou(ipointer+j)
!
            index=ipompc(mpc)
            sum=0.d0
            coefnorm=0.d0
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
     &                 *vfa(nodempc(2,index),neifa(ipnei(iel)+iface))
                  coefnorm=coefnorm+coefmpc(index)**2
               endif
               index=nodempc(3,index)
            enddo
!
!           distribute the sum across all terms which are not
!           fixed by boundary conditions
!
            index=ipompc(mpc)
            do
               if(index.eq.0) exit
               if(nodempc(1,index).gt.0) then
                  ielorig=int(nodempc(1,index)/10.d0)
                  iel=nactdoh(ielorig)
                  iface=nodempc(1,index)-10*ielorig
                  vfa(nodempc(2,index),neifa(ipnei(iel)+iface))=
     &               vfa(nodempc(2,index),neifa(ipnei(iel)+iface))-
     &               sum*coefmpc(index)/coefnorm
               endif
               index=nodempc(3,index)
            enddo
         enddo
      enddo
!     
      return
      end
