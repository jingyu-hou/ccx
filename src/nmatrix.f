!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine nmatrix(ad,au,jqs,irows,ndesi,nodedesi,dgdxglob,
     &   nactive,nobject,nnlconst,ipoacti,nk)         
!
!     calculates the values of the expression: N^(T)N
!
      implicit none
!
      integer jqs(*),irows(*),ndesi,nodedesi(*),idof,i,j,jdof,
     &   nactive,nobject,nnlconst,ipos,jpos,ipoacti(*),nk,node 
!
      real*8 ad(*),au(*),dgdxglob(2,nk,nobject)
!
      do idof=1,nactive
         if(idof.le.nnlconst) then
            ipos=ipoacti(idof)
            do i=1,ndesi
               node=nodedesi(i)
               ad(idof)=ad(idof)+dgdxglob(2,node,ipos)**2
            enddo
            do i=jqs(idof),jqs(idof+1)-1
               jdof=irows(i)
               if(jdof.le.nnlconst) then
                  jpos=ipoacti(i)
                  do j=1,ndesi
                     node=nodedesi(j)
                     au(i)=au(i)+dgdxglob(2,node,ipos)
     &                          *dgdxglob(2,node,jpos)
                  enddo
               else
                  node=nodedesi(ipoacti(i))
                  au(i)=dgdxglob(2,node,ipos)
               endif
            enddo
         else
            ad(idof)=1 
         endif
      enddo
!
      return        
      end




