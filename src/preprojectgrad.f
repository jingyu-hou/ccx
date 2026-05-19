!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine preprojectgrad(vector,ndesi,nodedesi,dgdxglob,nactive,
     &   nobject,nnlconst,ipoacti,nk,rhs,iconst,objectset,xtf)         
!
!     calculates the projected gradient
!
      implicit none
!
      character*81 objectset(4,*)
!
      integer ndesi,nodedesi(*),irow,icol,nactive,nobject,nnlconst,
     &   ipoacti(*),nk,ipos,node,iconst,i
!
      real*8 dgdxglob(2,nk,nobject),vector(ndesi),rhs(*),scalar,dd,
     &   len,xtf(*),brauch,nutz
!
!     initialization of enlarged field dgdxglob and 
!     calculate the second part of xlambd
!
      do irow=1,nk
         dgdxglob(2,irow,nobject)=0.d0
         dgdxglob(1,irow,nobject)=0.d0
      enddo
!     
      do icol=1,nactive
         if(icol.le.nnlconst) then   
            do irow=1,ndesi      
               ipos=ipoacti(icol)   
               node=nodedesi(irow)
               xtf(icol)=xtf(icol)+dgdxglob(2,node,1)
     &                   *dgdxglob(2,node,ipos)
            enddo
         else
            ipos=ipoacti(icol)
            node=nodedesi(ipos)
            xtf(icol)=dgdxglob(2,node,1)
         endif
      enddo
!
      return        
      end




