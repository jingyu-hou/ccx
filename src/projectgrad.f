!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine projectgrad(vector,ndesi,nodedesi,dgdxglob,nactive,
     &   nobject,nnlconst,ipoacti,nk,rhs,iconst,objectset,xlambd,xtf,
     &   objnorm)         
!
!     calculates the projected gradient, the lagrange multipliers and 
!     the correction term due to nonlinear constraints
!
      implicit none
!
      character*81 objectset(4,*)
!
      integer ndesi,nodedesi(*),irow,icol,nactive,nobject,nnlconst,
     &   ipoacti(*),nk,ipos,node,iconst,i
!
      real*8 dgdxglob(2,nk,nobject),vector(ndesi),rhs(*),scalar,dd,
     &   len,xlambd(*),xtf(*),objnorm(nobject)
!
!     calculation of the vector
!
      do irow=1,ndesi
         do icol=1,nactive
            if(icol.le.nnlconst) then
               ipos=ipoacti(icol)
               node=nodedesi(irow)
               vector(irow)=vector(irow)
     &                      +rhs(icol)*dgdxglob(2,node,ipos)
            else
               if(ipoacti(icol).eq.irow) then
                  vector(irow)=vector(irow)+rhs(icol)
               endif
            endif
         enddo
      enddo
!
!     calculation of the scalar value
!
      scalar=0.d0
      if(iconst.le.nnlconst) then
         do irow=1,ndesi
            ipos=ipoacti(iconst)      
            node=nodedesi(irow)
            scalar=scalar+dgdxglob(2,node,1)*dgdxglob(2,node,ipos)
         enddo
      else
         node=nodedesi(ipoacti(iconst))
         scalar=dgdxglob(2,node,1)
      endif            
!
!     multiplication of scalar and vector
!
      do irow=1,ndesi
         node=nodedesi(irow)
         dgdxglob(2,node,nobject)=dgdxglob(2,node,nobject)
     &                            +vector(irow)*scalar
      enddo
!
!     calculation of lagrange multiplier
! 
      do icol=1,nactive
         xlambd(icol)=xlambd(icol)+(-1)*rhs(icol)*xtf(iconst)
      enddo
!
!     calculation of correction term
!
      if(iconst.le.nnlconst) then
         do irow=1,ndesi
            node=nodedesi(irow)
             if(abs(objnorm(iconst)).gt.0.d0) then
                dgdxglob(1,node,nobject)=dgdxglob(1,node,nobject)
     &                   +vector(irow)*objnorm(iconst)
     &                   /abs(objnorm(iconst))
             else
                dgdxglob(1,node,nobject)=dgdxglob(1,node,nobject)
     &                   +vector(irow)*objnorm(iconst)       
             endif
         enddo
      endif
!
      return        
      end
