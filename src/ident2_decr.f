!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
!                                                                             
!     identifies the position id of px in an ordered array
!     x of real numbers; The numbers in x are at positions
!     1, 1+ninc, 1+2*ninc, 1+3*ninc... up to 1+(n-1)*ninc
!  
!     id is such that x(id+1).le.px and x(id).gt.px
!     ---2021.3.25 by Ye
!                                                                             
      subroutine ident2_decr(x,px,n,ninc,id)
      implicit none   
! 
      integer n,id,n2,m,ninc
!
      real*8 x(n*ninc),px
!
      intent(in) x,px,n,ninc
!
      intent(out) id
!
      id=0
      if(n.eq.0) return
      n2=n+1
      do                                                            
         m=(n2+id)/2
         if(px.le.x(1+ninc*(m-1))) then
            id=m      
         else
            n2=m            
         endif
         if((n2-id).eq.1) return
      enddo
      end
                                                                               
