!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
   
      subroutine pk_y0_yg(p2p1,beta,kappa,y0,yg)
!
      implicit none
!
      real*8 p2p1,beta,kappa,y0,yg,pcrit
!
!     adiabatic expansion factor y0 measured (eq.15-17) 
!
!     author: Yannick Muller
!     
      pcrit=(2.d0/(kappa+1.d0))**(kappa/(kappa-1.d0))

      if(p2p1.ge.0.63d0) then
         y0=1d0-(0.41d0+0.35d0*beta**4.d0)/kappa*(1.d0-p2p1)
      else
         y0=1d0-(0.41d0+0.35d0*beta**4.d0)/kappa*(1.d0-0.63d0)
     &        -(0.3475d0+0.1207d0*beta**2.d0-0.3177d0*beta**4.d0)
     &        *(0.63d0-p2p1)
!         
      endif
!    
!     adiabatic expension factor yg isentropic eq 18
!
      if(p2p1.ge.1d0) then
         yg=1.d0
!         
      elseif (p2p1.ge.pcrit) then
         yg=p2p1**(1.d0/kappa)*dsqrt(kappa/(kappa-1.d0)
     &        *(1.d0-p2p1**((kappa-1.d0)/kappa)))/dsqrt(1.d0-p2p1)
!      
      else
!     critical pressure ratio
         yg=(2.d0/(kappa+1.d0))**(1.d0/(kappa-1.d0))
     &       *dsqrt(kappa/(kappa+1.d0))/dsqrt(1.d0-p2p1)
      endif
!     
      return
!      
      end
