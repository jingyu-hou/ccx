!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine pk_cdc_cl1(lqd,reynolds,p2p1,beta,kappa,cdc_cl1)
!
      implicit none
!     
      real*8 lqd,reynolds,p2p1,beta,kappa,cdi_noz,cdi_r,cdi_se,
     &     y0,yg,cdc_cl1,rqd,cdqcv_noz,cdqcv_r
!      
      rqd=lqd
!     cd incompresssible nozzle eq. 4a 4b
      call pk_cdi_noz(reynolds,cdi_noz)
!     cdr eq.5
      call pk_cdi_r(rqd,reynolds,beta,cdi_r)
!     cd incompressible sharp edge eq.3
      call pk_cdi_se(reynolds,beta,cdi_se)
!     y0 and yg , eq.15-17 , eq.18
      call pk_y0_yg(p2p1,beta,kappa,y0,yg)
!      
      cdqcv_noz=cdi_noz/(0.0718d0*cdi_noz+0.9282d0)
      cdqcv_r=cdi_r/(0.0718d0*cdi_r+0.9282d0)
!     eq.25 
      cdc_cl1=cdi_r*((cdqcv_noz-cdqcv_r)
     &     /(cdqcv_noz-cdi_se/0.971d0)
     &     *(y0/yg-1d0)+1d0)
!     
      return
!
      end
