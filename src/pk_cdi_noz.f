!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
     
!cd  incompressible for ASME nozzles eq 4a 4b   
!
!     author: Yannick Muller
!
      subroutine pk_cdi_noz(reynolds,cdi_noz)
!
      implicit none
!      
      real*8 reynolds,cdi_noz,ln_reynolds,cdi_noz_lr,
     &     cdi_noz_hr,e,reynolds_cor
!
      if (reynolds.lt.40000d0) then
!
! formerly pk_cdi_noz_lr : for low Reynolds nsumber
!         
         if (reynolds.eq.0d0) then
            reynolds_cor=1.d0
         else
            reynolds_cor=reynolds
         endif
         e=2.718281828459045d0
         ln_reynolds=log(reynolds_cor)/log(e)
!     
         cdi_noz_lr=0.19436d0+0.152884d0*ln_reynolds
     &        -0.0097785d0*ln_reynolds**2d0+0.00020903d0
     &        *ln_reynolds**3d0
!
         cdi_noz=cdi_noz_lr
!
      elseif (reynolds.lt.50000d0) then
!     
         if (reynolds.eq.0) then
            reynolds_cor=1
         else
            reynolds_cor=reynolds
         endif
!
         e=2.718281828459045d0
         ln_reynolds=log(reynolds_cor)/log(e)
!     
         cdi_noz_lr=0.19436d0+0.152884d0*ln_reynolds
     &        -0.0097785d0*ln_reynolds**2+0.00020903d0
     &        *ln_reynolds**3d0
!     
         cdi_noz_hr=0.9975d0-0.00653d0*dsqrt(1000000d0/50000d0)
         
!     linear interpolation in order to achieve continuity
!     
         cdi_noz=cdi_noz_lr+(cdi_noz_hr-cdi_noz_lr)
     &        *(reynolds-40000d0)/(50000d0-40000d0)
      else
!
!     formerly pk_cdi_noz_hr for high Reynolds numbers
!
         cdi_noz=0.9975d0-0.00653d0*dsqrt(1000000d0/reynolds)
      endif
      
      return
      end
