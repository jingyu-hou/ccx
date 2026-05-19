!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
     
      subroutine friction_coefficient(l,d,ks,reynolds,form_fact,lambda)
!     
      implicit none
!     
      real*8 l,d,ks,reynolds,form_fact,lambda,alfa2,
     &     rey_turb_min,rey_lam_max,lzd,dd,ds,friction,dfriction,
     &     lambda_kr,lambda_turb,ksd
!     
      intent(in) l,d,ks,reynolds,form_fact
!
      intent(inout) lambda
!
      rey_turb_min=4000
      rey_lam_max=2000
      lzd=l/d
      ksd=ks/d
!     
!     transition laminar turbulent domain
!     
      if((reynolds.gt.rey_lam_max).and.(reynolds.lt.rey_turb_min))then
!     
         lambda_kr=64.d0/rey_lam_max
!     
!     Solving the implicit White-Colebrook equation
!     1/dsqrt(friction)=-2*log10(2.51/(Reynolds*dsqrt(friction)+0.27*Ks))
!     
!     Using Haaland explicit relationship for the initial friction value
!     S.E. Haaland 1983 (Source en.Wikipwedia.org)
!     
         friction=(-1.8d0*dlog10(6.9d0/4000.d0+(ksd/3.7d0)**1.11d0))
     &        **(-2d0)
!     
         do
            ds=dsqrt(friction)
            dd=2.51d0/(4000.d0*ds)+0.27d0*ksd
            dfriction=(1.d0/ds+2.d0*dlog10(dd))*2.d0*friction*ds/
     &           (1.d0+2.51d0/(4000.d0*dd))
            if(dfriction.le.friction*1.d-3) then
               friction=friction+dfriction
               exit
            endif
            friction=friction+dfriction
         enddo
         lambda_turb=friction
         
!     
!     logarithmic interpolation in the trans laminar turbulent domain
!     
         lambda=lambda_kr*(lambda_turb/lambda_kr)
     &       **(log(reynolds/rey_lam_max)/log(rey_turb_min/rey_lam_max))
!     
!     laminar flow
!     using Couette-Poiseuille formula
!     the form factor for non round section can be found in works such as
!     Bohl,W
!     "Technische Strömungslehre Stoffeigenschaften von Flüssigkeiten und
!     Gasen, hydrostatik,aerostatik,incompressible Strömungen,
!     Strömungsmesstechnik
!     Vogel Würzburg Verlag 1980
!     
      elseif(reynolds.lt.rey_lam_max) then
         lambda=64.d0/reynolds
         lambda=form_fact*lambda
!     
!     turbulent
!     
      else
!     Solving the implicit White-Colebrook equation
!     1/dsqrt(friction)=-2*log10(2.51/(Reynolds*dsqrt(friction)+0.27*Ks))
!     
!     Using Haaland explicit relationship for the initial friction value
!     S.E. Haaland 1983 (Source en.Wikipwedia.org)
!     
         friction=(-1.8d0*dlog10(6.9d0/reynolds+(ksd/3.7d0)
     &        **1.11d0))**-2d0
!     
         do
            ds=dsqrt(friction)
            dd=2.51d0/(reynolds*ds)+0.27d0*ksd
            dfriction=(1.d0/ds+2.d0*dlog10(dd))*2.d0*friction*ds/
     &           (1.d0+2.51d0/(reynolds*dd))
            if(dfriction.le.friction*1.d-3) then
               friction=friction+dfriction
               exit
            endif
            friction=friction+dfriction
         enddo
         lambda=friction
      endif
!     
      call interpol_alfa2(lzd,reynolds,alfa2)
!     
      return
!     
      end
      
      
      
