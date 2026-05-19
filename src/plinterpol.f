!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine plinterpol(plcon,nplcon,itemp,f,df,npmat_,ntmat_,
     &  imat,nelem,epl)
!
      implicit none
!
!     interpolation of isotropic or kinematic hardening data
!     input: hardening data plcon and nplcon, temperature itemp,
!       size parameters npmat_ and ntmat_, material number imat
!       and equivalent plastic strain at which the coefficients
!       are to be determined
!     output: hardening coefficient and its local derivative f and df
!
      integer npmat_,ntmat_,nplcon(0:ntmat_,*),itemp,ndata,imat,j,
     &  nelem
!
      real*8 plcon(0:2*npmat_,ntmat_,*),f,df,epl
!
      intent(in) plcon,nplcon,itemp,npmat_,ntmat_,
     &  imat,nelem,epl
!
      intent(inout) f,df
!
      ndata=nplcon(itemp,imat)
!
      do j=1,ndata
         if(epl.lt.plcon(2*j,itemp,imat)) exit
      enddo
!
      if((j.eq.1).or.(j.gt.ndata)) then
         if(j.eq.1) then
            f=plcon(1,itemp,imat)
            df=0.d0
         else
            f=plcon(2*ndata-1,itemp,imat)
            df=0.d0
         endif
         write(*,*) '*WARNING in plinterpol: plastic strain ',epl
         write(*,*) '         outside material plastic strain range'
         write(*,*) '         in element ',nelem,' and material ',imat
         write(*,*) '         for temperature ',plcon(0,itemp,imat)
      else
         df=(plcon(2*j-1,itemp,imat)-plcon(2*j-3,itemp,imat))/
     &      (plcon(2*j,itemp,imat)-plcon(2*j-2,itemp,imat))
         f=plcon(2*j-3,itemp,imat)+
     &      df*(epl-plcon(2*j-2,itemp,imat))
      endif
!
      return
      end
