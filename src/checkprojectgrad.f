!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine checkprojectgrad(nactiveold,nactive,ipoacti,ipoactiold,
     &   objectset,xlambd,nnlconst,iconstacti,iconstactiold,inameacti)         
!
!     checks the lagrange multipliers and reduces the number of active
!     constraints if possible
!
      implicit none
!
      character*81 objectset(4,*)
!
      integer nactiveold,nactive,ipoacti(*),ipoactiold(*),i,j,iconst,
     &   nnlconst,nnlconstold,iconstacti(*),iconstactiold(*),
     &   inameacti(*)
!
      real*8 xlambd(*)
!
!
      nactiveold=nactive
      nnlconstold=nnlconst
      do i=1,nactive
         ipoactiold(i)=ipoacti(i)
         iconstactiold(i)=iconstacti(i)
      enddo
!
      do i=1,nactiveold
         iconst=ipoactiold(i)
         if(i.le.nnlconstold) then           
            if(iconstactiold(i).eq.-1) then
               if(xlambd(i).lt.0.d0) then
                  nactive=nactive-1
                  nnlconst=nnlconst-1
                  do j=i+1,nactiveold
                     ipoacti(j-1)=ipoacti(j)
                     iconstacti(j-1)=iconstacti(j)
                     inameacti(j-1)=inameacti(j)
                  enddo
               endif
            elseif(iconstactiold(i).eq.1) then
               if(xlambd(i).gt.0.d0) then
                  nactive=nactive-1
                  nnlconst=nnlconst-1
                  do j=i+1,nactiveold
                     ipoacti(j-1)=ipoacti(j)
                     iconstacti(j-1)=iconstacti(j)
                     inameacti(j-1)=inameacti(j)
                  enddo
               endif
            endif
         else
            if(iconstactiold(i).eq.-1) then      
               if(xlambd(i).lt.0.d0) then
                  nactive=nactive-1
                  do j=i+1,nactiveold
                     ipoacti(j-1)=ipoacti(j)
                     iconstacti(j-1)=iconstacti(j)
                     inameacti(j-1)=inameacti(j)
                  enddo
               endif
            elseif(iconstactiold(i).eq.1) then
               if(xlambd(i).gt.0.d0) then
                  nactive=nactive-1
                  do j=i+1,nactiveold
                     ipoacti(j-1)=ipoacti(j)
                     iconstacti(j-1)=iconstacti(j)
                     inameacti(j-1)=inameacti(j)
                  enddo
               endif      
            endif
         endif
      enddo
!
      return        
      end
