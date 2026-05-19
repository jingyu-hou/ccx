!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine fill_neiel(nef,ipnei,neiel,neielcp)
!
!     copy neiel into neielcp, thereby substituting the zero's by
!     neighboring values
!
      implicit none
!
      integer nef,ipnei(*),neiel(*),neielcp(*) ,i,j,indexf
!
      intent(in) nef,ipnei,neiel
!
      intent(inout) neielcp
!
      do i=1,nef
         do indexf=ipnei(i)+1,ipnei(i+1)
            if(neiel(indexf).eq.0) then
               if(indexf.eq.ipnei(i)+1) then
                  do j=ipnei(i)+2,ipnei(i+1)
                     if(neiel(j).ne.0) then
                        neielcp(indexf)=neiel(j)
                        exit
                     endif
                  enddo
               else
                  neielcp(indexf)=neielcp(indexf-1)
               endif
            else
               neielcp(indexf)=neiel(indexf)
            endif
         enddo
      enddo
!     
      return
      end
