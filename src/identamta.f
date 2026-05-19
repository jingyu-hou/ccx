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
!     identifies the position id of reftime in an ordered array
!     amta(1,istart...iend) of real numbers; amta is defined as amta(2,*)
!  
!     id is such that amta(1,id).le.reftime and amta(1,id+1).gt.reftime
!                                                                             
      subroutine identamta(amta,reftime,istart,iend,id)
!
      implicit none
!
      integer id,istart,iend,n2,m
      real*8 amta(2,*),reftime
      id=istart-1
      if(iend.lt.istart) return
      n2=iend+1
      do                                                         
         m=(n2+id)/2
         if(reftime.ge.amta(1,m)) then
            id=m     
         else
            n2=m  
         endif
         if((n2-id).eq.1) return
      enddo
      end
