!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine writeobj(objectset,iobject,g0)
!
!     writes the results of the objective function in the .dat file
!
      implicit none
!
      character*81 objectset(4,*)
      integer iobject,i
      real*8 g0(*)
!          
      i=iobject+1
!
      if(objectset(1,i)(1:12).eq.'DISPLACEMENT') then
         write(5,*)
         write(5,*)'OBJECTIVE: DISPLACEMENT'  
         write(5,*)
         write(5,'(7x,e14.7)') g0(i)
!
      elseif(objectset(1,i)(1:14).eq.'EIGENFREQUENCY') then
         write(5,*)
         write(5,*)'OBJECTIVE: EIGENFREQUENCY'  
         write(5,*)
         write(5,'(7x,e14.7)') g0(i)
!
      elseif(objectset(1,i)(1:4).eq.'MASS') then
         write(5,*)
         write(5,*)'OBJECTIVE: MASS'  
         write(5,*)
         write(5,'(7x,e14.7)') g0(i)
!
      elseif(objectset(1,i)(1:11).eq.'SHAPEENERGY') then
         write(5,*)
         write(5,*)'OBJECTIVE: SHAPE ENERGY' 
         write(5,*) 
         write(5,'(7x,e14.7)') g0(i)
!
      elseif(objectset(1,i)(1:6).eq.'STRESS') then
         write(5,*)
         write(5,*)'OBJECTIVE: STRESS'  
         write(5,*)
         write(5,'(7x,e14.7)') g0(i)
!
      elseif(objectset(1,i)(1:9).eq.'THICKNESS') then
         write(5,*)
         write(5,*)'OBJECTIVE: THICKNESS'  
         write(5,*)
         write(5,'(7x,e14.7)') g0(i)
!
      elseif(objectset(1,i)(1:9).eq.'FIXGROWTH') then
         write(5,*)
         write(5,*)'OBJECTIVE: FIX GROWTH'  
         write(5,*)
         write(5,'(7x,e14.7)') g0(i)
!
      elseif(objectset(1,i)(1:12).eq.'FIXSHRINKAGE') then
         write(5,*)
         write(5,*)'OBJECTIVE: FIX SHRINKAGE'  
         write(5,*)
         write(5,'(7x,e14.7)') g0(i)
      endif
!      
      return
      end

