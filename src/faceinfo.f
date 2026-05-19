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
!     subroutine to find the right # nodes for element and surface  
!     based on current element number nelem and face number jface
!
      subroutine faceinfo(nelem,jface,lakon,nope,nopes,mint2d)
!
!     autor: Saskia Sitzmann 
!
      implicit none
!     
      character*8 lakon(*)
!     
      integer  nopes,nope,nelem,jface,mint2d
!     
      if(lakon(nelem)(4:5).eq.'8R') then
         mint2d=1
         nopes=4
         nope=8
      elseif(lakon(nelem)(4:4).eq.'8') then
         mint2d=4
         nopes=4
         nope=8
      elseif(lakon(nelem)(4:6).eq.'20R') then
         mint2d=4
         nopes=8
         nope=20
      elseif(lakon(nelem)(4:4).eq.'2') then
         mint2d=9
         nopes=8
         nope=20
      elseif(lakon(nelem)(4:5).eq.'10') then
         mint2d=3
         nopes=6
         nope=10
      elseif(lakon(nelem)(4:4).eq.'4') then
         mint2d=1
         nopes=3
         nope=4
      endif
!     
!     treatment of wedge faces
!     
      if(lakon(nelem)(4:4).eq.'6') then
         mint2d=1
         nope=6
         if(jface.le.2) then
            nopes=3
         else
            nopes=4
         endif
      endif
      if(lakon(nelem)(4:5).eq.'15') then
         nope=15
         if(jface.le.2) then
            mint2d=3
            nopes=6
         else
            mint2d=4
            nopes=8
         endif
      endif
!     
      return
      end
      
