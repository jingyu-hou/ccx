!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine propertynet(ieg,nflow,prop,ielprop,lakon,iin,
     &       prop_store,ttime,time,nam,amname,namta,amta)
!
!     user subroutine propertynet
!
!
!     INPUT:
!
!     ieg(i)             global element number corresponding to
!                        network element i (i=1,...,nflow)
!     nflow              number of network elements
!     ielprop(i)         property to the position in fields prop and
!                        prop_store after which the properties for 
!                        element i start (prop(ielprop(i)+1),
!                        prop(ielprop(i)+2).....). The number is dictated
!                        by the type of element.
!     lakon(i)           label of element i
!     iin                gas network iteration number
!     prop_store         property values as specified in the
!                        input deck
!     ttime              total time
!     time               step time
!     nam                number of amplitudes
!     amname(i)          amplitude name of amplitude i
!     namta(1,i)         location of first (time,amplitude) pair in
!                        field amta
!     namta(2,i)         location of last (time,amplitude) pair in
!                        field amta
!     namta(3,i)         in absolute value the amplitude it refers to; if
!                        abs(namta(3,i))=i it refers to itself. If
!                        abs(namta(3,i))=j, amplitude i is a time delay
!                        of amplitude j; in that case the value of the
!                        time delay is stored in
!                        amta(1,namta(1,i)); in the latter case
!                        amta(2,namta(1,i)) is without meaning; if
!                        namta(3,i)>0 the time in amta for amplitude i is
!                        step time, else it is total time.
!     amta(1,i)          time of (time,amplitude)-pair i
!     amta(2,i)          amplitude of (time,amplitude)-pair i
!
!     OUTPUT:
!
!     prop               actual property values
!           
      implicit none
!
      character*8 lakon(*)
      character*80 amname(*)
!
      integer ieg(*),nflow,ielprop(*),iin,nam,namta(3,*)
!
      real*8 prop(*),prop_store(*),ttime,time,amta(2,*)
!
      intent(in) ieg,nflow,ielprop,lakon,iin,
     &           prop_store,ttime,time,nam,amname,namta,amta
!
      intent(out) prop
!
      return
      end

