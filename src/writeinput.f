!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine writeinput(inpc,ipoinp,inp,nline,ninp,ipoinpc)
!
      implicit none
!
      integer nentries
      parameter(nentries=17)
!
      character*1 inpc(*)
      character*20 nameref(nentries)
!
      integer nline,i,j,ninp,ipoinp(2,nentries),inp(3,ninp),
     &  ipoinpc(0:*)
!
      data nameref /'RESTART,READ','NODE','USERELEMENT','ELEMENT',
     &              'NSET',
     &              'ELSET','SURFACE','TRANSFORM','MATERIAL',
     &              'ORIENTATION','TIE','INTERACTION',
     &              'INITIALCONDITIONS','AMPLITUDE',
     &              'CONTACTPAIR','COUPLING','REST'/
!
      open(16,file='input.inpc',status='unknown',err=161)
      do i=1,nline
         write(16,'(1x,i6,1x,1320a1)') i,
     &       (inpc(j),j=ipoinpc(i-1)+1,ipoinpc(i))
      enddo
      close(16)
!
      open(16,file='input.ipoinp',status='unknown',err=162)
      do i=1,nentries
         write(16,'(1x,a20,1x,i6,1x,i6)') nameref(i),(ipoinp(j,i),j=1,2)
      enddo
      close(16)
!
      open(16,file='input.inp',status='unknown',err=163)
      do i=1,ninp
         write(16,'(1x,i3,1x,i6,1x,i6,1x,i6)') i,(inp(j,i),j=1,3)
      enddo
      close(16)
!
      return
!
 161  write(*,*) '*ERROR in writeinput: could not open file input.inpc'
      call exit(201)
!
 162  write(*,*) 
     &    '*ERROR in writeinput: could not open file input.ipoinp'
      call exit(201)
!
 163  write(*,*) '*ERROR in writeinput: could not open file input.inp'
      call exit(201)
      end
