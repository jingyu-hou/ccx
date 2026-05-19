!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine writesubmatrix(submatrix,noderetain,ndirretain,
     &    nretain,jobnamec)
!
!     writing the matrix of a substructure to a .mtx-file
!
      implicit none
!
      character*132 jobnamec(*),fn
!
      integer i,j,nretain,ilen,noderetain(*),ndirretain(*)
!
      real*8 submatrix(nretain,nretain)
!
      ilen=index(jobnamec(5),' ')
      if(ilen.gt.129) then
         write(*,*) '*ERROR in writesubmatrix:'
         write(*,*) '       name of file for storing the submatrix'
         write(*,*) '       is too long (> 128 char); name = '
         write(*,*) jobnamec(5)(1:132)
         call exit(201)
      else
         fn(1:ilen-1)=jobnamec(5)(1:ilen-1)
         fn(ilen:ilen+3)='.mtx'
         do i=ilen+4,132
            fn(i:i)=' '
         enddo
      endif
!
      open(12,file=fn,status='unknown')
      write(12,100)
 100  format('**')
      write(12,101)
 101  format('** GENERATION OF SUBSTRUCTURE')
      write(12,102) nretain
 102  format('*USER ELEMENT,NODES= ',i10,',LINEAR')
      write(12,103)
 103  format('** ELEMENT NODES')
      write(12,104) (noderetain(i),i=1,nretain)
 104  format('**',i10,',',i10,',',i10,',',i10,',',i10,',',i10,',',
     &i10,',',i10,',',i10,',',i10,',')
      write(12,105) ndirretain(1)
 105  format(i10)
      write(12,106) (i,ndirretain(i),i=2,nretain)
 106  format(i10,',',i10)
      write(12,107)
 107  format('*MATRIX,TYPE=STIFFNESS')
!
      do j=1,nretain
         write(12,108) (submatrix(i,j),i=1,j)
      enddo
 108  format(e20.13,',',e20.13,',',e20.13,',',e20.13,',')
      close(12)
!
      return
      end

