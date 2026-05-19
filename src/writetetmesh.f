!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine writetetmesh(kontet,netet,cotet,nktet,field,nfield)
!
      implicit none
!
      character*1 c
      character*3 m1,m2,m3,m4,m5
      character*5 p0,p1,p2,p3,p4,p9999
      character*132 text
!
      integer kontet(4,*),netet,one,i,j,nsom,
     &  nktet,index,node,kode,nfield
!
      real*8 cotet(3,*),field(*),time
!
      c='C'
!
      m1=' -1'
      m2=' -2'
      m3=' -3'
      m4=' -4'
      m5=' -5'
!
      p0='    0'
      p1='    1'
      p2='    2'
      p3='    3'
      p4='    4'
      p9999=' 9999'
!
      one=1
      kode=1
      time=0.d0
!
      open(9,file='TetMasterSubmodel.frd',status='unknown')
!
      write(9,'(a5,a1)') p1,c
!
!       storing the coordinates of the nodes
!
      write(9,'(a5,a1,67x,i1)') p2,c,one
!
      do i=1,nktet
         write(9,100) m1,i,(cotet(j,i),j=1,3)
      enddo
!
      write(9,'(a3)') m3
!
!       storing the element topology
!
      write(9,'(a5,a1,67x,i1)') p3,c,one
!
      do i=1,netet
         if(kontet(1,i).eq.0) cycle
         write(9,'(a3,i10,3a5)') m1,i,p3,p0,p0
         write(9,'(a3,10i10)') m2,(kontet(j,i),j=1,4)
      enddo
      write(*,*) 'number of tetrahedra = ',netet
!
      write(9,'(a3)') m3
      write(9,'(a5)') p9999
!
      close(9)
!
 100  format(a3,i10,1p,3e12.5)
!
      return
      end
