!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine forcesolve(zc,nev,a,b,x,eiga,eigb,eigxx,iter,d,
     &  neq,z,istartnmd,iendnmd,nmd,cyclicsymmetry,neqact,
     &  igeneralizedforce)
!
!     solves for the complex eigenfrequencies due to Coriolis 
!     forces
!
      implicit none
!
      logical wantx
!
      integer nev,neq,iter(*),i,j,k,l,istartnmd(*),iendnmd(*),nmd,
     &  cyclicsymmetry,neqact,igeneralizedforce
!
      real*8 z(neq,*),d(*)
!
      complex*16 a(nev,*),b(nev,*),x(nev,*),eiga(*),eigb(*),eigxx(*),
     &  zc(neqact,*)
!
      if(igeneralizedforce.eq.0) then
!
!        no generalized force: multiplication with the eigenmodes
!        is necessary
!
         if(cyclicsymmetry==0) then
            do i=1,nev
               do j=1,nev
                  do k=1,neq
                     a(i,j)=a(i,j)+z(k,i)*zc(k,j)
                  enddo
               enddo
               write(*,*) 
     &              'aerodynamic stiffness/structural stiffness = ',
     &              a(i,i)/d(i)
               a(i,i)=a(i,i)+d(i)
               b(i,i)=(1.d0,0.d0)
            enddo
         else
!     
!     cyclic symmetry
!     
            do l=1,nmd
               do i=istartnmd(l),iendnmd(l)
                  do j=istartnmd(l),iendnmd(l)
                     do k=1,neqact
                        a(i,j)=a(i,j)+z(k,i)*zc(k,j)-
     &                       z(k+neqact,i)*zc(k,j)*(0.d0,1.d0)
                     enddo
                  enddo
                  write(*,*) 
     &                 'aerodynamic stiffness/structural stiffness = ',
     &                 a(i,i)/d(i)
                  a(i,i)=a(i,i)+d(i)
                  b(i,i)=(1.d0,0.d0)
               enddo
            enddo
         endif
      else
!
!        generalized force: the a-matrix is (apart from the diagonal)
!        known
!
         if(cyclicsymmetry==0) then
            do i=1,nev
               write(*,*) 
     &              'aerodynamic stiffness/structural stiffness = ',
     &              a(i,i)/d(i)
               a(i,i)=a(i,i)+d(i)
               b(i,i)=(1.d0,0.d0)
            enddo
         else
!     
!     cyclic symmetry
!     
            do l=1,nmd
               do i=istartnmd(l),iendnmd(l)
                  write(*,*) 
     &                 'aerodynamic stiffness/structural stiffness = ',
     &                 a(i,i)/d(i)
                  a(i,i)=a(i,i)+d(i)
                  b(i,i)=(1.d0,0.d0)
               enddo
            enddo
         endif
      endif
!     
      wantx=.true.
!
!     solving for the complex eigenvalues
!     
      call dlzhes(nev,a,nev,b,nev,x,nev,wantx)
      call dlzit(nev,a,nev,b,nev,x,nev,wantx,iter,eiga,eigb)
!     
      do i=1,nev
         if(iter(i).eq.-1) then
            write(*,*) '*ERROR in coriolissolve: fatal error'
            write(*,*) '       in dlzit'
            call exit(201)
         elseif(cdabs(eigb(i)).lt.1.d-10) then
            write(*,*) '*ERROR in coriolissolve: eigenvalue'
            write(*,*) '       out of bounds'
            call exit(201)
         else
            eigxx(i)=cdsqrt(eiga(i)/eigb(i))
         endif
      enddo
!     
      return
      end
      
      
