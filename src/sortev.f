!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
    
!     ADDITIONAL INPUT Parameters:
!     ipev:     Position of Eigenvalues, saves original Position of 
!            Eigenvalues before sorting
!     eigxr:   Real Part of Eigenvalues out of eigxx, used for
!            sorting Eigenvalues in increasing order
!
      subroutine sortev(nev,nmd,eigxx,cyclicsymmetry,x,eigxr,ipev,
     &     istartnmd,iendnmd,a,b)
!     
!     sorts the eigenvalues and eigenvectors of complex frequency
!     
      implicit none
!     
      integer nev,i,j,k,l,m,istartnmd(*),iendnmd(*),nmd,
     &     cyclicsymmetry,ipev(*)
!     
      real*8 eigxr(*)
!     
      complex*16 eigxx(*),a(nev),b(nev,*),x(nev,*)
!     
      if (cyclicsymmetry.eq.0)then
!     
!     sorting the eigenvalues according to their size
!     
         do i=1,nev
            ipev(i)=i
            eigxr(i)=cdabs(eigxx(i))
         enddo
         call dsort(eigxr,ipev,nev,2)
!     
!     sorting the eigenvectors
!     
         do i=1,nev
            a(i)=eigxx(ipev(i))
            do j=1,nev
               b(j,i)=x(j,ipev(i))
            enddo
         enddo
!     
!     copying in the original fields
!     
         do i=1,nev
            eigxx(i)=a(i)
            do j=1,nev
               x(i,j)=b(i,j)
               enddo
            enddo
         else
!     
!     Cyclic Symmetry
!     
            do l=1,nmd
!     
!     sorting the eigenvalues according to their size
!     
!     
            do i=istartnmd(l),iendnmd(l)
               if (l.eq.1) then
                  ipev(i)=i
                  eigxr(i)=cdabs(eigxx(i))
                  k=i
               else
                  k=i-istartnmd(l)+1
                  ipev(k)=i
                  eigxr(i)=cdabs(eigxx(i))
               endif
            enddo
            call dsort (eigxr,ipev,k,2)
!     
!     sorting the eigenvectors
!     
            do i=istartnmd(l),iendnmd(l)
               if (l.eq.1) then
               m=ipev(i)
               a(i)=eigxx(m)
               do j=istartnmd(l),iendnmd(l)
                  b(j,i)=x(j,m)
               enddo
            else
               k=i-istartnmd(l)+1
               a(i)=eigxx(ipev(k))
               do j=istartnmd(l),iendnmd(l)
                  b(j,i)=x(j,m)
               enddo
            endif
            enddo
         enddo
!     
!     copying in the original fields
!     
         do l=1,nmd
            if((a(istartnmd(l)).ne.0).and.
     &           (b(istartnmd(l),istartnmd(l)).ne.0))then
               do i=istartnmd(l),iendnmd(l)
               eigxx(i)=a(i)
               do j=istartnmd(l),iendnmd(l)
                  x(i,j)=b(i,j)
               enddo
            enddo
         endif
      enddo   
      endif
!     
      return
      end
      
