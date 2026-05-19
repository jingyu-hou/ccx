!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine calceth(alcon,nalcon,imat,iorien,pgauss,
     &  orab,ntmat_,alzero,t0l,t1l,eth)
!
      implicit none
!
!     determines the material data for element iel
!
!     istiff=0: only interpolation of material data
!     istiff=1: copy the consistent tangent matrix from the field
!               xstiff and check for zero entries
!
!
      integer nalcon(2,*),imat,iorien,j,k,kal(2,6),j1,j2,j3,j4,jj,
     &  id,two,seven,ntmat_
!
      real*8 alcon(0:6,ntmat_,*),eth(6),
     &  orab(7,*),alph(6),alzero(*),t0l,t1l,
     &  skl(3,3),xa(3,3),pgauss(3)
!
      intent(in) alcon,nalcon,imat,iorien,pgauss,orab,ntmat_,
     &  alzero,t0l,t1l
!
      intent(inout) eth
!
      kal=reshape((/1,1,2,2,3,3,1,2,1,3,2,3/),(/2,6/))
!
      two=2
      seven=7
!
!     calculating the expansion coefficients
!     
      call ident2(alcon(0,1,imat),t1l,nalcon(2,imat),seven,id)
      if(nalcon(2,imat).eq.0) then
         do k=1,6
            alph(k)=0.d0
         enddo
         continue
      elseif(nalcon(2,imat).eq.1) then
         do k=1,nalcon(1,imat)
            alph(k)=alcon(k,1,imat)*(t1l-alzero(imat))
         enddo
      elseif(id.eq.0) then
         do k=1,nalcon(1,imat)
            alph(k)=alcon(k,1,imat)*(t1l-alzero(imat))
         enddo
      elseif(id.eq.nalcon(2,imat)) then
         do k=1,nalcon(1,imat)
            alph(k)=alcon(k,id,imat)*(t1l-alzero(imat))
         enddo
      else
         do k=1,nalcon(1,imat)
            alph(k)=(alcon(k,id,imat)+
     &           (alcon(k,id+1,imat)-alcon(k,id,imat))*
     &           (t1l-alcon(0,id,imat))/
     &           (alcon(0,id+1,imat)-alcon(0,id,imat)))
     &           *(t1l-alzero(imat))
         enddo
      endif
!     
!     subtracting the initial temperature influence       
!     
      call ident2(alcon(0,1,imat),t0l,nalcon(2,imat),seven,id)
      if(nalcon(2,imat).eq.0) then
         continue
      elseif(nalcon(2,imat).eq.1) then
         do k=1,nalcon(1,imat)
            alph(k)=alph(k)-alcon(k,1,imat)*(t0l-alzero(imat))
         enddo
      elseif(id.eq.0) then
         do k=1,nalcon(1,imat)
            alph(k)=alph(k)-alcon(k,1,imat)*(t0l-alzero(imat))
         enddo
      elseif(id.eq.nalcon(2,imat)) then
         do k=1,nalcon(1,imat)
            alph(k)=alph(k)-alcon(k,id,imat)*(t0l-alzero(imat))
         enddo
      else
         do k=1,nalcon(1,imat)
            alph(k)=alph(k)-(alcon(k,id,imat)+
     &           (alcon(k,id+1,imat)-alcon(k,id,imat))*
     &           (t0l-alcon(0,id,imat))/
     &           (alcon(0,id+1,imat)-alcon(0,id,imat)))
     &           *(t0l-alzero(imat))
         enddo
      endif
!     
!     storing the thermal strains
!     
      if(nalcon(1,imat).eq.1) then
         do k=1,3
            eth(k)=alph(1)
         enddo
         do k=4,6
            eth(k)=0.d0
         enddo
      elseif(nalcon(1,imat).eq.3) then
         do k=1,3
            eth(k)=alph(k)
         enddo
         do k=4,6
            eth(k)=0.d0
         enddo
      else
         do k=1,6
            eth(k)=alph(k)
         enddo
      endif
!     
!
!     modifying the thermal constants if anisotropic and
!     a transformation was defined
!
      if((iorien.ne.0).and.(nalcon(1,imat).gt.1)) then
!     
!        calculating the transformation matrix
!     
         call transformatrix(orab(1,iorien),pgauss,skl)
!     
!        transforming the thermal strain
!     
         xa(1,1)=eth(1)
         xa(1,2)=eth(4)
         xa(1,3)=eth(5)
         xa(2,1)=eth(4)
         xa(2,2)=eth(2)
         xa(2,3)=eth(6)
         xa(3,1)=eth(5)
         xa(3,2)=eth(6)
         xa(3,3)=eth(3)
!     
         do jj=1,6
            eth(jj)=0.d0
            j1=kal(1,jj)
            j2=kal(2,jj)
            do j3=1,3
               do j4=1,3
                  eth(jj)=eth(jj)+
     &                 xa(j3,j4)*skl(j1,j3)*skl(j2,j4)
               enddo
            enddo
         enddo
      endif
!
      return
      end
