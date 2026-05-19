!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine postinitialnet(ieg,lakon,v,ipkon,kon,nflow,prop,
     &     ielprop,ielmat,ntmat_,shcon,nshcon,rhcon,nrhcon,mi,iponoel,
     &     inoel,itg,ntg)
!
!     this routine only applies to compressible networks
!
!     determination of initial values based on the boundary conditions
!     and the initial values given by the user by propagating these
!     through the network (using information on the mass flow direction
!     derived from unidirectional network elements or mass flow given
!     by the user (boundary conditions or initial conditions)) 
!
!     it is assumed that mass flows cannot
!     be identically zero (a zero mass flow leads to convergence problems).
!
!     This routine is used for elements in which the pressure gradient
!     does not allow to determine the mass flow, e.g. the free vortex
!     and the forced vortex
!     
      implicit none
!
      character*8 lakon(*)
!           
      integer mi(*),ieg(*),nflow,i,ielmat(mi(3),*),ntmat_,node1,node2,
     &     nelem,index,nshcon(*),ipkon(*),kon(*),nodem,imat,ielprop(*),
     &     nrhcon(*),neighbor,ichange,iponoel(*),inoel(2,*),indexe,
     &     itg(*),ntg,j
!     
      real*8 prop(*),shcon(0:3,ntmat_,*),xflow,v(0:mi(2),*),cp,r,
     &     dvi,rho,rhcon(0:1,ntmat_,*),kappa,cti,Ti,ri,ro,p1zp2,omega,
     &     p2zp1
!
c      write(*,*) 'postinitialnet '
c      do i=1,ntg
c         write(*,'(i10,3(1x,e11.4))') itg(i),(v(j,itg(i)),j=0,2)
c      enddo
!
      do
         ichange=0
!
!        propagation of the mass flow through the network
!
         do i=1,nflow
            nelem=ieg(i)
            indexe=ipkon(nelem)
            nodem=kon(indexe+2)
!
            if(dabs(v(1,nodem)).le.0.d0) then
!
!              no initial mass flow given yet
!              check neighbors for mass flow (only if not
!              branch nor joint)
!
!              first end node
!
               node1=kon(indexe+1)
!
               if(node1.ne.0) then
                  index=iponoel(node1)
!
                  if(inoel(2,inoel(2,index)).eq.0) then
!
!                 no branch nor joint; determine neighboring element
!
                     if(inoel(1,index).eq.nelem) then
                        neighbor=inoel(1,inoel(2,index))
                     else
                        neighbor=inoel(1,index)
                     endif
!
!                 initial mass flow in neighboring element
!
                     xflow=v(1,kon(ipkon(neighbor)+2))
!
                     if(dabs(v(1,nodem)).gt.0.d0) then
!
!                    propagate initial mass flow
!
                        if(dabs(xflow).gt.0.d0) then
                           v(1,nodem)=xflow
                           ichange=1
                           cycle
                        endif
                     else
!
!                    propagate only the sign of the mass flow
!
                        if(dabs(xflow).gt.0.d0) then
                           v(1,nodem)=xflow
                           ichange=1
                           cycle
                        endif
                     endif
                  endif
               endif
!
!              second end node
!
               node2=kon(indexe+3)
!
               if(node2.ne.0) then
                  index=iponoel(node2)
!
                  if(inoel(2,inoel(2,index)).eq.0) then
!
!                 no branch nor joint; determine neighboring element
!
                     if(inoel(1,index).eq.nelem) then
                        neighbor=inoel(1,inoel(2,index))
                     else
                        neighbor=inoel(1,index)
                     endif
!
!                 initial mass flow in neighboring element
!
                     xflow=v(1,kon(ipkon(neighbor)+2))
!
                     if(dabs(v(1,nodem)).gt.0.d0) then
!
!                    propagate initial mass flow
!
                        if(dabs(xflow).gt.0.d0) then
                           v(1,nodem)=xflow
                           ichange=1
                           cycle
                        endif
                     else
!
!                    propagate only the sign of the mass flow
!
                        if(dabs(xflow).gt.0.d0) then
                           v(1,nodem)=xflow
                           ichange=1
                           cycle
                        endif
                     endif
                  endif
               endif
            endif
         enddo
c         write(*,*) 'postinitialnet '
c         do i=1,ntg
c            write(*,'(i10,3(1x,e11.4))') itg(i),(v(j,itg(i)),j=0,2)
c         enddo
         if(ichange.eq.0) exit
      enddo
!     
      return
      end
      
      
