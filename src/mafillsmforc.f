!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine mafillsmforc(nforc,ndirforc,nodeforc,xforc,nactdof,
     &  fext,ipompc,nodempc,coefmpc,mi,rhsi,fnext,
     &  nmethod,ntrans,inotr,trab,co)
!
!     including point forces into the external force vector
!
      implicit none
!
      integer nforc,ndirforc(*),nodeforc(2,*),mi(*),nactdof(0:mi(2),*),
     &  ipompc(*),nodempc(3,*),i,jdof,id,ist,
     &  index,nmethod,rhsi,ntrans,inotr(2,*),node,idir,j,itr
!
      real*8 xforc(*),fext(*),coefmpc(*),fnext(0:mi(2),*),trab(7,*),
     &  a(3,3),co(3,*)
!
      if(rhsi.eq.1) then
!
!        point forces
!      
         do i=1,nforc
            if(ndirforc(i).gt.mi(2)) cycle
!
            node=nodeforc(1,i)
!
!           check for transformation
!
            if(ntrans.eq.0) then
               itr=0
            else
               itr=inotr(1,node)
            endif
!
            if(itr.eq.0) then
!
!              no transformation
!
!              updating the external force vector for dynamic
!              calculations
!
               if(nmethod.eq.4) fnext(ndirforc(i),nodeforc(1,i))=
     &              fnext(ndirforc(i),nodeforc(1,i))+xforc(i)
!     
               jdof=nactdof(ndirforc(i),nodeforc(1,i))
               if(jdof.gt.0) then
                  fext(jdof)=fext(jdof)+xforc(i)
               else
!     
!                 node is a dependent node of a MPC: distribute
!                 the forces among the independent nodes
!                 (proportional to their coefficients)
!     
                  if(jdof.ne.2*(jdof/2)) then
                     id=(-jdof+1)/2
                     ist=ipompc(id)
                     index=nodempc(3,ist)
                     if(index.eq.0) cycle
                     do
                        jdof=nactdof(nodempc(2,index),nodempc(1,index))
                        if(jdof.gt.0) then
                           fext(jdof)=fext(jdof)-
     &                          coefmpc(index)*xforc(i)/coefmpc(ist)
                        endif
                        index=nodempc(3,index)
                        if(index.eq.0) exit
                     enddo
                  endif
               endif
            else
!
!     transformation
!
               call transformatrix(trab(1,itr),co(1,node),a)
!     
               idir=ndirforc(i)
!     
!     updating the external force vector for dynamic
!     calculations
!     
               if(nmethod.eq.4) then
                  do j=1,3
                     fnext(j,node)=fnext(j,node)+a(j,idir)*xforc(i)
                  enddo
               endif
!     
               do j=1,3
                  jdof=nactdof(j,node)
                  if(jdof.gt.0) then
                     fext(jdof)=fext(jdof)+a(j,idir)*xforc(i)
                  else
!     
!     node is a dependent node of a MPC: distribute
!     the forces among the independent nodes
!     (proportional to their coefficients)
!     
                     if(jdof.ne.2*(jdof/2)) then
                        id=(-jdof+1)/2
                        ist=ipompc(id)
                        index=nodempc(3,ist)
                        if(index.eq.0) cycle
                        do
                           jdof=nactdof(nodempc(2,index),
     &                                  nodempc(1,index))
                           if(jdof.gt.0) then
                              fext(jdof)=fext(jdof)-
     &                             coefmpc(index)*a(j,idir)*xforc(i)
     &                             /coefmpc(ist)
                           endif
                           index=nodempc(3,index)
                           if(index.eq.0) exit
                        enddo
                     endif
                  endif
               enddo
            endif
         enddo
!     
      endif
!     
      return
      end
      
