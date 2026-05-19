!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine correctflux(nef,ipnei,neifa,neiel,flux,vfa,advfa,area,
     &  vel,alet,ielfa,ale,ifabou,nefa,nefb,xxnj,gradpcfa)
!
!     correction of v due to the balance of mass
!     the correction is in normal direction to the face
!
      implicit none
!
      integer i,nef,indexf,ipnei(*),neifa(*),neiel(*),ielfa(4,*),
     &  iel,ifa,ifabou(*),indexb,nefa,nefb
!
      real*8 flux(*),vfa(0:7,*),advfa(*),area(*),vel(nef,0:7),alet(*),
     &  ale(*),xxnj(3,*),gradpcfa(3,*)
!
!
!
      do i=nefa,nefb
c         totflux=0.d0
         do indexf=ipnei(i)+1,ipnei(i+1)
            ifa=neifa(indexf)
            iel=neiel(indexf)
            if(iel.gt.0) then
!
!              internal face
!
               flux(indexf)=flux(indexf)+vfa(5,ifa)*advfa(ifa)*
     &                      ((vel(i,4)-vel(iel,4))*alet(indexf)
     &                 -(gradpcfa(1,ifa)*xxnj(1,indexf)+
     &                   gradpcfa(2,ifa)*xxnj(2,indexf)+
     &                   gradpcfa(3,ifa)*xxnj(3,indexf)))
            else
               indexb=-ielfa(2,ifa)
               if(indexb.gt.0) then
                  if(((ifabou(indexb+1).eq.0).or.
     &                 (ifabou(indexb+2).eq.0).or.
     &                 (ifabou(indexb+3).eq.0)).and.
     &                 (ifabou(indexb+4).ne.0)) then
!
!                    external face with pressure boundary conditions
!
                     flux(indexf)=flux(indexf)
     &                      +vfa(5,ifa)*advfa(ifa)*(vel(i,4)*ale(indexf)
     &                             -(gradpcfa(1,ifa)*xxnj(1,indexf)+
     &                               gradpcfa(2,ifa)*xxnj(2,indexf)+
     &                               gradpcfa(3,ifa)*xxnj(3,indexf)))
                  endif
               endif
            endif
c            totflux=totflux+flux(indexf)
         enddo
      enddo
! 
      return
      end
