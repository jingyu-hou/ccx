!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
!
!
      subroutine checkconstraint(nobject,objectset,g0,nactive,
     &   nnlconst,ipoacti,ndesi,dgdxglob,nk,nodedesi,iconstacti,
     &   objnorm,inameacti)               
!
!     check which constraints are active on the basis of the 
!     function values of the constraints     
!
      implicit none
!
      character*81 objectset(4,*)
      character*20 empty
!
      integer iobject,nobject,istat,nactive,nnlconst,inameacti(*),
     &   ipoacti(*),ifree,i,ndesi,nk,node,nodedesi(ndesi),
     &   iconstacti(*)
!
      real*8 g0(nobject),bounds(20),scale,bound,objnorm(nobject),
     &   dgdxglob(2,nk,nobject)
!
!     determine bounds of constraints
!
      empty='                    '
!
      do iobject=1,nobject
         if((objectset(1,iobject)(19:20).eq.'LE').or.
     &        (objectset(1,iobject)(19:20).eq.'GE')) then
            if((objectset(1,iobject)(1:9).eq.'THICKNESS').or.
     &         (objectset(1,iobject)(1:9).eq.'FIXGROWTH').or.
     &         (objectset(1,iobject)(1:12).eq.'FIXSHRINKAGE')) then
               do i=1,ndesi
                  node=nodedesi(i)
                  if(dgdxglob(2,node,iobject).gt.0) then
                     g0(iobject)=1.d0+g0(iobject)
                  endif
               enddo
            else 
               if(objectset(1,iobject)(61:80).eq.empty) then
                  read(objectset(1,iobject)(41:60),'(f20.0)',
     &                 iostat=istat) scale         
                  bounds(iobject)=g0(iobject)*scale
               else
                  read(objectset(1,iobject)(41:60),'(f20.0)',
     &                 iostat=istat) scale         
                  read(objectset(1,iobject)(61:80),'(f20.0)',
     &                 iostat=istat) bound
                  bounds(iobject)=bound*scale
               endif
            endif
         endif
      enddo
!
!     determine active constraints
!
      nactive=0
      nnlconst=0
      ifree=1
!
!     determine all nonlinear constraints
!
      do iobject=1,nobject
         if(objectset(1,iobject)(19:20).eq.'LE') then
            if(objectset(1,iobject)(1:9).eq.'THICKNESS') cycle 
            if(objectset(1,iobject)(1:9).eq.'FIXGROWTH') cycle 
            if(objectset(1,iobject)(1:12).eq.'FIXSHRINKAGE') cycle 
            objnorm(ifree)=g0(iobject)/bounds(iobject)-1
            if(objnorm(ifree).gt.-0.02) then
               nactive=nactive+1
               nnlconst=nnlconst+1
               ipoacti(ifree)=iobject
         inameacti(ifree)=iobject
         iconstacti(ifree)=-1
               ifree=ifree+1
            endif
         elseif(objectset(1,iobject)(19:20).eq.'GE') then
            if(objectset(1,iobject)(1:9).eq.'THICKNESS') cycle 
            if(objectset(1,iobject)(1:9).eq.'FIXGROWTH') cycle 
            if(objectset(1,iobject)(1:12).eq.'FIXSHRINKAGE') cycle 
            objnorm(ifree)=-1*(g0(iobject)/bounds(iobject))+1
            if(objnorm(ifree).gt.-0.02) then
               nactive=nactive+1
               nnlconst=nnlconst+1
               ipoacti(ifree)=iobject
               inameacti(ifree)=iobject
               iconstacti(ifree)=1
               ifree=ifree+1
            endif
         endif
      enddo
!
!     determine all linear constraints
!
      do iobject=1,nobject
         if(objectset(1,iobject)(19:20).eq.'LE') then
            if((objectset(1,iobject)(1:9).eq.'THICKNESS').or.
     &         (objectset(1,iobject)(1:9).eq.'FIXGROWTH')) then
               if(g0(iobject)>0) then
                  do i=1,ndesi
                     node=nodedesi(i)
                     if(dgdxglob(2,node,iobject).eq.1) then
                        ipoacti(ifree)=i
                        inameacti(ifree)=iobject
                        iconstacti(ifree)=-1
                        ifree=ifree+1
                        nactive=nactive+1
                     endif
                  enddo
               endif   
            endif
         elseif(objectset(1,iobject)(19:20).eq.'GE') then
            if((objectset(1,iobject)(1:9).eq.'THICKNESS').or.
     &         (objectset(1,iobject)(1:12).eq.'FIXSHRINKAGE')) then
               if(g0(iobject)>0) then
                  do i=1,ndesi
                     node=nodedesi(i)
                     if(dgdxglob(2,node,iobject).eq.1) then
                        ipoacti(ifree)=i
                        inameacti(ifree)=iobject
                        iconstacti(ifree)=1
                        ifree=ifree+1
                        nactive=nactive+1
                     endif
                  enddo
               endif   
            endif
         endif
      enddo
!
      return        
      end




