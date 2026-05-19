!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine thickness(dgdx,nobject,nodedesiboun,ndesiboun,
     &           objectset,xo,yo,zo,x,y,z,nx,ny,nz,co,ifree,ndesia,
     &           ndesib,iobject,ndesi,dgdxglob,nk)                       
!
!     calcualtion of the actual wall thickness      
!
      implicit none
!
      character*81 objectset(4,*)

      integer nobject,nodedesiboun(*),nnodesinside,i,ndesi,
     &        ndesiboun,j,k,neighbor(10),nx(*),ny(*),nz(*),
     &        istat,ndesia,ndesib,ifree,nnodes,irefnode,
     &        iactnode,iobject,nk
!
      real*8 dgdx(ndesi,nobject),xo(*),yo(*),zo(*),x(*),
     &       y(*),z(*),refdist,co(3,*),xdesi,ydesi,zdesi,      
     &       scale,actdist,dd,xrefnode,yrefnode,zrefnode,
     &       numericdist,dgdxglob(2,nk,nobject)
!   
      read(objectset(1,iobject)(61:80),'(f20.0)',iostat=istat) refdist
!      
      nnodes=1
      numericdist=0.02d0
!
!     check if upper or lower constraint is defined
!      
      if(objectset(1,iobject)(19:20).eq.'LE') then
         refdist=refdist*(1-numericdist)
      else
         refdist=refdist*(1+numericdist)
      endif
!
!     find minimum distance 
!
      do j=ndesia,ndesib
!         
         iactnode=nodedesiboun(j)    
         xdesi=co(1,iactnode)
         ydesi=co(2,iactnode)
         zdesi=co(3,iactnode)
!
         call near3d(xo,yo,zo,x,y,z,nx,ny,nz,xdesi,ydesi,zdesi,
     &        ifree,neighbor,nnodes)
!  
!        Calculate distance and check if node is active
!
         irefnode=neighbor(1)
         xrefnode=xo(irefnode)
         yrefnode=yo(irefnode)
         zrefnode=zo(irefnode)
         dd=(xrefnode-xdesi)**2+(yrefnode-ydesi)**2
     &      +(zrefnode-zdesi)**2
         actdist=dsqrt(dd)
!
         if(objectset(1,iobject)(19:20).eq.'LE') then
            if(actdist.gt.refdist) then
               dgdxglob(1,iactnode,iobject)=1.0d0
	       dgdxglob(2,iactnode,iobject)=1.0d0
c               dgdx(j,iobject)=1.d0
            endif
         else
            if(actdist.lt.refdist) then
               dgdxglob(1,iactnode,iobject)=1.0d0
	       dgdxglob(2,iactnode,iobject)=1.0d0
c               dgdx(j,iobject)=1.d0
            endif
         endif      
      enddo
!
      return        
      end
