!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine transition(dgdxglob,nobject,nk,nodedesi,ndesi,
     &           objectset,xo,yo,zo,x,y,z,nx,ny,nz,co,ifree,
     &           ndesia,ndesib)                         
!
!     scaling of sensitivitites between the designspace
!     and non-designspace      
!
      implicit none
!
      character*81 objectset(4,*)

      integer nobject,nk,nodedesi(*),
     &        ndesi,j,m,neighbor(10),nx(*),ny(*),nz(*),
     &        istat,ndesia,ndesib,ifree,nnodes,irefnode,
     &        iactnode
!
      real*8 dgdxglob(2,nk,nobject),xo(*),yo(*),zo(*),x(*),
     &       y(*),z(*),trans,co(3,*),xdesi,ydesi,zdesi,      
     &       scale,actdist,dd,xrefnode,yrefnode,zrefnode
!
!     Read transition distance.
!   
      read(objectset(1,1)(21:40),'(f20.0)',iostat=istat) trans
!      
      nnodes=1
!
      do j=ndesia,ndesib
!         
         iactnode=nodedesi(j)
         xdesi=co(1,iactnode)
         ydesi=co(2,iactnode)
         zdesi=co(3,iactnode)
!
         call near3d(xo,yo,zo,x,y,z,nx,ny,nz,xdesi,ydesi,zdesi,
     &        ifree,neighbor,nnodes)
!  
!        Calculate scaled sensitivity 
!
         irefnode=neighbor(1)
         xrefnode=xo(irefnode)
         yrefnode=yo(irefnode)
         zrefnode=zo(irefnode)
         dd=(xrefnode-xdesi)**2+(yrefnode-ydesi)**2
     &      +(zrefnode-zdesi)**2
         actdist=dsqrt(dd)
         if(actdist.ge.trans) cycle
!
!        Linear scaling
         scale=actdist/trans
!        Exponential scaling
!        scale=1/(1+dexp(-actdist/trans*10+5))
         do m=1,nobject
            if(objectset(1,m)(1:9).eq.'THICKNESS') cycle
            if(objectset(1,m)(1:9).eq.'FIXGROWTH') cycle
            if(objectset(1,m)(1:12).eq.'FIXSHRINKAGE') cycle
            dgdxglob(2,iactnode,m)=dgdxglob(2,iactnode,m)*scale
         enddo
      enddo
!
      return        
      end
