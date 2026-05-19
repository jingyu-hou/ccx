!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine createelemneigh(nk,iponoel,inoel,istartnneigh,
     &   ialnneigh,icheckelems,istarteneigh,ialeneigh)
!
      implicit none
!
      integer nk,iponoel(*),inoel(2,*),istartnneigh(*),ialnneigh(*),
     &   ifree,index,i,j,ipos,na,nb,node,istarteneigh(*),ialeneigh(*),
     &   icheckelems(*)
!
!     determining all the elems of the neighbouring 
!     nodes of node nk.
!     They are stored in ialeneigh(istarteneigh(i))..
!     ...up to..... ialeneigh(istarteneigh(i+1)-1)
!
      ifree=1
      do i=1,nk
!
         istarteneigh(i)=ifree 
         index=iponoel(i)
         if(index.eq.0) cycle 
         na=istartnneigh(i)
         nb=istartnneigh(i+1)-1
!   
         do j=na,nb
!   
            node=ialnneigh(j)
            index=iponoel(node)
!      
            do
               if(index.eq.0) exit
               ipos=inoel(1,index)
               if(icheckelems(ipos).ne.i) then
                  ialeneigh(ifree)=inoel(1,index)
                  ifree=ifree+1
                  icheckelems(ipos)=i
               endif
               index=inoel(2,index)
            enddo
         enddo
      enddo   
      istarteneigh(nk+1)=ifree
!
      return
      end
