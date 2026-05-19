!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine trianeighbor(ipe,ime,imastop,ncont,koncont,
     &  ifreeme)
!
!     Catalogueing the neighboring triangles for a given master
!     triangle
!
!     Authors: Li,Yang; Rakotonanahary, Samoela; 
!
      implicit none
!
      integer j,k,node,ipe(*),ime(4,*),imastop(3,*),ipos,node1,node2,
     &  index1,index1old,ifreeme,ncont,koncont(4,*)
!
!     catalogueing the edges in the triangulation
!     determining neighboring triangles
!
      ifreeme=0
      do j=1,ncont
         do k=1,3
            node1=koncont(k,j)
            if(k.eq.3) then
               node2=koncont(1,j)
            else
               node2=koncont(k+1,j)
            endif
!
            if(k.eq.1) then
               ipos=3
            else
               ipos=k-1
            endif
!
!           making sure that node1 < node2
!
            if(node1.gt.node2) then
               node=node1
               node1=node2
               node2=node
            endif    
            if(ipe(node1).eq.0) then
               ifreeme=ifreeme+1
               ipe(node1)=ifreeme
               ime(1,ifreeme)=node2
               ime(2,ifreeme)=j
               ime(3,ifreeme)=ipos
            else
               index1=ipe(node1)
               if(ime(1,index1).eq.node2) then
                  imastop(ipos,j)=ime(2,index1)
                  imastop(ime(3,index1),ime(2,index1))=j
                  cycle
               endif
!
               index1old=index1
               index1=ime(4,index1)
               do
                  if(index1.eq.0) then
                     ifreeme=ifreeme+1
                     ime(4,index1old)=ifreeme
                     ime(1,ifreeme)=node2
                     ime(2,ifreeme)=j
                     ime(3,ifreeme)=ipos 
                     exit
                  endif
                  if(ime(1,index1).eq.node2) then
                     imastop(ipos,j)=ime(2,index1)
                     imastop(ime(3,index1),ime(2,index1))=j
c                     ime(4,index1old)=ime(4,index1)
                     exit
                  endif
                  index1old=index1
                  index1=ime(4,index1)
               enddo
            endif
         enddo
      enddo
!
      return
      end
