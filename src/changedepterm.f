!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine changedepterm(ikmpc,ilmpc,nmpc,mpc,idofrem,idofins)
!
!     changes the dependent term in ikmpc and ilmpc for MPC mpc.
!
      implicit none
!
      integer ikmpc(*),ilmpc(*),nmpc,idofrem,idofins,id,k,mpc
!
!     remove MPC from ikmpc
!
      call nident(ikmpc,idofrem,nmpc,id)
      if(id.gt.0) then
         if(ikmpc(id).eq.idofrem) then
            do k=id+1,nmpc
               ikmpc(k-1)=ikmpc(k)
               ilmpc(k-1)=ilmpc(k)
            enddo
         else
            write(*,*) '*ERROR in changedepterm'
            write(*,*) '       ikmpc database corrupted'
            call exit(201)
         endif
      else
         write(*,*) '*ERROR in changedepterm'
         write(*,*) '       ikmpc database corrupted'
         call exit(201)
      endif
!
!     insert new MPC
!
      call nident(ikmpc,idofins,nmpc-1,id)
      if((id.gt.0).and.(ikmpc(id).eq.idofins)) then
         write(*,*) '*ERROR in changedepterm: dependent DOF'
         write(*,*) '       of nonlinear MPC cannot be changed'
         write(*,*) '       since new dependent DOF is already'
         write(*,*) '       used in another MPC'
         call exit(201)
      else
         do k=nmpc,id+2,-1
            ikmpc(k)=ikmpc(k-1)
            ilmpc(k)=ilmpc(k-1)
         enddo
         ikmpc(id+1)=idofins
         ilmpc(id+1)=mpc
      endif
!
      return
      end











