!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine temploadmodal(amta,namta,nam,ampli,time,ttime,dtime,
     &  xbounold,xboun,xbounact,iamboun,nboun,nodeboun,ndirboun,
     &  amname,reltime)
!
!     calculates the SPC boundary conditions at a given time for
!     a modal dynamic procedure; used to calculate the velocity and
!     acceleration by use of finite differences
!
      implicit none
!
      character*80 amname(*)
!
      integer nam,i,istart,iend,id,namta(3,*),
     &  iamboun(*),nboun,iambouni,nodeboun(*),ndirboun(*)
!
      real*8 amta(2,*),ampli(*),time,reltime,
     &  xbounold(*),xboun(*),xbounact(*),ttime,dtime,reftime
!
!     if an amplitude is active, the loading is scaled according to
!     the actual time. If no amplitude is active, then the load is
!     applied as a step loading
!
!     calculating all amplitude values for the current time
!
      do i=1,nam
         if(namta(3,i).lt.0) then
            reftime=ttime+time
         else
            reftime=time
         endif
         if(abs(namta(3,i)).ne.i) then
            reftime=reftime-amta(1,namta(1,i))
            istart=namta(1,abs(namta(3,i)))
            iend=namta(2,abs(namta(3,i)))
            if(istart.eq.0) then
               call uamplitude(reftime,amname(namta(3,i)),ampli(i))
               cycle
            endif
         else
            istart=namta(1,i)
            iend=namta(2,i)
            if(istart.eq.0) then
               call uamplitude(reftime,amname(i),ampli(i))
               cycle
            endif
         endif
         call identamta(amta,reftime,istart,iend,id)
         if(id.lt.istart) then
            ampli(i)=amta(2,istart)
         elseif(id.eq.iend) then
            ampli(i)=amta(2,iend)
         else
            ampli(i)=amta(2,id)+(amta(2,id+1)-amta(2,id))
     &           *(reftime-amta(1,id))/(amta(1,id+1)-amta(1,id))
         endif
      enddo
!
!     scaling the boundary conditions
!
      do i=1,nboun
         if(nam.gt.0) then
            iambouni=iamboun(i)
         else
            iambouni=0
         endif
         if(iambouni.gt.0) then
            if(amname(iambouni).eq.'RAMP12357111317') then
               xbounact(i)=xbounold(i)+
     &              (xboun(i)-xbounold(i))*reltime
            else
               xbounact(i)=xboun(i)*ampli(iambouni)
            endif
         else
            xbounact(i)=xboun(i)
         endif
      enddo
!
      return
      end
