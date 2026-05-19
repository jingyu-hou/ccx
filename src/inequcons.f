!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine inequcons(inpc,textpart,istep,istat,n,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc,inequcon,inequbou)        
!
!     reading the input deck: *INEQUALITY CONSTRAINTS
!
!     criteria: MASS
!               STRESS
!               DISPLACEMENT
!            
      implicit none
!
      character*1 inpc(*)
      character*132 textpart(16)
      character*81 inequcon(3)
!
      integer istep,istat,n,key,i,iline,ipol,inl,ipoinp(2,*),
     &  inp(3,*),ipoinpc(0:*)
!
      real*8 inequbou
!
      if(istep.lt.1) then
         write(*,*) '*ERROR reading *INEQUALITY CONSTRAINTS: 
     &*INEQUALITY CONSTRAINTS can only be used within a 
     &SENSITIVITY STEP'     
         call exit(201)
      endif
!
      do i=2,n
         if(textpart(i)(1:9).eq.'CRITERIA=') then
            read(textpart(i)(10:85),'(a80)',iostat=istat) 
     &           inequcon(1)(1:80)
         elseif(textpart(i)(1:7).eq.'ENTITY=') then
            read(textpart(i)(8:85),'(a80)',iostat=istat) 
     &           inequcon(2)(1:80) 
         elseif(textpart(i)(1:3).eq.'MIN') then
            read(textpart(i)(1:85),'(a80)',iostat=istat) 
     &           inequcon(3)(1:80) 
         endif
      enddo
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
!
      read(textpart(1)(1:20),'(f20.0)',iostat=istat) inequbou
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
!
      return
      end

