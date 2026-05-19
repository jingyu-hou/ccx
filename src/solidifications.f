!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine solidifications(inpc,textpart,nmethod,iperturb,isolver,
     &  istep,istat,n,tinc,tper,tmin,tmax,idrct,ithermal,iline,ipol,
     &  inl,ipoinp,inp,alpha,mei,fei,ipoinpc,ctrl,ttime,ier,enthactrl)
!
!     reading the input deck: *SOLIDIFICATION
!
!     isolver=0: SPOOLES
!             2: iterative solver with diagonal scaling
!             3: iterative solver with Cholesky preconditioning
!             4: sgi solver
!             5: TAUCS
!             7: pardiso
!
      implicit none
!
      logical timereset
!
      character*1 inpc(*)
      character*20 solver
      character*132 textpart(16)
!
      integer nmethod,iperturb,isolver,istep,istat,n,key,i,idrct,nev,
     &  ithermal,iline,ipol,inl,ipoinp(2,*),inp(3,*),mei(4),ncv,mxiter,
     &  ipoinpc(0:*),idirect,ier
!
      real*8 tinc,tper,tmin,tmax,alpha,fei(3),tol,fmin,fmax,ctrl(*),
     &  ttime,enthactrl(*)
!
      tmin=0.d0
      tmax=0.d0
      nmethod=4
      alpha=0.d0
      mei(4)=0
      timereset=.false.
!
!     enthactrl(1): value 1 means enthalpy based method
!     enthactrl(2): temporal discretization scheme
!                   1:explicit
!                   2:implicit with frozen marix element
!                   3:implicit with Newton-Raphson
!     enthactrl(3): criterion deciding when ny should be calculated
!     enthactrl(4): criterion deciding when ny* should be calculated
!     enthactrl(5): C_lambda
!     enthactrl(6): viscosity of liquidus
!     
      enthactrl(1)=1.d0
      enthactrl(2)=-1.d0
      enthactrl(3)=-1.d0
      enthactrl(4)=-1.d0
      enthactrl(5)=-1.d0
      enthactrl(6)=-1.d0
!
!     defaults for fmin and fmax
!
      fmin=-1.d0
      fmax=-1.d0
!
      if(iperturb.eq.0) then
         iperturb=2
      elseif((iperturb.eq.1).and.(istep.gt.1)) then
         write(*,*) 
     &     '*ERROR reading *SOLIDIFICATION: perturbation analysis is'
         write(*,*) '       not provided in a *SOLIDIFICATION step.'
         ier=1
         return
      endif
!
      if(istep.lt.1) then
         write(*,*) 
     &     '*ERROR reading *SOLIDIFICATION: *SOLIDIFICATION can only'
         write(*,*) '       be used within a STEP'
         ier=1
         return
      endif
!
!     default solver
!
      solver='                    '
      if(isolver.eq.0) then
         solver(1:7)='SPOOLES'
      elseif(isolver.eq.2) then
         solver(1:16)='ITERATIVESCALING'
      elseif(isolver.eq.3) then
         solver(1:17)='ITERATIVECHOLESKY'
      elseif(isolver.eq.4) then
         solver(1:3)='SGI'
      elseif(isolver.eq.5) then
         solver(1:5)='TAUCS'
      elseif(isolver.eq.7) then
         solver(1:7)='PARDISO'
      endif
!
      idirect=2
      do i=2,n
         if(textpart(i)(1:7).eq.'SOLVER=') then
            read(textpart(i)(8:27),'(a20)') solver
         elseif((textpart(i)(1:6).eq.'DIRECT').and.
     &          (textpart(i)(1:9).ne.'DIRECT=NO')) then
            idirect=1
         elseif(textpart(i)(1:9).eq.'DIRECT=NO') then
            idirect=0
         elseif(textpart(i)(1:11).eq.'STEADYSTATE') then
            nmethod=1
         elseif(textpart(i)(1:9).eq.'FREQUENCY') then
            nmethod=2
         elseif(textpart(i)(1:12).eq.'MODALDYNAMIC') then
            iperturb=0
         elseif(textpart(i)(1:11).eq.'STORAGE=YES') then
            mei(4)=1
         elseif(textpart(i)(1:7).eq.'DELTMX=') then
            read(textpart(i)(8:27),'(f20.0)',iostat=istat) ctrl(27)
         elseif(textpart(i)(1:9).eq.'TIMERESET') then
            timereset=.true.
         elseif(textpart(i)(1:17).eq.'TOTALTIMEATSTART=') then
            read(textpart(i)(18:37),'(f20.0)',iostat=istat) ttime
         elseif(textpart(i)(1:8).eq.'EXPLICIT') then
            enthactrl(2)=1.d0
         elseif(textpart(i)(1:10).eq.'IMPLICIT') then
            enthactrl(2)=2.d0
         elseif(textpart(i)(1:10).eq.'NRIMPLICIT') then
            enthactrl(2)=3.d0
         else
            write(*,*) 
     &   '*WARNING reading *SOLIDIFICATION: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*SOLIDIFICATION%")
         endif
      enddo
!
      if(enthactrl(2).lt.0.d0) then
         enthactrl(2)=1.d0
         write(*,*) 
     &     '*WARNING reading *SOLIDIFICATION: time scheme is not '
         write(*,*) '       specified, options are EXPLICIT, '
         write(*,*) '       IMPLICIT, NRIMPLICIT. EXPLICIT is'
         write(*,*) '       autolly selected for simulation now'
      endif
!
      if(nmethod.eq.1) ctrl(27)=1.d30
!
!     default for modal dynamic calculations is DIRECT,
!     for static or dynamic calculations DIRECT=NO
!
      if(iperturb.eq.0) then
         idrct=1
         if(idirect.eq.0)idrct=0
      else
         idrct=0
         if(idirect.eq.1)idrct=1
      endif
!
      if((ithermal.eq.0).and.(nmethod.ne.1).and.
     &   (nmethod.ne.2).and.(iperturb.ne.0)) then
         write(*,*) 
     &     '*ERROR reading *SOLIDIFICATION: please define initial '
         write(*,*) '       conditions for the temperature'
         ier=1
         return
      else
         ithermal=2
      endif
!
      if((nmethod.ne.2).and.(iperturb.ne.0)) then
!
!        static or dynamic thermal analysis
!
         if(solver(1:7).eq.'SPOOLES') then
            isolver=0
         elseif(solver(1:16).eq.'ITERATIVESCALING') then
            isolver=2
         elseif(solver(1:17).eq.'ITERATIVECHOLESKY') then
            isolver=3
         elseif(solver(1:3).eq.'SGI') then
            isolver=4
         elseif(solver(1:5).eq.'TAUCS') then
            isolver=5
         elseif(solver(1:7).eq.'PARDISO') then
            isolver=7
         else
            write(*,*) 
     &          '*WARNING reading *SOLIDIFICATION: unknown solver;'
            write(*,*) '         the default solver is used'
         endif
!
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) then
            if(iperturb.ge.2) then
               write(*,*) 
     &             '*WARNING reading *SOLIDIFICATION: a nonlinear geomet
     &ric analysis is requested'
               write(*,*) '         but no time increment nor step is sp
     &ecified'
               write(*,*) '         the defaults (1,1) are used'
               tinc=1.d0
               tper=1.d0
               tmin=1.d-5
               tmax=1.d+30
            endif
            if(timereset)ttime=ttime-tper
            return
         endif
!
         read(textpart(1)(1:20),'(f20.0)',iostat=istat) tinc
         if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*SOLIDIFICATION%",ier)
            return
         endif
         read(textpart(2)(1:20),'(f20.0)',iostat=istat) tper
         if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*SOLIDIFICATION%",ier)
            return
         endif
         read(textpart(3)(1:20),'(f20.0)',iostat=istat) tmin
         if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*SOLIDIFICATION%",ier)
            return
         endif
         read(textpart(4)(1:20),'(f20.0)',iostat=istat) tmax
         if(istat.gt.0) then
            call inputerror(inpc,ipoinpc,iline,
     &           "*SOLIDIFICATION%",ier)
            return
         endif
!
         if(tinc.le.0.d0) then
            write(*,*) 
     &        '*ERROR reading *SOLIDIFICATION: initial increment size 
     &is negative'
         endif
         if(tper.le.0.d0) then
            write(*,*) 
     &        '*ERROR reading *SOLIDIFICATION: step size is negative'
         endif
         if(tinc.gt.tper) then
            write(*,*) 
     &         '*ERROR reading *SOLIDIFICATION: initial increment size 
     &exceeds step size'
         endif
!      
         if(idrct.ne.1) then
c            if(dabs(tmin).lt.1.d-10) then
c               tmin=min(tinc,1.d-5*tper)
            if(dabs(tmin).lt.1.d-16*tper) then
               tmin=min(tinc,1.d-16*tper)
            endif
            if(dabs(tmax).lt.1.d-10) then
               tmax=1.d+30
            endif
         endif
      else
!        if..else.. here maybe not necessary,
!        reserved once bugs appears ---20210830 by Ye
         write(*,*) 
     &     '*ERROR reading *SOLIDIFICATION: nmethod '
         write(*,*) '       and iperturb is not suitable'
         call exit(201)
      endif
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
      if((istat.lt.0).or.(key.eq.1)) then
         write(*,*) '*WARNING reading *SOLIDIFICATION: Niyama'
         write(*,*) '       criterion is missing, 0.9 is specifed'
         write(*,*) '       for both ny and ny*'
         enthactrl(3)=0.9d0
         enthactrl(4)=0.9d0
         enthactrl(5)=4.09d-5
         enthactrl(6)=1.47d-3
         call inputerror(inpc,ipoinpc,iline,
     &        "*SOLIDIFICATION%",ier)
         return
      endif
!
      read(textpart(1)(1:20),'(f20.0)',iostat=istat) enthactrl(3)
      if(istat.gt.0) then
         call inputerror(inpc,ipoinpc,iline,
     &        "*SOLIDIFICATION%",ier)
c         return
      endif
      read(textpart(2)(1:20),'(f20.0)',iostat=istat) enthactrl(4)
      if(istat.gt.0) then
         call inputerror(inpc,ipoinpc,iline,
     &        "*SOLIDIFICATION%",ier)
c         return
      endif
      read(textpart(3)(1:20),'(f20.0)',iostat=istat) enthactrl(5)
      if(istat.gt.0) then
         call inputerror(inpc,ipoinpc,iline,
     &        "*SOLIDIFICATION%",ier)
c         return
      endif
      read(textpart(4)(1:20),'(f20.0)',iostat=istat) enthactrl(6)
      if(istat.gt.0) then
         call inputerror(inpc,ipoinpc,iline,
     &        "*SOLIDIFICATION%",ier)
c         return
      endif
!
      if(enthactrl(3).lt.1.d-10) then
         enthactrl(3)=0.9d0
         write(*,*) 
     &     '*WARNING reading *SOLIDIFICATION: Niyama criterion'
         write(*,*) '       is missing or negative, 0.9 is specifed'
      elseif(enthactrl(3).gt.1.d0) then
         enthactrl(3)=0.9d0
         write(*,*) 
     &     '*ERROR reading *SOLIDIFICATION: Niyama criterion'
         write(*,*) '       larger than 1, 0.9 is specifed'
      endif
!
      if(enthactrl(4).lt.1.d-10) then
         enthactrl(4)=0.9d0
         write(*,*) 
     &     '*WARNING reading *SOLIDIFICATION: Niyama star criterion'
         write(*,*) '       is missing or negative, 0.9 is specifed'
      elseif(enthactrl(4).gt.1.d0) then
         enthactrl(4)=0.9d0
         write(*,*) 
     &     '*ERROR reading *SOLIDIFICATION: Niyama star criterion'
         write(*,*) '       larger than 1, 0.9 is specifed'
      endif
!
      if(enthactrl(5).lt.1.d-10) then
         enthactrl(5)=4.09d-5
         write(*,*) 
     &     '*WARNING reading *SOLIDIFICATION: C_lambda of alloy'
         write(*,*) '       is missing or negative, 4.09d-5 is specifed'
      endif
!
      if(enthactrl(6).lt.1.d-10) then
         enthactrl(6)=1.47d-3
         write(*,*) 
     &     '*WARNING reading *SOLIDIFICATION: viscosity of alloy'
         write(*,*) '       is missing or negative'
      endif
!
      if(timereset)ttime=ttime-tper
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end








