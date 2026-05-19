!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine springstiff_n2f_entha(xl,voldl,s,imat,elcon,nelcon,
     &  ncmat_,ntmat_,nope,kode,plkcon,
     &  nplkcon,npmat_,iperturb,springarea,mi,timeend,matname,
     &  node,noel,istep,iinc,ithermal,enthaoldl,htcon,nhtcon,ipairmat)
!
!     calculates the stiffness of a spring
!
      implicit none
!
      character*80 matname(*),slname,msname
!
      integer i,j,imat,ncmat_,ntmat_,k,nope,nterms,iflag,
     &  kode,niso,id,nplkcon(0:ntmat_,*),npmat_,nelcon(2,*),
     &  iperturb(*),mi(*),node,noel,istep,iinc,npred,ithermal(2),
     &  nhtcon(*),ipairmat(0:2,*),imat2
!
      real*8 xl(3,10),ratio(9),q(3),val,shp2(7,9),ak(5),
     &  al(3),s(60,60),voldl(0:mi(2),10),pl(3,10),xn(3),dm,
     &  alpha,beta,elcon(0:ncmat_,ntmat_,*),xm(3),pressure,
     &  xi,et,xs2(3,7),t1l,elconloc(21),plconloc(802),xk,
     &  xiso(200),yiso(200),plkcon(0:2*npmat_,ntmat_,*),
     &  springarea,dist,eps,pi,constant,conductance,dtemp,temp(2),
     &  predef(2),coords(3),tmean,d(2),timeend(2),flowm(2),
     &  enthaoldl(*),htcon(0:1,nhtcon(1),*),entha_n,entha_f,dtdh,
     &  entha
!
      intent(in) xl,voldl,imat,elcon,nelcon,
     &  ncmat_,ntmat_,nope,kode,plkcon,
     &  nplkcon,npmat_,iperturb,mi,timeend,matname,
     &  node,noel,istep,iinc
!
      intent(inout) s,springarea
!
      iflag=4
!
!     actual positions of the nodes belonging to the contact spring
!
      if(iperturb(2).eq.0) then
         do i=1,nope
            do j=1,3
               pl(j,i)=xl(j,i)
            enddo
         enddo
      else
         do i=1,nope
            do j=1,3
               pl(j,i)=xl(j,i)+voldl(j,i)
            enddo
         enddo
      endif
!
!     contact springs
!
      nterms=nope-1
!
!     vector al connects the actual position of the slave node 
!     with its projection on the master face
!
      do i=1,3
         q(i)=pl(i,nope)
      enddo
      call attach(pl,q,nterms,ratio,dist,xi,et)
      do i=1,3
         al(i)=pl(i,nope)-q(i)
      enddo
!
!     determining the jacobian vector on the surface 
!
c      if(nterms.eq.9) then
c         call shape9q(xi,et,pl,xm,xs2,shp2,iflag)
      if(nterms.eq.8) then
         call shape8q(xi,et,pl,xm,xs2,shp2,iflag)
      elseif(nterms.eq.4) then
         call shape4q(xi,et,pl,xm,xs2,shp2,iflag)
      elseif(nterms.eq.6) then
         call shape6tri(xi,et,pl,xm,xs2,shp2,iflag)
c      elseif(nterms.eq.7) then
c         call shape7tri(xi,et,pl,xm,xs2,shp2,iflag)
      else
         call shape3tri(xi,et,pl,xm,xs2,shp2,iflag)
      endif
!
!
!     normal on the surface
!
      dm=dsqrt(xm(1)*xm(1)+xm(2)*xm(2)+xm(3)*xm(3))
      do i=1,3
         xn(i)=xm(i)/dm
      enddo
!
!     distance from surface along normal
!
      val=al(1)*xn(1)+al(2)*xn(2)+al(3)*xn(3)
!
!     representative area: usually the slave surface stored in
!     springarea; however, if no area was assigned because the
!     node does not belong to any element, the master surface
!     is used
!
      if(springarea.le.0.d0) then
         if(nterms.eq.3) then
            springarea=dm/2.d0
         else
            springarea=dm*4.d0
         endif
      endif
!
!     alpha and beta, taking the representative area into account
!     (conversion of pressure into force)
!
      if(elcon(1,1,imat).gt.0.d0) then
!
!        exponential overclosure
!
         if(dabs(elcon(2,1,imat)).lt.1.d-30) then
            pressure=0.d0
         else
            alpha=elcon(2,1,imat)
            beta=elcon(1,1,imat)
            if(-beta*val.gt.23.d0-dlog(alpha)) then
               beta=(dlog(alpha)-23.d0)/val
            endif
            pressure=dexp(-beta*val+dlog(alpha))
         endif
      else
!     
!        linear overclosure
!
         pi=4.d0*datan(1.d0)
         eps=elcon(1,1,imat)*pi/elcon(2,1,imat)
         pressure=-elcon(2,1,imat)*val*
     &            (0.5d0+datan(-val/eps)/pi)
      endif
!
!     calculating the temperature difference across the contact
!     area and the mean temperature for the determination of the
!     conductance
!
      t1l=0.d0
      entha_f=0.d0
      do j=1,nterms
         t1l=t1l+ratio(j)*voldl(0,j)
         entha_f=entha_f+ratio(j)*enthaoldl(j)
      enddo
      dtemp=t1l-voldl(0,nope)
c      tmean=(voldl(0,nope)+t1l)/2.d0
      tmean=t1l
      entha_n=enthaoldl(nope)
!
!     interpolating the material data according to temperature
!
      call materialdata_sp(elcon,nelcon,imat,ntmat_,i,tmean,
     &     elconloc,kode,plkcon,nplkcon,npmat_,plconloc,ncmat_)
!
!     interpolating the material data according to pressure
!
      niso=int(plconloc(801))
!
      if(niso.eq.0) then
!
!        user subroutine for the conductance
!
         d(1)=val
         d(2)=pressure
         temp(1)=voldl(0,nope)
         temp(2)=t1l
         do j=1,3
            coords(j)=xl(j,nope)
         enddo
         call gapcon(ak,d,flowm,temp,predef,timeend,matname(imat),
     &               slname,msname,coords,noel,node,npred,istep,iinc,
     &               springarea)
         conductance=ak(1)
      else
         do i=1,niso
            xiso(i)=plconloc(2*i-1)
            yiso(i)=plconloc(2*i)
         enddo
         call ident(xiso,pressure,niso,id)
         if(id.eq.0) then
            xk=0.d0
            conductance=yiso(1)
         elseif(id.eq.niso) then
            xk=0.d0
            conductance=yiso(niso)
         else
            xk=(yiso(id+1)-yiso(id))/(xiso(id+1)-xiso(id))
            conductance=yiso(id)+xk*(pressure-xiso(id))
         endif
      endif
!
!     assembling the upper triangle of the element matrix
!
      constant=conductance*springarea
!
c      call materialdata_dtdh(htcon,nhtcon,imat,dtdh,
c     &                 entha_n,ntmat_,ithermal)
!
c      s(nope,nope)=constant*dtdh
c      do k=1,nterms
c         dtdh=dabs(enthaoldl(k)-enthaoldl(nope))
c         if(dtdh.lt.1.d-2) then
c            entha=enthaoldl(k)
c            call materialdata_dtdh(htcon,nhtcon,imat,dtdh,
c     &                 entha,ntmat_,ithermal)
c         else
c            dtdh=dabs(voldl(0,k)-voldl(0,nope))/dtdh
c         endif
c         s(k,nope)=-shp2(4,k)*constant*dtdh
c      enddo
c      do i=1,nterms
c         do j=i,nterms
c            dtdh=dabs(enthaoldl(j)-enthaoldl(i))
c            if(dtdh.lt.1.d-2) then
c               entha=enthaoldl(j)
c               call materialdata_dtdh(htcon,nhtcon,imat,dtdh,
c     &                    entha,ntmat_,ithermal)
c            else
c               dtdh=dabs(voldl(0,j)-voldl(0,i))/dtdh
c            endif            
c            s(i,j)=shp2(4,i)*shp2(4,j)*constant*dtdh
c         enddo
c      enddo
!
      dtdh=0.d0
      imat2=ipairmat(2,imat)
      call materialdata_dtdh(htcon,nhtcon,imat2,dtdh,
     &                 entha_n,ntmat_,ithermal)
c      print *,"springstiff_n2f_entha",imat,imat2,dtdh,entha_n
!
      s(nope,nope)=constant*dtdh
      do k=1,nterms
         s(k,nope)=-shp2(4,k)*constant*dtdh
      enddo
!
      imat2=ipairmat(1,imat)
      call materialdata_dtdh(htcon,nhtcon,imat2,dtdh,
     &                 entha_f,ntmat_,ithermal)
c      print *,"springstiff_n2f_entha",imat,imat2,dtdh,entha_f
!
      do i=1,nterms
         do j=i,nterms           
            s(i,j)=shp2(4,i)*shp2(4,j)*constant*dtdh
         enddo
      enddo
!
      return
      end

