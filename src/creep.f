!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine creep(decra,deswa,statev,serd,ec,esw,p,qtild,
     &  temp,dtemp,predef,dpred,time,dtime,cmname,leximp,lend,
     &  coords,nstatv,noel,npt,layer,kspt,kstep,kinc)
!
!     user creep routine
!
!     INPUT (general):
!
!     statev(1..nstatv)  internal variables
!     serd               not used
!     ec(1)              equivalent creep at the start of the increment
!     ec(2)              not used
!     esw(1..2)          not used
!     p                  not used
!     temp               temperature at the end of the increment
!     dtemp              not used
!     predef             not used
!     dpred              not used
!     time(1)            value of the step time at the end of the increment
!     time(2)            value of the total time at the end of the increment
!     dtime              time increment
!     cmname             material name
!     leximp             not used
!     lend               if = 2: isotropic creep
!                        if = 3: anisotropic creep
!     coords(1..3)       coordinates of the current integration point
!     nstatv             number of internal variables
!     noel               element number
!     npt                integration point number
!     layer              not used
!     kspt               not used
!     kstep              not used
!     kinc               not used
!
!    INPUT only for elastic isotropic materials:
!     qtild              von Mises stress
!
!    INPUT only for elastic anisotropic materials:
!     decra(1)           equivalent deviatoric creep strain increment
!
!
!     OUTPUT (general):
!
!     decra(5)           derivative of the equivalent deviatoric
!                        creep strain increment w.r.t. the von Mises
!                        stress
!
!     OUTPUT only for elastic isotropic materials:
!     decra(1)           equivalent deviatoric creep strain increment
!
!     OUTPUT only for elastic anisotropic materials:
!     qtild              von Mises stress
!
      implicit none
!   
      character*80 cmname
!
      integer leximp,lend,nstatv,noel,npt,layer,kspt,kstep,kinc
!
      real*8 decra(5),deswa(5),statev(*),serd,ec(2),esw(2),p,qtild,
     &  temp,dtemp,predef(*),dpred(*),time(2),dtime,coords(*)
!
      qtild=(1.d10*decra(1)/dtime)**0.2d0
      decra(5)=5.d-10*dtime*qtild**4
!
      return
      end
