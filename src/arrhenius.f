      
C     WeICME (Wedge Integrated Computational Materials Engineering)
C                 - A 3-dimensional finite element program.
C     
C     Developed and maintained by Shenzhen Wedge Central 
C     South Research Institute co., Ltd., Shenzhen, China
C     
C     Copy Right 2019-2023.
C
      subroutine uhard(syield,hard,eqplas,eqplasrt,temp,numprops,
     $   props)
C
C       INCLUDE 'ABA_PARAM.INC'
C
      integer numprops
      real*8 syield,hard(3),eqplas,eqplasrt,temp,props(numprops)
c     
      real*8 alpha,a,q,n,x,t,z
      integer i,j
      real*8, parameter :: r=8.314d0,one=1.0d0

      
      alpha=props(1)
      a=props(2)
      q=props(3)
      n=props(4)
      t=temp+273.d0
c
        
      if (eqplasrt .lt. 0.001d0) then
          epsilon_rate=0.001d0
          x=(epsilon_rate*exp(q/r/t)/a)**(one/n)
          hard(2)=0.d0
      elseif (eqplasrt .gt. 1.d0) then
          epsilon_rate=1.d0
          x=(epsilon_rate*exp(q/r/t)/a)**(one/n)
          hard(2)=0.d0
      else
          epsilon_rate=eqplasrt
          x=(epsilon_rate*exp(q/r/t)/a)**(one/n)
          hard(2)=x/(epsilon_rate*alpha*n*(1.0d0+x*x)**(0.5d0))*1.d06
      endif

      syield=one/alpha*asinh(x)*1.d06
      hard(1)=0.0d0
      hard(3)=-q*x/(r*t**2*alpha*n*(1.0d0+x*x)**(0.5d0))*1.d06

      return
      end

