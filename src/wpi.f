!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine wpi(W, PI, Q, SQTT,kappa,RGAS)
!
!-----------------------------------------------------------------------
!                                                                      |
!     Dieses Unterprogramm berechnet die Stroemungs-Geschwindigkeit    |
!     fuer das eingegebene Druckverhaeltnis PI.                        |
!                                                                      |
!     Eingabe-Groessen:                                                |
!       PI     = Druckverhaeltnis PS/PT                                |
!       Q      = reduzierter Durchsatz                                 |
!       SQTT   = SQRT (Totaltemperatur)                                |
!                                                                      |
!     Ausgabe-Groessen:                                                |
!       W      = Stroemungs-Geschwindigkeit                            |
!                                                                      |
!-----------------------------------------------------------------------
!
      IMPLICIT CHARACTER*1 (A-Z)
!       INCLUDE 'comkapfk.inc'
      real*8    W, PI, Q, SQTT,kappaq,kappa,RGAS,pikrit,kappah,wkritf
!
!-----------------------------------------------------------------------
!
      kappaq = 1.d0/kappa
      PIKRIT = (2.d0/(KAPPA+1.d0)) ** (KAPPA/(KAPPA-1.d0))
!
      KAPPAH = 2.d0 * KAPPA / (KAPPA + 1.d0)
      WKRITF = SQRT( KAPPAH * RGAS )
!
      IF (PI.GE.1.d0) THEN
!       Druckverhaeltnis groesser gleich 1
        W    = 0.d0
      ELSEIF (PI.GT.PIKRIT) THEN
!       Druckverhaeltnis unterkritisch
        IF (Q.GT.0.d0) THEN
          W    = Q * RGAS * SQTT * PI**(-KAPPAQ)
        ELSE
          W    = 0.d0
        ENDIF
      ELSEIF (PI.GT.0.d0) THEN
!       Druckverhaeltnis ueberkritisch
        W    = WKRITF * SQTT
      ELSE
!       Druckverhaeltnis ungueltig
        W    = 1.d20
      ENDIF
!
      RETURN
      END
      
      
      
