!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine machpi (MACH, PI,kappa, rgas)
!
!-----------------------------------------------------------------------
!                                                                      |
!     Dieses Unterprogramm berechnet die Mach-Zahl fuer das            |
!     eingegebene Druckverhaeltnis PI.                                 |
!                                                                      |
!     Eingabe-Groessen:                                                |
!       PI     = Druckverhaeltnis PS/PT                                |
!                                                                      |
!     Ausgabe-Groessen:                                                |
!       MACH   = Mach-Zahl                                             |
!                                                                      |
!-----------------------------------------------------------------------
!
      IMPLICIT CHARACTER*1 (A-Z)
      real*8    PI, MACH, MA2, kappa, rgas, kappam,kappax,pikrit
!
!-----------------------------------------------------------------------
!
      kappax = (kappa-1)/kappa
      KAPPAM = 2.d0 / (KAPPA - 1.d0)
      PIKRIT = (2.d0/(KAPPA+1.d0)) ** (KAPPA/(KAPPA-1.d0))
!
      IF (PI.GE.1.d0) THEN
!       Druckverhaeltnis groesser gleich 1
        MACH = 0.d0
      ELSEIF (PI.GT.PIKRIT) THEN
!       Druckverhaeltnis unterkritisch
        MA2  = KAPPAM * (PI**(-KAPPAX) - 1.d0)
        IF (MA2.GT.0.d0) THEN
          MACH = SQRT (MA2)
        ELSE
          MACH = 0.d0
        ENDIF
      ELSEIF (PI.GT.0.d0) THEN
!       Druckverhaeltnis ueberkritisch
        MACH = 1.d0
      ELSE
!       Druckverhaeltnis ungueltig
        MACH = 1.d20
      ENDIF
!
      RETURN
      END
!
