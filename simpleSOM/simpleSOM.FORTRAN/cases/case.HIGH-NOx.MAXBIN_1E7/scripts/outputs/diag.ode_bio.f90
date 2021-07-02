!     ==========================================================================
!                 THIS SUBROUTINE DEFINES ODEs FOR BIO-REGIME CHEMISTRY
!     ==========================================================================

      SUBROUTINE ODE_BIO
  
      USE mod_MAIN
      USE mod_GAS

      IMPLICIT NONE
      
      REAL(R8) :: OXY_AVE

!FLAG1
      p_BIO(iCN3) = &
      1.02e-03*r_BIO(1) + &
      2.04e-03*r_BIO(2) + &
      2.51e-03*r_BIO(3) + &
      3.39e-03*r_BIO(4) + &
      6.20e-03*r_BIO(5) + &
      6.13e-03*r_BIO(6) + &
      1.40e-03*r_BIO(7) + &
      4.74e-05*r_BIO(8) + &
      2.20e-07*r_BIO(9) + &
      0.00e+00*r_BIO(10) + &
      0.00e+00*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      d_BIO(iCN3) = r_BIO(1)

      p_BIO(iCN2) = &
      4.99e-01*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      1.54e-03*r_BIO(3) + &
      1.61e-03*r_BIO(4) + &
      5.40e-04*r_BIO(5) + &
      6.21e-03*r_BIO(6) + &
      1.08e-02*r_BIO(7) + &
      2.73e-03*r_BIO(8) + &
      9.41e-05*r_BIO(9) + &
      4.39e-07*r_BIO(10) + &
      0.00e+00*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      d_BIO(iCN2) = r_BIO(2)

      p_BIO(iCN1) = &
      4.99e-01*r_BIO(1) + &
      4.99e-01*r_BIO(2) + &
      0.00e+00*r_BIO(3) + &
      3.07e-03*r_BIO(4) + &
      3.21e-03*r_BIO(5) + &
      1.08e-03*r_BIO(6) + &
      1.24e-02*r_BIO(7) + &
      2.15e-02*r_BIO(8) + &
      5.44e-03*r_BIO(9) + &
      1.87e-04*r_BIO(10) + &
      8.74e-07*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      d_BIO(iCN1) = r_BIO(3)

      p_BIO(iC0) = &
      0.00e+00*r_BIO(1) + &
      4.99e-01*r_BIO(2) + &
      4.98e-01*r_BIO(3) + &
      0.00e+00*r_BIO(4) + &
      6.12e-03*r_BIO(5) + &
      6.40e-03*r_BIO(6) + &
      2.14e-03*r_BIO(7) + &
      2.46e-02*r_BIO(8) + &
      4.28e-02*r_BIO(9) + &
      1.08e-02*r_BIO(10) + &
      3.73e-04*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      d_BIO(iC0) = r_BIO(4)

      p_BIO(iC1) = &
      0.00e+00*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      4.98e-01*r_BIO(3) + &
      4.96e-01*r_BIO(4) + &
      0.00e+00*r_BIO(5) + &
      1.22e-02*r_BIO(6) + &
      1.27e-02*r_BIO(7) + &
      4.26e-03*r_BIO(8) + &
      4.90e-02*r_BIO(9) + &
      8.53e-02*r_BIO(10) + &
      2.15e-02*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      d_BIO(iC1) = r_BIO(5)

      p_BIO(iC2) = &
      0.00e+00*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      0.00e+00*r_BIO(3) + &
      4.96e-01*r_BIO(4) + &
      4.92e-01*r_BIO(5) + &
      0.00e+00*r_BIO(6) + &
      2.42e-02*r_BIO(7) + &
      2.53e-02*r_BIO(8) + &
      8.49e-03*r_BIO(9) + &
      9.75e-02*r_BIO(10) + &
      1.70e-01*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      d_BIO(iC2) = r_BIO(6)

      p_BIO(iC3) = &
      0.00e+00*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      0.00e+00*r_BIO(3) + &
      0.00e+00*r_BIO(4) + &
      4.92e-01*r_BIO(5) + &
      4.84e-01*r_BIO(6) + &
      0.00e+00*r_BIO(7) + &
      4.83e-02*r_BIO(8) + &
      5.05e-02*r_BIO(9) + &
      1.69e-02*r_BIO(10) + &
      1.94e-01*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      d_BIO(iC3) = r_BIO(7)

      p_BIO(iC4) = &
      0.00e+00*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      0.00e+00*r_BIO(3) + &
      0.00e+00*r_BIO(4) + &
      0.00e+00*r_BIO(5) + &
      4.84e-01*r_BIO(6) + &
      4.68e-01*r_BIO(7) + &
      0.00e+00*r_BIO(8) + &
      9.61e-02*r_BIO(9) + &
      1.00e-01*r_BIO(10) + &
      3.36e-02*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      d_BIO(iC4) = r_BIO(8)

      p_BIO(iC5) = &
      0.00e+00*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      0.00e+00*r_BIO(3) + &
      0.00e+00*r_BIO(4) + &
      0.00e+00*r_BIO(5) + &
      0.00e+00*r_BIO(6) + &
      4.68e-01*r_BIO(7) + &
      4.37e-01*r_BIO(8) + &
      0.00e+00*r_BIO(9) + &
      1.91e-01*r_BIO(10) + &
      2.00e-01*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      d_BIO(iC5) = r_BIO(9)

      p_BIO(iC6) = &
      0.00e+00*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      0.00e+00*r_BIO(3) + &
      0.00e+00*r_BIO(4) + &
      0.00e+00*r_BIO(5) + &
      0.00e+00*r_BIO(6) + &
      0.00e+00*r_BIO(7) + &
      4.37e-01*r_BIO(8) + &
      3.74e-01*r_BIO(9) + &
      0.00e+00*r_BIO(10) + &
      3.81e-01*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      d_BIO(iC6) = r_BIO(10)

      p_BIO(iC7) = &
      0.00e+00*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      0.00e+00*r_BIO(3) + &
      0.00e+00*r_BIO(4) + &
      0.00e+00*r_BIO(5) + &
      0.00e+00*r_BIO(6) + &
      0.00e+00*r_BIO(7) + &
      0.00e+00*r_BIO(8) + &
      3.74e-01*r_BIO(9) + &
      2.49e-01*r_BIO(10) + &
      0.00e+00*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      d_BIO(iC7) = r_BIO(11)

      p_BIO(iC8) = &
      0.00e+00*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      0.00e+00*r_BIO(3) + &
      0.00e+00*r_BIO(4) + &
      0.00e+00*r_BIO(5) + &
      0.00e+00*r_BIO(6) + &
      0.00e+00*r_BIO(7) + &
      0.00e+00*r_BIO(8) + &
      0.00e+00*r_BIO(9) + &
      0.00e+00*r_BIO(10) + &
      0.00e+00*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      d_BIO(iC8) = r_BIO(12)

      p_BIO(iC9) = &
      0.00e+00*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      0.00e+00*r_BIO(3) + &
      0.00e+00*r_BIO(4) + &
      0.00e+00*r_BIO(5) + &
      0.00e+00*r_BIO(6) + &
      0.00e+00*r_BIO(7) + &
      0.00e+00*r_BIO(8) + &
      0.00e+00*r_BIO(9) + &
      0.00e+00*r_BIO(10) + &
      0.00e+00*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      d_BIO(iC9) = r_BIO(13)

      p_OXY(1) = &
      3.32e-03*r_BIO(1) + &
      6.61e-03*r_BIO(2) + &
      1.14e-02*r_BIO(3) + &
      2.66e-02*r_BIO(4) + &
      3.72e-02*r_BIO(5) + &
      2.73e-02*r_BIO(6) + &
      5.68e-03*r_BIO(7) + &
      1.90e-04*r_BIO(8) + &
      8.82e-07*r_BIO(9) + &
      0.00e+00*r_BIO(10) + &
      0.00e+00*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      OXY_AVE = j_GASMOLE_OXY(1)/j_GASMOLE_SOM(1)

      IF (OXY_AVE.GT.10.0) OXY_AVE = 1.0

      IF (j_GASMOLE_SOM(1).LT.1.) OXY_AVE = 1.0

      d_OXY(1) = r_BIO(1)*OXY_AVE

      p_OXY(2) = &
      4.99e-01*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      1.54e-03*r_BIO(3) + &
      1.80e-03*r_BIO(4) + &
      1.41e-03*r_BIO(5) + &
      2.43e-02*r_BIO(6) + &
      4.31e-02*r_BIO(7) + &
      1.09e-02*r_BIO(8) + &
      3.76e-04*r_BIO(9) + &
      1.76e-06*r_BIO(10) + &
      0.00e+00*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      OXY_AVE = j_GASMOLE_OXY(2)/j_GASMOLE_SOM(2)

      IF (OXY_AVE.GT.10.0) OXY_AVE = 1.0

      IF (j_GASMOLE_SOM(2).LT.1.) OXY_AVE = 1.0

      d_OXY(2) = r_BIO(2)*OXY_AVE

      p_OXY(3) = &
      4.99e-01*r_BIO(1) + &
      4.99e-01*r_BIO(2) + &
      0.00e+00*r_BIO(3) + &
      3.07e-03*r_BIO(4) + &
      3.58e-03*r_BIO(5) + &
      2.80e-03*r_BIO(6) + &
      4.84e-02*r_BIO(7) + &
      8.59e-02*r_BIO(8) + &
      2.17e-02*r_BIO(9) + &
      7.49e-04*r_BIO(10) + &
      3.49e-06*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      OXY_AVE = j_GASMOLE_OXY(3)/j_GASMOLE_SOM(3)

      IF (OXY_AVE.GT.10.0) OXY_AVE = 1.0

      IF (j_GASMOLE_SOM(3).LT.1.) OXY_AVE = 1.0

      d_OXY(3) = r_BIO(3)*OXY_AVE

      p_OXY(4) = &
      0.00e+00*r_BIO(1) + &
      4.99e-01*r_BIO(2) + &
      4.98e-01*r_BIO(3) + &
      0.00e+00*r_BIO(4) + &
      6.12e-03*r_BIO(5) + &
      7.12e-03*r_BIO(6) + &
      5.58e-03*r_BIO(7) + &
      9.62e-02*r_BIO(8) + &
      1.71e-01*r_BIO(9) + &
      4.33e-02*r_BIO(10) + &
      1.49e-03*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      OXY_AVE = j_GASMOLE_OXY(4)/j_GASMOLE_SOM(4)

      IF (OXY_AVE.GT.10.0) OXY_AVE = 1.0

      IF (j_GASMOLE_SOM(4).LT.1.) OXY_AVE = 1.0

      d_OXY(4) = r_BIO(4)*OXY_AVE

      p_OXY(5) = &
      0.00e+00*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      4.98e-01*r_BIO(3) + &
      4.96e-01*r_BIO(4) + &
      0.00e+00*r_BIO(5) + &
      1.22e-02*r_BIO(6) + &
      1.42e-02*r_BIO(7) + &
      1.11e-02*r_BIO(8) + &
      1.92e-01*r_BIO(9) + &
      3.40e-01*r_BIO(10) + &
      8.62e-02*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      OXY_AVE = j_GASMOLE_OXY(5)/j_GASMOLE_SOM(5)

      IF (OXY_AVE.GT.10.0) OXY_AVE = 1.0

      IF (j_GASMOLE_SOM(5).LT.1.) OXY_AVE = 1.0

      d_OXY(5) = r_BIO(5)*OXY_AVE

      p_OXY(6) = &
      0.00e+00*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      0.00e+00*r_BIO(3) + &
      4.96e-01*r_BIO(4) + &
      4.92e-01*r_BIO(5) + &
      0.00e+00*r_BIO(6) + &
      2.42e-02*r_BIO(7) + &
      2.82e-02*r_BIO(8) + &
      2.21e-02*r_BIO(9) + &
      3.81e-01*r_BIO(10) + &
      6.77e-01*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      OXY_AVE = j_GASMOLE_OXY(6)/j_GASMOLE_SOM(6)

      IF (OXY_AVE.GT.10.0) OXY_AVE = 1.0

      IF (j_GASMOLE_SOM(6).LT.1.) OXY_AVE = 1.0

      d_OXY(6) = r_BIO(6)*OXY_AVE

      p_OXY(7) = &
      0.00e+00*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      0.00e+00*r_BIO(3) + &
      0.00e+00*r_BIO(4) + &
      4.92e-01*r_BIO(5) + &
      4.84e-01*r_BIO(6) + &
      0.00e+00*r_BIO(7) + &
      4.83e-02*r_BIO(8) + &
      5.62e-02*r_BIO(9) + &
      4.40e-02*r_BIO(10) + &
      7.59e-01*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      OXY_AVE = j_GASMOLE_OXY(7)/j_GASMOLE_SOM(7)

      IF (OXY_AVE.GT.10.0) OXY_AVE = 1.0

      IF (j_GASMOLE_SOM(7).LT.1.) OXY_AVE = 1.0

      d_OXY(7) = r_BIO(7)*OXY_AVE

      p_OXY(8) = &
      0.00e+00*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      0.00e+00*r_BIO(3) + &
      0.00e+00*r_BIO(4) + &
      0.00e+00*r_BIO(5) + &
      4.84e-01*r_BIO(6) + &
      4.68e-01*r_BIO(7) + &
      0.00e+00*r_BIO(8) + &
      9.61e-02*r_BIO(9) + &
      1.12e-01*r_BIO(10) + &
      8.75e-02*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      OXY_AVE = j_GASMOLE_OXY(8)/j_GASMOLE_SOM(8)

      IF (OXY_AVE.GT.10.0) OXY_AVE = 1.0

      IF (j_GASMOLE_SOM(8).LT.1.) OXY_AVE = 1.0

      d_OXY(8) = r_BIO(8)*OXY_AVE

      p_OXY(9) = &
      0.00e+00*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      0.00e+00*r_BIO(3) + &
      0.00e+00*r_BIO(4) + &
      0.00e+00*r_BIO(5) + &
      0.00e+00*r_BIO(6) + &
      4.68e-01*r_BIO(7) + &
      4.37e-01*r_BIO(8) + &
      0.00e+00*r_BIO(9) + &
      1.91e-01*r_BIO(10) + &
      2.23e-01*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      OXY_AVE = j_GASMOLE_OXY(9)/j_GASMOLE_SOM(9)

      IF (OXY_AVE.GT.10.0) OXY_AVE = 1.0

      IF (j_GASMOLE_SOM(9).LT.1.) OXY_AVE = 1.0

      d_OXY(9) = r_BIO(9)*OXY_AVE

      p_OXY(10) = &
      0.00e+00*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      0.00e+00*r_BIO(3) + &
      0.00e+00*r_BIO(4) + &
      0.00e+00*r_BIO(5) + &
      0.00e+00*r_BIO(6) + &
      0.00e+00*r_BIO(7) + &
      4.37e-01*r_BIO(8) + &
      3.74e-01*r_BIO(9) + &
      0.00e+00*r_BIO(10) + &
      3.81e-01*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      OXY_AVE = j_GASMOLE_OXY(10)/j_GASMOLE_SOM(10)

      IF (OXY_AVE.GT.10.0) OXY_AVE = 1.0

      IF (j_GASMOLE_SOM(10).LT.1.) OXY_AVE = 1.0

      d_OXY(10) = r_BIO(10)*OXY_AVE

      p_OXY(11) = &
      0.00e+00*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      0.00e+00*r_BIO(3) + &
      0.00e+00*r_BIO(4) + &
      0.00e+00*r_BIO(5) + &
      0.00e+00*r_BIO(6) + &
      0.00e+00*r_BIO(7) + &
      0.00e+00*r_BIO(8) + &
      3.74e-01*r_BIO(9) + &
      2.49e-01*r_BIO(10) + &
      0.00e+00*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      OXY_AVE = j_GASMOLE_OXY(11)/j_GASMOLE_SOM(11)

      IF (OXY_AVE.GT.10.0) OXY_AVE = 1.0

      IF (j_GASMOLE_SOM(11).LT.1.) OXY_AVE = 1.0

      d_OXY(11) = r_BIO(11)*OXY_AVE

      p_OXY(12) = &
      0.00e+00*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      0.00e+00*r_BIO(3) + &
      0.00e+00*r_BIO(4) + &
      0.00e+00*r_BIO(5) + &
      0.00e+00*r_BIO(6) + &
      0.00e+00*r_BIO(7) + &
      0.00e+00*r_BIO(8) + &
      0.00e+00*r_BIO(9) + &
      0.00e+00*r_BIO(10) + &
      0.00e+00*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      OXY_AVE = j_GASMOLE_OXY(12)/j_GASMOLE_SOM(12)

      IF (OXY_AVE.GT.10.0) OXY_AVE = 1.0

      IF (j_GASMOLE_SOM(12).LT.1.) OXY_AVE = 1.0

      d_OXY(12) = r_BIO(12)*OXY_AVE

      p_OXY(13) = &
      0.00e+00*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      0.00e+00*r_BIO(3) + &
      0.00e+00*r_BIO(4) + &
      0.00e+00*r_BIO(5) + &
      0.00e+00*r_BIO(6) + &
      0.00e+00*r_BIO(7) + &
      0.00e+00*r_BIO(8) + &
      0.00e+00*r_BIO(9) + &
      0.00e+00*r_BIO(10) + &
      0.00e+00*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      OXY_AVE = j_GASMOLE_OXY(13)/j_GASMOLE_SOM(13)

      IF (OXY_AVE.GT.10.0) OXY_AVE = 1.0

      IF (j_GASMOLE_SOM(13).LT.1.) OXY_AVE = 1.0

      d_OXY(13) = r_BIO(13)*OXY_AVE

      

      p_BIO(iOH)= 0.0
      d_BIO(iOH)= 0.0

      p_BIO(iISOP) = 0.0
      d_BIO(iISOP) = r_BIO(14)

      RETURN
      END SUBROUTINE ODE_BIO






