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
      6.61e-03*r_BIO(1) + &
      1.09e-02*r_BIO(2) + &
      1.80e-02*r_BIO(3) + &
      2.73e-02*r_BIO(4) + &
      2.60e-02*r_BIO(5) + &
      1.82e-02*r_BIO(6) + &
      1.05e-02*r_BIO(7) + &
      7.66e-03*r_BIO(8) + &
      1.25e-02*r_BIO(9) + &
      2.06e-02*r_BIO(10) + &
      3.40e-02*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      d_BIO(iCN3) = r_BIO(1)

      p_BIO(iCN2) = &
      5.48e-03*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      7.86e-06*r_BIO(3) + &
      2.51e-03*r_BIO(4) + &
      1.91e-02*r_BIO(5) + &
      2.47e-02*r_BIO(6) + &
      1.96e-02*r_BIO(7) + &
      9.65e-03*r_BIO(8) + &
      1.64e-04*r_BIO(9) + &
      3.83e-05*r_BIO(10) + &
      0.00e+00*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      d_BIO(iCN2) = r_BIO(2)

      p_BIO(iCN1) = &
      5.48e-03*r_BIO(1) + &
      5.46e-03*r_BIO(2) + &
      0.00e+00*r_BIO(3) + &
      1.30e-05*r_BIO(4) + &
      4.15e-03*r_BIO(5) + &
      3.15e-02*r_BIO(6) + &
      4.08e-02*r_BIO(7) + &
      3.24e-02*r_BIO(8) + &
      1.59e-02*r_BIO(9) + &
      2.72e-04*r_BIO(10) + &
      6.32e-05*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      d_BIO(iCN1) = r_BIO(3)

      p_BIO(iC0) = &
      0.00e+00*r_BIO(1) + &
      5.46e-03*r_BIO(2) + &
      5.42e-03*r_BIO(3) + &
      0.00e+00*r_BIO(4) + &
      2.14e-05*r_BIO(5) + &
      6.86e-03*r_BIO(6) + &
      5.20e-02*r_BIO(7) + &
      6.74e-02*r_BIO(8) + &
      5.35e-02*r_BIO(9) + &
      2.63e-02*r_BIO(10) + &
      4.49e-04*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      d_BIO(iC0) = r_BIO(4)

      p_BIO(iC1) = &
      0.00e+00*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      5.42e-03*r_BIO(3) + &
      5.35e-03*r_BIO(4) + &
      0.00e+00*r_BIO(5) + &
      3.54e-05*r_BIO(6) + &
      1.13e-02*r_BIO(7) + &
      8.59e-02*r_BIO(8) + &
      1.11e-01*r_BIO(9) + &
      8.84e-02*r_BIO(10) + &
      4.35e-02*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      d_BIO(iC1) = r_BIO(5)

      p_BIO(iC2) = &
      0.00e+00*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      0.00e+00*r_BIO(3) + &
      5.35e-03*r_BIO(4) + &
      5.25e-03*r_BIO(5) + &
      0.00e+00*r_BIO(6) + &
      5.85e-05*r_BIO(7) + &
      1.87e-02*r_BIO(8) + &
      1.42e-01*r_BIO(9) + &
      1.84e-01*r_BIO(10) + &
      1.46e-01*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      d_BIO(iC2) = r_BIO(6)

      p_BIO(iC3) = &
      0.00e+00*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      0.00e+00*r_BIO(3) + &
      0.00e+00*r_BIO(4) + &
      5.25e-03*r_BIO(5) + &
      5.07e-03*r_BIO(6) + &
      0.00e+00*r_BIO(7) + &
      9.66e-05*r_BIO(8) + &
      3.09e-02*r_BIO(9) + &
      2.34e-01*r_BIO(10) + &
      3.04e-01*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      d_BIO(iC3) = r_BIO(7)

      p_BIO(iC4) = &
      0.00e+00*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      0.00e+00*r_BIO(3) + &
      0.00e+00*r_BIO(4) + &
      0.00e+00*r_BIO(5) + &
      5.07e-03*r_BIO(6) + &
      4.78e-03*r_BIO(7) + &
      0.00e+00*r_BIO(8) + &
      1.60e-04*r_BIO(9) + &
      5.11e-02*r_BIO(10) + &
      3.87e-01*r_BIO(11) + &
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
      4.78e-03*r_BIO(7) + &
      4.29e-03*r_BIO(8) + &
      0.00e+00*r_BIO(9) + &
      2.64e-04*r_BIO(10) + &
      8.44e-02*r_BIO(11) + &
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
      4.29e-03*r_BIO(8) + &
      3.50e-03*r_BIO(9) + &
      0.00e+00*r_BIO(10) + &
      4.35e-04*r_BIO(11) + &
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
      3.50e-03*r_BIO(9) + &
      2.18e-03*r_BIO(10) + &
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
      4.04e-02*r_BIO(1) + &
      6.68e-02*r_BIO(2) + &
      9.25e-02*r_BIO(3) + &
      9.99e-02*r_BIO(4) + &
      9.00e-02*r_BIO(5) + &
      6.74e-02*r_BIO(6) + &
      3.20e-02*r_BIO(7) + &
      3.44e-02*r_BIO(8) + &
      4.99e-02*r_BIO(9) + &
      0.00e+00*r_BIO(10) + &
      0.00e+00*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      OXY_AVE = j_GASMOLE_OXY(1)/j_GASMOLE_SOM(1)

      IF (OXY_AVE.GT.10.0) OXY_AVE = 1.0

      IF (j_GASMOLE_SOM(1).LT.1.) OXY_AVE = 1.0

      d_OXY(1) = r_BIO(1)*OXY_AVE

      p_OXY(2) = &
      5.48e-03*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      7.86e-06*r_BIO(3) + &
      5.01e-03*r_BIO(4) + &
      3.83e-02*r_BIO(5) + &
      5.49e-02*r_BIO(6) + &
      5.88e-02*r_BIO(7) + &
      2.91e-02*r_BIO(8) + &
      6.58e-04*r_BIO(9) + &
      1.53e-04*r_BIO(10) + &
      0.00e+00*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      OXY_AVE = j_GASMOLE_OXY(2)/j_GASMOLE_SOM(2)

      IF (OXY_AVE.GT.10.0) OXY_AVE = 1.0

      IF (j_GASMOLE_SOM(2).LT.1.) OXY_AVE = 1.0

      d_OXY(2) = r_BIO(2)*OXY_AVE

      p_OXY(3) = &
      5.48e-03*r_BIO(1) + &
      5.46e-03*r_BIO(2) + &
      0.00e+00*r_BIO(3) + &
      1.30e-05*r_BIO(4) + &
      8.28e-03*r_BIO(5) + &
      6.33e-02*r_BIO(6) + &
      9.06e-02*r_BIO(7) + &
      9.72e-02*r_BIO(8) + &
      4.80e-02*r_BIO(9) + &
      1.09e-03*r_BIO(10) + &
      2.53e-04*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      OXY_AVE = j_GASMOLE_OXY(3)/j_GASMOLE_SOM(3)

      IF (OXY_AVE.GT.10.0) OXY_AVE = 1.0

      IF (j_GASMOLE_SOM(3).LT.1.) OXY_AVE = 1.0

      d_OXY(3) = r_BIO(3)*OXY_AVE

      p_OXY(4) = &
      0.00e+00*r_BIO(1) + &
      5.46e-03*r_BIO(2) + &
      5.42e-03*r_BIO(3) + &
      0.00e+00*r_BIO(4) + &
      2.14e-05*r_BIO(5) + &
      1.37e-02*r_BIO(6) + &
      1.05e-01*r_BIO(7) + &
      1.50e-01*r_BIO(8) + &
      1.61e-01*r_BIO(9) + &
      7.93e-02*r_BIO(10) + &
      1.79e-03*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      OXY_AVE = j_GASMOLE_OXY(4)/j_GASMOLE_SOM(4)

      IF (OXY_AVE.GT.10.0) OXY_AVE = 1.0

      IF (j_GASMOLE_SOM(4).LT.1.) OXY_AVE = 1.0

      d_OXY(4) = r_BIO(4)*OXY_AVE

      p_OXY(5) = &
      0.00e+00*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      5.42e-03*r_BIO(3) + &
      5.35e-03*r_BIO(4) + &
      0.00e+00*r_BIO(5) + &
      3.54e-05*r_BIO(6) + &
      2.26e-02*r_BIO(7) + &
      1.73e-01*r_BIO(8) + &
      2.47e-01*r_BIO(9) + &
      2.65e-01*r_BIO(10) + &
      1.31e-01*r_BIO(11) + &
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
      5.35e-03*r_BIO(4) + &
      5.25e-03*r_BIO(5) + &
      0.00e+00*r_BIO(6) + &
      5.85e-05*r_BIO(7) + &
      3.73e-02*r_BIO(8) + &
      2.85e-01*r_BIO(9) + &
      4.09e-01*r_BIO(10) + &
      4.38e-01*r_BIO(11) + &
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
      5.25e-03*r_BIO(5) + &
      5.07e-03*r_BIO(6) + &
      0.00e+00*r_BIO(7) + &
      9.66e-05*r_BIO(8) + &
      6.16e-02*r_BIO(9) + &
      4.71e-01*r_BIO(10) + &
      6.75e-01*r_BIO(11) + &
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
      5.07e-03*r_BIO(6) + &
      4.78e-03*r_BIO(7) + &
      0.00e+00*r_BIO(8) + &
      1.60e-04*r_BIO(9) + &
      1.02e-01*r_BIO(10) + &
      7.79e-01*r_BIO(11) + &
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
      4.78e-03*r_BIO(7) + &
      4.29e-03*r_BIO(8) + &
      0.00e+00*r_BIO(9) + &
      2.64e-04*r_BIO(10) + &
      1.68e-01*r_BIO(11) + &
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
      4.29e-03*r_BIO(8) + &
      3.50e-03*r_BIO(9) + &
      0.00e+00*r_BIO(10) + &
      4.35e-04*r_BIO(11) + &
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
      3.50e-03*r_BIO(9) + &
      2.18e-03*r_BIO(10) + &
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






