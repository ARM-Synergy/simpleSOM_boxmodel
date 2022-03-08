!     ==========================================================================
!                 THIS SUBROUTINE DEFINES ODEs FOR BIO-REGIME CHEMISTRY
!     ==========================================================================

      SUBROUTINE ODE_BIO
  
      USE mod_MAIN
      USE mod_GAS

      IMPLICIT NONE

!FLAG1
      p_BIO(iCN3) = &
      9.24e-03*r_BIO(1) + &
      1.36e-02*r_BIO(2) + &
      2.02e-02*r_BIO(3) + &
      2.73e-02*r_BIO(4) + &
      2.33e-02*r_BIO(5) + &
      1.46e-02*r_BIO(6) + &
      7.50e-03*r_BIO(7) + &
      4.89e-03*r_BIO(8) + &
      7.14e-03*r_BIO(9) + &
      1.05e-02*r_BIO(10) + &
      1.56e-02*r_BIO(11) + &
      2.30e-02*r_BIO(12) + &
      3.40e-02*r_BIO(13)

      d_BIO(iCN3) = r_BIO(1)

      p_BIO(iCN2) = &
      5.22e-03*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      7.40e-06*r_BIO(3) + &
      2.51e-03*r_BIO(4) + &
      1.70e-02*r_BIO(5) + &
      1.98e-02*r_BIO(6) + &
      1.40e-02*r_BIO(7) + &
      6.19e-03*r_BIO(8) + &
      7.95e-05*r_BIO(9) + &
      1.66e-05*r_BIO(10) + &
      0.00e+00*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      d_BIO(iCN2) = r_BIO(2)

      p_BIO(iCN1) = &
      5.22e-03*r_BIO(1) + &
      5.20e-03*r_BIO(2) + &
      0.00e+00*r_BIO(3) + &
      1.09e-05*r_BIO(4) + &
      3.70e-03*r_BIO(5) + &
      2.52e-02*r_BIO(6) + &
      2.92e-02*r_BIO(7) + &
      2.08e-02*r_BIO(8) + &
      9.14e-03*r_BIO(9) + &
      1.17e-04*r_BIO(10) + &
      2.45e-05*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      d_BIO(iCN1) = r_BIO(3)

      p_BIO(iC0) = &
      0.00e+00*r_BIO(1) + &
      5.20e-03*r_BIO(2) + &
      5.17e-03*r_BIO(3) + &
      0.00e+00*r_BIO(4) + &
      1.62e-05*r_BIO(5) + &
      5.47e-03*r_BIO(6) + &
      3.72e-02*r_BIO(7) + &
      4.32e-02*r_BIO(8) + &
      3.07e-02*r_BIO(9) + &
      1.35e-02*r_BIO(10) + &
      1.74e-04*r_BIO(11) + &
      3.62e-05*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      d_BIO(iC0) = r_BIO(4)

      p_BIO(iC1) = &
      0.00e+00*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      5.17e-03*r_BIO(3) + &
      5.11e-03*r_BIO(4) + &
      0.00e+00*r_BIO(5) + &
      2.39e-05*r_BIO(6) + &
      8.09e-03*r_BIO(7) + &
      5.50e-02*r_BIO(8) + &
      6.38e-02*r_BIO(9) + &
      4.53e-02*r_BIO(10) + &
      2.00e-02*r_BIO(11) + &
      2.56e-04*r_BIO(12) + &
      5.35e-05*r_BIO(13)

      d_BIO(iC1) = r_BIO(5)

      p_BIO(iC2) = &
      0.00e+00*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      0.00e+00*r_BIO(3) + &
      5.11e-03*r_BIO(4) + &
      5.04e-03*r_BIO(5) + &
      0.00e+00*r_BIO(6) + &
      3.53e-05*r_BIO(7) + &
      1.19e-02*r_BIO(8) + &
      8.12e-02*r_BIO(9) + &
      9.43e-02*r_BIO(10) + &
      6.69e-02*r_BIO(11) + &
      2.95e-02*r_BIO(12) + &
      3.79e-04*r_BIO(13)

      d_BIO(iC2) = r_BIO(6)

      p_BIO(iC3) = &
      0.00e+00*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      0.00e+00*r_BIO(3) + &
      0.00e+00*r_BIO(4) + &
      5.04e-03*r_BIO(5) + &
      4.93e-03*r_BIO(6) + &
      0.00e+00*r_BIO(7) + &
      5.21e-05*r_BIO(8) + &
      1.77e-02*r_BIO(9) + &
      1.20e-01*r_BIO(10) + &
      1.39e-01*r_BIO(11) + &
      9.89e-02*r_BIO(12) + &
      4.36e-02*r_BIO(13)

      d_BIO(iC3) = r_BIO(7)

      p_BIO(iC4) = &
      0.00e+00*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      0.00e+00*r_BIO(3) + &
      0.00e+00*r_BIO(4) + &
      0.00e+00*r_BIO(5) + &
      4.93e-03*r_BIO(6) + &
      4.77e-03*r_BIO(7) + &
      0.00e+00*r_BIO(8) + &
      7.70e-05*r_BIO(9) + &
      2.61e-02*r_BIO(10) + &
      1.77e-01*r_BIO(11) + &
      2.06e-01*r_BIO(12) + &
      1.46e-01*r_BIO(13)

      d_BIO(iC4) = r_BIO(8)

      p_BIO(iC5) = &
      0.00e+00*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      0.00e+00*r_BIO(3) + &
      0.00e+00*r_BIO(4) + &
      0.00e+00*r_BIO(5) + &
      0.00e+00*r_BIO(6) + &
      4.77e-03*r_BIO(7) + &
      4.52e-03*r_BIO(8) + &
      0.00e+00*r_BIO(9) + &
      1.14e-04*r_BIO(10) + &
      3.85e-02*r_BIO(11) + &
      2.62e-01*r_BIO(12) + &
      3.04e-01*r_BIO(13)

      d_BIO(iC5) = r_BIO(9)

      p_BIO(iC6) = &
      0.00e+00*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      0.00e+00*r_BIO(3) + &
      0.00e+00*r_BIO(4) + &
      0.00e+00*r_BIO(5) + &
      0.00e+00*r_BIO(6) + &
      0.00e+00*r_BIO(7) + &
      4.52e-03*r_BIO(8) + &
      4.17e-03*r_BIO(9) + &
      0.00e+00*r_BIO(10) + &
      1.68e-04*r_BIO(11) + &
      5.70e-02*r_BIO(12) + &
      3.87e-01*r_BIO(13)

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
      4.17e-03*r_BIO(9) + &
      3.64e-03*r_BIO(10) + &
      0.00e+00*r_BIO(11) + &
      2.48e-04*r_BIO(12) + &
      8.41e-02*r_BIO(13)

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
      3.64e-03*r_BIO(10) + &
      2.86e-03*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      3.67e-04*r_BIO(13)

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
      2.86e-03*r_BIO(11) + &
      1.70e-03*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      d_BIO(iC9) = r_BIO(13)

      p_OXY(1) = &
      5.65e-02*r_BIO(1) + &
      8.35e-02*r_BIO(2) + &
      1.03e-01*r_BIO(3) + &
      9.99e-02*r_BIO(4) + &
      8.05e-02*r_BIO(5) + &
      5.40e-02*r_BIO(6) + &
      2.28e-02*r_BIO(7) + &
      2.20e-02*r_BIO(8) + &
      2.86e-02*r_BIO(9) + &
      0.00e+00*r_BIO(10) + &
      0.00e+00*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      d_OXY(1) = r_BIO(1)*j_GASMOLE_OXY(1)/ &
      MAX(j_GASMOLE_SOM(1),j_GASMOLE_OXY(1)+0.001)

      p_OXY(2) = &
      5.22e-03*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      7.40e-06*r_BIO(3) + &
      5.00e-03*r_BIO(4) + &
      3.43e-02*r_BIO(5) + &
      4.39e-02*r_BIO(6) + &
      4.21e-02*r_BIO(7) + &
      1.86e-02*r_BIO(8) + &
      3.18e-04*r_BIO(9) + &
      6.63e-05*r_BIO(10) + &
      0.00e+00*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      d_OXY(2) = r_BIO(2)*j_GASMOLE_OXY(2)/ &
      MAX(j_GASMOLE_SOM(2),j_GASMOLE_OXY(2)+0.001)

      p_OXY(3) = &
      5.22e-03*r_BIO(1) + &
      5.20e-03*r_BIO(2) + &
      0.00e+00*r_BIO(3) + &
      1.09e-05*r_BIO(4) + &
      7.39e-03*r_BIO(5) + &
      5.06e-02*r_BIO(6) + &
      6.49e-02*r_BIO(7) + &
      6.23e-02*r_BIO(8) + &
      2.75e-02*r_BIO(9) + &
      4.70e-04*r_BIO(10) + &
      9.80e-05*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      d_OXY(3) = r_BIO(3)*j_GASMOLE_OXY(3)/ &
      MAX(j_GASMOLE_SOM(3),j_GASMOLE_OXY(3)+0.001)

      p_OXY(4) = &
      0.00e+00*r_BIO(1) + &
      5.20e-03*r_BIO(2) + &
      5.17e-03*r_BIO(3) + &
      0.00e+00*r_BIO(4) + &
      1.62e-05*r_BIO(5) + &
      1.09e-02*r_BIO(6) + &
      7.48e-02*r_BIO(7) + &
      9.59e-02*r_BIO(8) + &
      9.20e-02*r_BIO(9) + &
      4.06e-02*r_BIO(10) + &
      6.94e-04*r_BIO(11) + &
      1.45e-04*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      d_OXY(4) = r_BIO(4)*j_GASMOLE_OXY(4)/ &
      MAX(j_GASMOLE_SOM(4),j_GASMOLE_OXY(4)+0.001)

      p_OXY(5) = &
      0.00e+00*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      5.17e-03*r_BIO(3) + &
      5.11e-03*r_BIO(4) + &
      0.00e+00*r_BIO(5) + &
      2.39e-05*r_BIO(6) + &
      1.61e-02*r_BIO(7) + &
      1.11e-01*r_BIO(8) + &
      1.42e-01*r_BIO(9) + &
      1.36e-01*r_BIO(10) + &
      6.00e-02*r_BIO(11) + &
      1.03e-03*r_BIO(12) + &
      2.14e-04*r_BIO(13)

      d_OXY(5) = r_BIO(5)*j_GASMOLE_OXY(5)/ &
      MAX(j_GASMOLE_SOM(5),j_GASMOLE_OXY(5)+0.001)

      p_OXY(6) = &
      0.00e+00*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      0.00e+00*r_BIO(3) + &
      5.11e-03*r_BIO(4) + &
      5.04e-03*r_BIO(5) + &
      0.00e+00*r_BIO(6) + &
      3.53e-05*r_BIO(7) + &
      2.38e-02*r_BIO(8) + &
      1.63e-01*r_BIO(9) + &
      2.09e-01*r_BIO(10) + &
      2.01e-01*r_BIO(11) + &
      8.87e-02*r_BIO(12) + &
      1.52e-03*r_BIO(13)

      d_OXY(6) = r_BIO(6)*j_GASMOLE_OXY(6)/ &
      MAX(j_GASMOLE_SOM(6),j_GASMOLE_OXY(6)+0.001)

      p_OXY(7) = &
      0.00e+00*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      0.00e+00*r_BIO(3) + &
      0.00e+00*r_BIO(4) + &
      5.04e-03*r_BIO(5) + &
      4.93e-03*r_BIO(6) + &
      0.00e+00*r_BIO(7) + &
      5.21e-05*r_BIO(8) + &
      3.52e-02*r_BIO(9) + &
      2.41e-01*r_BIO(10) + &
      3.09e-01*r_BIO(11) + &
      2.97e-01*r_BIO(12) + &
      1.31e-01*r_BIO(13)

      d_OXY(7) = r_BIO(7)*j_GASMOLE_OXY(7)/ &
      MAX(j_GASMOLE_SOM(7),j_GASMOLE_OXY(7)+0.001)

      p_OXY(8) = &
      0.00e+00*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      0.00e+00*r_BIO(3) + &
      0.00e+00*r_BIO(4) + &
      0.00e+00*r_BIO(5) + &
      4.93e-03*r_BIO(6) + &
      4.77e-03*r_BIO(7) + &
      0.00e+00*r_BIO(8) + &
      7.70e-05*r_BIO(9) + &
      5.20e-02*r_BIO(10) + &
      3.57e-01*r_BIO(11) + &
      4.57e-01*r_BIO(12) + &
      4.39e-01*r_BIO(13)

      d_OXY(8) = r_BIO(8)*j_GASMOLE_OXY(8)/ &
      MAX(j_GASMOLE_SOM(8),j_GASMOLE_OXY(8)+0.001)

      p_OXY(9) = &
      0.00e+00*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      0.00e+00*r_BIO(3) + &
      0.00e+00*r_BIO(4) + &
      0.00e+00*r_BIO(5) + &
      0.00e+00*r_BIO(6) + &
      4.77e-03*r_BIO(7) + &
      4.52e-03*r_BIO(8) + &
      0.00e+00*r_BIO(9) + &
      1.14e-04*r_BIO(10) + &
      7.69e-02*r_BIO(11) + &
      5.27e-01*r_BIO(12) + &
      6.75e-01*r_BIO(13)

      d_OXY(9) = r_BIO(9)*j_GASMOLE_OXY(9)/ &
      MAX(j_GASMOLE_SOM(9),j_GASMOLE_OXY(9)+0.001)

      p_OXY(10) = &
      0.00e+00*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      0.00e+00*r_BIO(3) + &
      0.00e+00*r_BIO(4) + &
      0.00e+00*r_BIO(5) + &
      0.00e+00*r_BIO(6) + &
      0.00e+00*r_BIO(7) + &
      4.52e-03*r_BIO(8) + &
      4.17e-03*r_BIO(9) + &
      0.00e+00*r_BIO(10) + &
      1.68e-04*r_BIO(11) + &
      1.14e-01*r_BIO(12) + &
      7.78e-01*r_BIO(13)

      d_OXY(10) = r_BIO(10)*j_GASMOLE_OXY(10)/ &
      MAX(j_GASMOLE_SOM(10),j_GASMOLE_OXY(10)+0.001)

      p_OXY(11) = &
      0.00e+00*r_BIO(1) + &
      0.00e+00*r_BIO(2) + &
      0.00e+00*r_BIO(3) + &
      0.00e+00*r_BIO(4) + &
      0.00e+00*r_BIO(5) + &
      0.00e+00*r_BIO(6) + &
      0.00e+00*r_BIO(7) + &
      0.00e+00*r_BIO(8) + &
      4.17e-03*r_BIO(9) + &
      3.64e-03*r_BIO(10) + &
      0.00e+00*r_BIO(11) + &
      2.48e-04*r_BIO(12) + &
      1.68e-01*r_BIO(13)

      d_OXY(11) = r_BIO(11)*j_GASMOLE_OXY(11)/ &
      MAX(j_GASMOLE_SOM(11),j_GASMOLE_OXY(11)+0.001)

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
      3.64e-03*r_BIO(10) + &
      2.86e-03*r_BIO(11) + &
      0.00e+00*r_BIO(12) + &
      3.67e-04*r_BIO(13)

      d_OXY(12) = r_BIO(12)*j_GASMOLE_OXY(12)/ &
      MAX(j_GASMOLE_SOM(12),j_GASMOLE_OXY(12)+0.001)

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
      2.86e-03*r_BIO(11) + &
      1.70e-03*r_BIO(12) + &
      0.00e+00*r_BIO(13)

      d_OXY(13) = r_BIO(13)*j_GASMOLE_OXY(13)/ &
      MAX(j_GASMOLE_SOM(13),j_GASMOLE_OXY(13)+0.001)

      

      p_BIO(iOH)= 0.0
      d_BIO(iOH)= 0.0

      p_BIO(iISOP) = 0.0
      d_BIO(iISOP) = r_BIO(14)

      RETURN
      END SUBROUTINE ODE_BIO






